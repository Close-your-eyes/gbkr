#' Convert Per-Base Wide Annotations to a GenBank Object
#'
#' Builds a DNA sequence from one base per row and converts annotation columns
#' to labelled GenBank features. Returns an object like [gtf_to_gbk()] and can
#' optionally write it directly with [write_gbk()].
#'
#' @details
#' Each annotation column represents one feature. Non-missing, non-empty cells
#' mark bases belonging to that feature; cell contents are not used as sequence
#' or qualifiers. This supports columns containing aligned bases or other
#' markers. In particular, zero, `FALSE`, and `"-"` are present markers, not
#' missing values. Entirely empty columns are skipped. Overlapping columns
#' produce overlapping features. Separated runs in a column become a `join()`
#' location; negative-strand features use `complement()` around the location.
#' No strand or coding status is inferred from column names or cell contents.
#'
#' The sequence column must contain one IUPAC DNA base per row. With
#' `position_col`, rows are sorted by position, which must contain each integer
#' from one through the number of rows exactly once. With `position_col = NULL`,
#' row order defines positions. Alignment gaps and missing sequence bases are
#' rejected rather than silently changing feature coordinates.
#'
#' A full-length `source` feature is added automatically. Annotation column
#' names become `label` qualifiers. Use [features_to_df()] and
#' [features_from_df()] to add further qualifiers to the returned record.
#'
#' @param wide_df A non-empty data frame with one row per DNA base, a sequence
#'   column, and one column per annotation.
#' @param seq_id A non-empty record name without whitespace. Defaults to
#'   `"sequence"`; used for LOCUS, ACCESSION, and VERSION.
#' @param organism An organism name, or `NA` to omit the organism qualifier.
#' @param topology One of `"circular"` (default) or `"linear"`.
#' @param sequence_col Name of the base column. Defaults to `"subject"`.
#' @param position_col Name of the position column, or `NULL` for row order.
#'   Defaults to `"position"`.
#' @param annotation_cols Character vector of annotation column names. `NULL`
#'   selects all columns except the sequence, position, and `ignore_cols`.
#'   Use `character()` for a sequence with only a source feature.
#' @param ignore_cols Columns excluded from automatic annotation selection.
#'   Defaults to `"subject.position"`. Names absent from the data are allowed.
#' @param feature_type A feature key applied to all annotations (default
#'   `"misc_feature"`), or a named character vector with exactly one key for
#'   every selected annotation column, for example `c(geneA = "CDS")`.
#' @param strand `"+"` (default) or `"-"` for all annotations, or a named
#'   character vector with exactly one strand for every selected column.
#' @param file Optional output path or writable connection passed to
#'   [write_gbk()]. `NULL` (default) constructs the object without writing.
#'
#' @return A `gbk_file` containing one `gbk_record`, also when `file` is supplied.
#' @seealso [gtf_to_gbk()], [write_gbk()], [features_to_df()]
#' @export
#' @examples
#' wide <- data.frame(
#'   subject = c("A", "T", "G", "C", "A", "T"),
#'   position = 1:6,
#'   geneA = c("A", "T", "G", NA, NA, NA),
#'   site = c(NA, NA, "G", "C", NA, NA)
#' )
#' gbk <- wide_to_gbk(wide, seq_id = "construct",
#'                    feature_type = c(geneA = "CDS", site = "misc_feature"),
#'                    strand = c(geneA = "+", site = "-"))
#' features_to_df(gbk)
#' path <- tempfile(fileext = ".gbk")
#' wide_to_gbk(wide, seq_id = "construct", file = path)
#' unlink(path)
wide_to_gbk <- function(wide_df,
                        seq_id = "sequence",
                        organism = NA,
                        topology = c("circular", "linear"),
                        sequence_col = "subject",
                        position_col = "position",
                        annotation_cols = NULL,
                        ignore_cols = "subject.position",
                        feature_type = "misc_feature",
                        strand = "+",
                        file = NULL) {
  topology <- match.arg(topology)
  if (!is.data.frame(wide_df) || !nrow(wide_df)) {
    stop("wide_df must be a non-empty data frame.", call. = FALSE)
  }
  if (is.null(names(wide_df)) || anyNA(names(wide_df)) ||
      any(!nzchar(names(wide_df))) || anyDuplicated(names(wide_df))) {
    stop("wide_df must have unique, non-empty column names.", call. = FALSE)
  }
  check_column <- function(x, argument) {
    if (!is.character(x) || length(x) != 1L || is.na(x) ||
        !x %in% names(wide_df)) {
      stop(argument, " must name one column in wide_df.", call. = FALSE)
    }
  }
  check_column(sequence_col, "sequence_col")
  if (!is.null(position_col)) {
    check_column(position_col, "position_col")
    if (identical(sequence_col, position_col)) {
      stop("sequence_col and position_col must differ.", call. = FALSE)
    }
    positions <- wide_df[[position_col]]
    if (!is.numeric(positions) || anyNA(positions) ||
        !identical(as.double(sort(positions)), as.double(seq_len(nrow(wide_df))))) {
      stop("position_col must contain each integer from 1 to nrow(wide_df) exactly once.",
           call. = FALSE)
    }
    wide_df <- wide_df[order(positions), , drop = FALSE]
  }
  if (!is.character(seq_id) || length(seq_id) != 1L || is.na(seq_id) ||
      !nzchar(seq_id) || grepl("[[:space:]]", seq_id)) {
    stop("seq_id must be one non-empty string without whitespace.", call. = FALSE)
  }
  if (length(organism) != 1L ||
      !(is.character(organism) || (is.logical(organism) && is.na(organism))) ||
      (!is.na(organism) && !nzchar(trimws(organism)))) {
    stop("organism must be one non-empty string or NA.", call. = FALSE)
  }
  if (!is.character(ignore_cols) || anyNA(ignore_cols)) {
    stop("ignore_cols must be a character vector without NA.", call. = FALSE)
  }
  if (is.null(annotation_cols)) {
    annotation_cols <- setdiff(names(wide_df),
                               c(sequence_col, position_col, ignore_cols))
  }
  if (!is.character(annotation_cols) || anyNA(annotation_cols) ||
      anyDuplicated(annotation_cols) ||
      !all(annotation_cols %in% names(wide_df)) ||
      any(annotation_cols %in% c(sequence_col, position_col))) {
    stop("annotation_cols must uniquely name annotation columns, excluding sequence and position.",
         call. = FALSE)
  }
  bases <- wide_df[[sequence_col]]
  if (!(is.character(bases) || is.factor(bases)) || anyNA(bases) ||
      any(!grepl("^[ACGTRYSWKMBDHVN]$", as.character(bases), ignore.case = TRUE))) {
    stop("sequence_col must contain one non-missing IUPAC DNA base per row.",
         call. = FALSE)
  }
  sequence <- tolower(paste0(as.character(bases), collapse = ""))

  per_annotation <- function(x, argument) {
    if (!is.character(x) || anyNA(x) || any(!nzchar(x))) {
      stop(argument, " must contain non-empty character values.", call. = FALSE)
    }
    if (length(x) == 1L && is.null(names(x))) {
      return(rep(x, length(annotation_cols)))
    }
    if (is.null(names(x)) || anyNA(names(x)) || anyDuplicated(names(x)) ||
        !setequal(names(x), annotation_cols)) {
      stop(argument, " must be one unnamed value or a vector named for every annotation column.",
           call. = FALSE)
    }
    unname(x[annotation_cols])
  }
  # Validate supplied values even if no annotation columns are selected.
  if (any(!grepl("^[A-Za-z][A-Za-z0-9_'*-]{0,14}$", feature_type))) {
    stop("feature_type must contain feature keys of at most 15 characters.", call. = FALSE)
  }
  types <- per_annotation(feature_type, "feature_type")
  if (anyNA(strand) || !all(strand %in% c("+", "-"))) {
    stop("strand must contain only '+' or '-'.", call. = FALSE)
  }
  strands <- per_annotation(strand, "strand")

  features <- features_from_df(data.frame(
    type = "source", location = paste0("1..", nrow(wide_df)),
    organism = organism, stringsAsFactors = FALSE
  ))
  for (i in seq_along(annotation_cols)) {
    column <- wide_df[[annotation_cols[[i]]]]
    if (!is.atomic(column) || !is.null(dim(column))) {
      stop("Annotation column '", annotation_cols[[i]],
           "' must be a one-dimensional atomic vector.", call. = FALSE)
    }
    present <- !is.na(column) & nzchar(trimws(as.character(column)))
    positions <- which(present)
    if (!length(positions)) next
    starts <- positions[c(TRUE, diff(positions) != 1L)]
    ends <- positions[c(diff(positions) != 1L, TRUE)]
    location <- make_transcript_location(data.frame(
      start = starts, end = ends, strand = strands[[i]]
    ))
    features <- c(features, features_from_df(data.frame(
      type = types[[i]], location = location, label = annotation_cols[[i]],
      stringsAsFactors = FALSE
    )))
  }
  definition <- paste(c(if (!is.na(organism)) organism,
                        "sequence with wide-table-derived annotations."),
                      collapse = " ")
  record <- structure(list(
    locus = list(name = seq_id, length = nchar(sequence), unit = "bp",
                 molecule_type = "DNA", topology = topology,
                 date = toupper(format(Sys.Date(), "%d-%b-%Y"))),
    headers = list(make_header("DEFINITION", definition),
                   make_header("ACCESSION", seq_id),
                   make_header("VERSION", seq_id)),
    features = features,
    sequence = sequence
  ), class = c("gbk_record", "list"))
  gbk <- structure(list(record), class = c("gbk_file", "list"))
  if (!is.null(file)) write_gbk(gbk, file)
  gbk
}
