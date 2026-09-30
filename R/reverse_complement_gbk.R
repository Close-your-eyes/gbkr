#' Reverse-Complement GenBank Records and Their Features
#'
#' Reverse-complement each DNA sequence and remap its feature locations so that
#' the annotated biological sequences retain their original orientation.
#'
#' @param x A GenBank file path, readable connection, `gbk_file`, single
#'   `gbk_record`, or non-empty list of `gbk_record` objects.
#'
#' @return A transformed object with the same class as the input object. A file
#'   path or connection returns a `gbk_file`. Use [write_gbk()] to save it.
#'
#' @details
#' Supports IUPAC DNA bases, single positions, ranges, fuzzy bounds marked with
#' `<` or `>`, between-base sites, `join()`, `order()`, `complement()`, and
#' descending ranges on circular records. Compound locations retain their
#' biological segment order. Remote-accession locations, unsupported location
#' syntax, empty sequences, and non-DNA sequences cause an error.
#'
#' Feature order is preserved. The `direction` qualifier swaps LEFT and RIGHT;
#' locations in `transl_except`, `anticodon`, and `rpt_unit_range` qualifiers
#' are remapped. Other qualifiers (including translations), headers, and locus
#' metadata are retained, except that locus length is set from the sequence.
#' Coordinates embedded in free-text notes or headers are not rewritten.
#' The input object and input file are not modified.
#'
#' @seealso [read_gbk()], [write_gbk()], [extract_feature_sequences()]
#' @export
#'
#' @examples
#' record <- structure(
#'   list(
#'     locus = list(name = "example", length = 8L, topology = "linear"),
#'     headers = list(),
#'     features = list(list(type = "gene", location = "2..5",
#'                          qualifiers = list())),
#'     sequence = "aaccttgg"
#'   ),
#'   class = c("gbk_record", "list")
#' )
#' reversed <- reverse_complement_gbk(record)
#' reversed$sequence
#' reversed$features[[1]]$location
#' identical(extract_feature_sequences(record),
#'           extract_feature_sequences(reversed))
#'
#' \dontrun{
#' reversed <- reverse_complement_gbk("input.gbk")
#' write_gbk(reversed, "reverse_complement.gbk")
#' }
reverse_complement_gbk <- function(x) {
  if (is.character(x) || inherits(x, "connection")) x <- read_gbk(x)
  if (inherits(x, "gbk_record")) return(.gbk_reverse_record(x))
  if (!is.list(x) || !length(x) ||
      !all(vapply(x, inherits, logical(1L), "gbk_record"))) {
    stop("x must be a GenBank file path, connection, gbk_record, or non-empty ",
         "list of gbk_record objects.", call. = FALSE)
  }
  for (i in seq_along(x)) {
    x[[i]] <- tryCatch(.gbk_reverse_record(x[[i]]), error = function(e) {
      stop("Could not reverse-complement record ", i, ": ",
           conditionMessage(e), call. = FALSE)
    })
  }
  x
}

.gbk_reverse_record <- function(record) {
  sequence <- record$sequence
  if (!is.character(sequence) || length(sequence) != 1L ||
      is.na(sequence) || !nzchar(sequence) ||
      grepl("[^ACGTRYSWKMBDHVNacgtryswkmbdhvn]", sequence)) {
    stop("Each record must contain a non-empty sequence of IUPAC DNA bases.",
         call. = FALSE)
  }
  if (identical(tolower(record$locus$unit %||% "bp"), "aa") ||
      grepl("RNA|protein", record$locus$molecule_type %||% "",
            ignore.case = TRUE)) {
    stop("Reverse complementation requires a DNA record.", call. = FALSE)
  }
  size <- nchar(sequence)
  circular <- identical(tolower(record$locus$topology %||% ""), "circular")
  for (i in seq_along(record$features)) {
    record$features[[i]] <- tryCatch({
      feature <- record$features[[i]]
      feature$location <- .gbk_reverse_location(feature$location, size, circular)
      feature$qualifiers <- lapply(feature$qualifiers, .gbk_reverse_qualifier,
                                   size = size, circular = circular)
      feature
    }, error = function(e) {
      stop("Could not remap feature ", i, " (",
           record$features[[i]]$location, "): ", conditionMessage(e),
           call. = FALSE)
    })
  }
  record$sequence <- .gbk_reverse_complement(sequence)
  record$locus$length <- size
  record
}

.gbk_reverse_location <- function(location, size, circular, flip = TRUE) {
  if (!is.character(location) || length(location) != 1L || is.na(location) ||
      !nzchar(location)) {
    stop("Feature location must be one non-empty string.", call. = FALSE)
  }
  location <- gsub("[[:space:]]+", "", location)
  mirrored <- .gbk_mirror_location(location, size, circular)
  if (!flip) return(mirrored)
  if (startsWith(mirrored, "complement(")) {
    return(substr(mirrored, 12L, nchar(mirrored) - 1L))
  }
  paste0("complement(", mirrored, ")")
}

# Mirror the coordinates and segment order before changing strand. Keeping the
# tree structure preserves the extraction order of mixed-strand joins as well.
.gbk_mirror_location <- function(location, size, circular) {
  parts <- .gbk_split_location(location)
  if (length(parts) != 1L) {
    stop("Multiple locations require join() or order().", call. = FALSE)
  }
  for (operator in c("complement", "join", "order")) {
    opener <- paste0(operator, "(")
    if (startsWith(location, opener) && endsWith(location, ")")) {
      inner <- substr(location, nchar(opener) + 1L, nchar(location) - 1L)
      parts <- .gbk_split_location(inner)
      if (!length(parts) || (operator == "complement" && length(parts) != 1L)) {
        stop("Invalid ", operator, " location.", call. = FALSE)
      }
      mapped <- vapply(rev(parts), .gbk_mirror_location, character(1L),
                       size = size, circular = circular, USE.NAMES = FALSE)
      return(paste0(opener, paste(mapped, collapse = ","), ")"))
    }
  }
  if (grepl(":", location, fixed = TRUE)) {
    stop("Remote-accession locations are not supported.", call. = FALSE)
  }
  if (!grepl("^[<>]?[0-9]+(\\.\\.[<>]?[0-9]+|\\^[0-9]+)?$", location)) {
    stop("Unsupported GenBank location expression: ", location, call. = FALSE)
  }
  separator <- if (grepl("..", location, fixed = TRUE)) ".." else
    if (grepl("^", location, fixed = TRUE)) "^" else ""
  bounds <- if (nzchar(separator)) strsplit(location, separator, fixed = TRUE)[[1L]] else location
  positions <- as.numeric(gsub("[<>]", "", bounds))
  if (any(!is.finite(positions) | positions < 1 | positions > size)) {
    stop("Feature location is outside the record sequence.", call. = FALSE)
  }
  if (separator == ".." && positions[[1L]] > positions[[2L]] && !circular) {
    stop("A descending range requires a circular record.", call. = FALSE)
  }
  if (separator == "^" &&
      !(positions[[2L]] == positions[[1L]] + 1 ||
        (circular && positions[[1L]] == size && positions[[2L]] == 1))) {
    stop("A between-base site must name adjacent bases.", call. = FALSE)
  }
  markers <- chartr("<>", "><", gsub("[0-9]", "", bounds))
  mapped <- paste0(markers, format(size + 1 - positions, scientific = FALSE,
                                  trim = TRUE))
  paste(rev(mapped), collapse = separator)
}

.gbk_reverse_qualifier <- function(q, size, circular) {
  if (is.null(q$value)) return(q)
  if (identical(q$name, "direction")) {
    replacement <- c(LEFT = "RIGHT", RIGHT = "LEFT", left = "right", right = "left")
    if (q$value %in% names(replacement)) q$value <- unname(replacement[q$value])
  } else if (identical(q$name, "rpt_unit_range")) {
    q$value <- .gbk_reverse_location(q$value, size, circular, flip = FALSE)
  } else if (q$name %in% c("transl_except", "anticodon")) {
    value <- trimws(q$value)
    if (!startsWith(value, "(") || !endsWith(value, ")")) {
      stop("Unsupported ", q$name, " qualifier: ", value, call. = FALSE)
    }
    fields <- trimws(.gbk_split_location(substr(value, 2L, nchar(value) - 1L)))
    at <- which(startsWith(fields, "pos:"))
    if (length(at) != 1L) {
      stop("Expected one pos: location in ", q$name, ".", call. = FALSE)
    }
    fields[at] <- paste0("pos:", .gbk_reverse_location(
      substring(fields[at], 5L), size, circular))
    q$value <- paste0("(", paste(fields, collapse = ","), ")")
  }
  q
}
