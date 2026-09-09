#' Convert GTF Annotations to a GenBank Object
#'
#' Combines a DNA sequence with annotations from a GTF data frame and returns
#' an in-memory GenBank file object suitable for writing with [write_gbk()].
#'
#' @details
#' `gtf_df` must contain the standard GTF columns `seqname`, `source`,
#' `feature`, `start`, `end`, `score`, `strand`, `frame`, and `attributes`.
#' It must also contain a `transcript_id` column and the column selected by
#' `label_attr`. Attributes such as `gene_name` must therefore be extracted
#' from the raw GTF `attributes` column before calling this function.
#'
#' Rows are grouped by `transcript_id`. Multiple intervals belonging to the
#' same transcript are represented with a GenBank `join()` location. Features
#' on the negative strand are wrapped in `complement()`. Interval order is
#' inherited from the row order of `gtf_df`.
#'
#' A `source` feature spanning the complete sequence is added automatically.
#' Remaining annotation columns are converted to GenBank qualifiers by
#' [features_from_df()]. The column selected by `label_attr` is also copied to
#' the `label` qualifier.
#'
#' The input must describe exactly one sequence: all values in `seqname` should
#' be identical. This function constructs an object in memory; it does not
#' write a file.
#'
#' @param gtf_df A data frame containing parsed GTF annotations. In addition to
#'   the standard nine GTF columns, it must contain `transcript_id` and the
#'   annotation column selected by `label_attr`.
#' @param sequence A single character string containing the DNA sequence.
#'   Supply the sequence alone, without a FASTA header, numbering, or
#'   whitespace.
#' @param organism A single character string containing the organism name.
#'   This is added to the `source` feature and the GenBank `DEFINITION` header.
#'   The default is `NA`.
#' @param label_attr A single character string naming a column in `gtf_df`
#'   whose values should be copied to the GenBank `label` qualifier. The column
#'   must remain available after the standard GTF metadata columns are removed.
#'   Defaults to `"gene_name"`.
#' @param topology Genome topology. One of `"circular"` or `"linear"`.
#'   Defaults to `"circular"`.
#'
#' @return A `gbk_file` containing one `gbk_record`. The record contains the
#'   sequence, a `LOCUS` definition, basic headers, a full-length `source`
#'   feature, and the features derived from `gtf_df`.
#'
#' @seealso [write_gbk()], [features_from_df()], [read_gbk()]
#' @export
#'
#' @examples
#'\dontrun{
#' gtf <- igsc::read_gtf("/Users/chris/Desktop/Lymphocryptovirus_humangamma4.gtf")$gtf
#' seq <- igsc::read_fasta("/Users/chris/Desktop/Lymphocryptovirus_humangamma4.fna")
#' organism <- "Lymphocryptovirus_humangamma4"
#' gbkout <- gtf_to_gbk(gtf_df = gtf |>
#'                        dplyr::mutate(dplyr::across(c(gene_id, transcript_id), ~gsub(organism, "", .x))) |>
#'                        dplyr::mutate(gene_name = gsub("Lymphocryptovirus-humangamma4-", "", gene_name)),
#'                      sequence = seq,
#'                      organism = organism)
#' write_gbk(gbkout, "/Users/chris/Desktop/Lymphocryptovirus_humangamma4.gbk")
#' }

gtf_to_gbk <- function(gtf_df,
                       sequence,
                       organism = NA,
                       label_attr = "gene_name",
                       topology = c("circular", "linear")) {

  topology <- rlang::arg_match(topology)
  label_attr <- rlang::arg_match(label_attr, values = names(gtf_df))

  sequence_length = nchar(sequence)
  seq_id <- unique(gtf_df$seqname)

  feature_df <- gtf_df_to_feature_df(gtf_df = gtf_df,
                                     sequence_length = sequence_length,
                                     organism = organism) |>
    dplyr::mutate(label = !!rlang::sym(label_attr))


  record <- structure(
    list(
      locus = list(
        name = seq_id,
        length = sequence_length,
        unit = "bp",
        molecule_type = "DNA",
        topology = topology,
        date = toupper(format(Sys.Date(), "%d-%b-%Y"))
      ),
      headers = list(
        make_header(
          "DEFINITION",
          paste(
            organism,
            "sequence with GTF-derived annotations."
          )
        ),
        make_header("ACCESSION", seq_id),
        make_header("VERSION", seq_id)
      ),
      features = features_from_df(feature_df),
      sequence = tolower(sequence)
    ),
    class = c("gbk_record", "list")
  )

  gbk <- structure(
    list(record),
    class = c("gbk_file", "list")
  )

  return(gbk)
}



gtf_df_to_feature_df <- function(gtf_df,
                                 sequence_length,
                                 organism = NA) {

  gtf_df <- gtf_df |>
    dplyr::select(-c(attributes, score, frame, seqname, source)) |>
    dplyr::rename("type" = feature)
  transcripts <- split(gtf_df, gtf_df$transcript_id)
  feature_df2 <- purrr::map_dfr(transcripts, function(x) {
    x$location <- make_transcript_location(x)
    x |> dplyr::select(-c(start, end)) |> dplyr::distinct()
  })
  feature_df1 <- data.frame(
    type = "source",
    location = paste0("1..", sequence_length),
    organism = organism
  )
  feature_df <- dplyr::bind_rows(feature_df1, feature_df2)
  return(feature_df)
}


make_transcript_location <- function(x) {
  strands <- unique(x$strand)

  if (length(strands) != 1L || !strands %in% c("+", "-")) {
    stop("Each transcript must have exactly one '+' or '-' strand.")
  }

  parts <- paste0(x$start, "..", x$end)

  location <- if (length(parts) == 1L) {
    parts
  } else {
    paste0("join(", paste(parts, collapse = ","), ")")
  }

  if (strands == "-") {
    location <- paste0("complement(", location, ")")
  }

  location
}



make_header <- function(key, value) {
  list(
    key = key,
    prefix = sprintf("%-12s", key),
    value = value
  )
}
