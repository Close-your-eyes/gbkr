#' Write a GenBank Flat File
#'
#' `write_gbk()` serializes one or more parsed or manually constructed records
#' in GenBank flat-file format.
#'
#' @param x A `gbk_file`, a single `gbk_record`, or a non-empty list of
#'   `gbk_record` objects.
#' @param file A character string naming the output file, or a writable
#'   [connection][base::connections].
#' @param update_length A logical value. If `TRUE`, replace each `LOCUS` length
#'   with the number of characters in the corresponding sequence.
#' @param validate A logical value. If `TRUE`, validate the record structure
#'   and sequence before formatting it.
#'
#' @return `file`, invisibly.
#'
#' @seealso [read_gbk()], [features_from_df()]
#' @export
#'
#' @examples
#' gbk_lines <- c(
#'   "LOCUS       TEST 8 bp DNA linear PLN 01-JAN-2000",
#'   "FEATURES             Location/Qualifiers",
#'   "     source          1..8",
#'   "ORIGIN",
#'   "        1 acgtacgt",
#'   "//"
#' )
#' con <- textConnection(gbk_lines)
#' record <- read_gbk(con)[[1]]
#' close(con)
#'
#' path <- tempfile(fileext = ".gbk")
#' write_gbk(record, path)
#' readLines(path)
#' unlink(path)
write_gbk <- function(x, file, update_length = TRUE, validate = TRUE) {
  records <- if (inherits(x, "gbk_record")) list(x) else unclass(x)
  if (!length(records)) stop("There are no records to write.", call. = FALSE)

  output <- unlist(lapply(records, .format_record,
                          update_length = update_length, validate = validate),
                   use.names = FALSE)
  writeLines(output, con = file, sep = "\n", useBytes = TRUE)
  invisible(file)
}
