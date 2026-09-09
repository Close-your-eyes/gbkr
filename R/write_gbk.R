#' Write one or more records to a GenBank flat file.
#'
#' @param x A `gbk_file`, a single `gbk_record`, or a list of record objects.
#' @param file Output path or writable connection.
#' @param update_length If TRUE, update each LOCUS length from its sequence.
#' @param validate If TRUE, check the record structure and sequence.
#' @return The output path or connection, invisibly.
write_gbk <- function(x, file, update_length = TRUE, validate = TRUE) {
  records <- if (inherits(x, "gbk_record")) list(x) else unclass(x)
  if (!length(records)) stop("There are no records to write.", call. = FALSE)

  output <- unlist(lapply(records, .format_record,
                          update_length = update_length, validate = validate),
                   use.names = FALSE)
  writeLines(output, con = file, sep = "\n", useBytes = TRUE)
  invisible(file)
}
