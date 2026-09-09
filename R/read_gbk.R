#' Read one or more records from a GenBank flat file.
#'
#' @param file File path or readable connection.
#' @return A `gbk_file` object: a list of `gbk_record` objects.
read_gbk <- function(file) {
  lines <- readLines(file, warn = FALSE)
  lines <- sub("\\r$", "", lines)
  terminators <- which(trimws(lines) == "//")
  if (!length(terminators)) {
    stop("No GenBank record terminator ('//') was found.", call. = FALSE)
  }

  starts <- c(1L, head(terminators, -1L) + 1L)
  records <- Map(function(from, to) {
    chunk <- lines[from:(to - 1L)]
    chunk <- chunk[nzchar(trimws(chunk))]
    .parse_record(chunk)
  }, starts, terminators)

  last_terminator <- tail(terminators, 1L)
  trailing <- if (last_terminator < length(lines))
    lines[(last_terminator + 1L):length(lines)] else character()
  if (length(trailing) && any(nzchar(trimws(trailing)))) {
    warning("Ignoring non-blank text after the final GenBank terminator.", call. = FALSE)
  }
  structure(records, class = c("gbk_file", "list"))
}
