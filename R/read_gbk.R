#' Read a GenBank Flat File
#'
#' `read_gbk()` parses one or more GenBank records from a text file or
#' connection. Records must be separated by the standard `//` terminator.
#'
#' @details
#' Each parsed record contains four editable components:
#'
#' * `locus`: the fields parsed from the `LOCUS` line.
#' * `headers`: the remaining header entries, in their original order.
#' * `features`: feature keys, locations, and ordered qualifiers.
#' * `sequence`: the lower-case sequence with numbering and whitespace removed.
#'
#' Repeated qualifiers are retained as separate entries and flag qualifiers,
#' such as `/pseudo`, have a `NULL` value.
#'
#' @param file A character string naming a GenBank flat file, or a readable
#'   [connection][base::connections].
#'
#' @return A `gbk_file` object, which is a list of `gbk_record` objects in file
#'   order.
#'
#' @seealso [write_gbk()], [features_to_df()]
#' @export
#'
#' @examples
#' gbk_lines <- c(
#'   "LOCUS       TEST 12 bp DNA linear PLN 01-JAN-2000",
#'   "DEFINITION  A small example record.",
#'   "FEATURES             Location/Qualifiers",
#'   "     source          1..12",
#'   "                     /organism=\"synthetic construct\"",
#'   "ORIGIN",
#'   "        1 acgtacgtacgt",
#'   "//"
#' )
#' con <- textConnection(gbk_lines)
#' gbk <- read_gbk(con)
#' close(con)
#'
#' gbk
#' gbk[[1]]$sequence
read_gbk <- function(file) {
  lines <- readLines(file, warn = FALSE)
  lines <- sub("\\r$", "", lines)
  terminators <- which(trimws(lines) == "//")
  if (!length(terminators)) {
    stop("No GenBank record terminator ('//') was found.", call. = FALSE)
  }

  starts <- c(1L, utils::head(terminators, -1L) + 1L)
  records <- Map(function(from, to) {
    chunk <- lines[from:(to - 1L)]
    chunk <- chunk[nzchar(trimws(chunk))]
    .parse_record(chunk)
  }, starts, terminators)

  last_terminator <- utils::tail(terminators, 1L)
  trailing <- if (last_terminator < length(lines))
    lines[(last_terminator + 1L):length(lines)] else character()
  if (length(trailing) && any(nzchar(trimws(trailing)))) {
    warning("Ignoring non-blank text after the final GenBank terminator.", call. = FALSE)
  }
  structure(records, class = c("gbk_file", "list"))
}
