#' Title
#'
#' @param x
#' @param ...
#'
#' @returns
#' @export
#'
#' @examples
print.gbk_file <- function(x, ...) {
  cat("<gbk_file>", length(x), "record(s)\n")
  for (i in seq_along(x)) {
    cat(sprintf("  [%d] %s: %d bp, %d feature(s)\n", i,
                x[[i]]$locus$name, nchar(x[[i]]$sequence),
                length(x[[i]]$features)))
  }
  invisible(x)
}

#' Title
#'
#' @param x
#' @param ...
#'
#' @returns
#' @export
#'
#' @examples
print.gbk_record <- function(x, ...) {
  cat(sprintf("<gbk_record> %s: %d bp, %d feature(s)\n",
              x$locus$name, nchar(x$sequence), length(x$features)))
  invisible(x)
}
