#' Print a GenBank File Summary
#'
#' Print a compact summary containing the number of records and, for each
#' record, its locus name, sequence length, and feature count.
#'
#' @param x A `gbk_file` object returned by [read_gbk()].
#' @param ... Additional arguments. Currently ignored.
#'
#' @return `x`, invisibly.
#'
#' @seealso [read_gbk()], [print.gbk_record()]
#' @export
#'
#' @examples
#' record <- structure(
#'   list(
#'     locus = list(name = "example"),
#'     features = list(list(type = "source", location = "1..4")),
#'     sequence = "acgt"
#'   ),
#'   class = c("gbk_record", "list")
#' )
#' gbk <- structure(list(record), class = c("gbk_file", "list"))
#' print(gbk)
print.gbk_file <- function(x, ...) {
  cat("<gbk_file>", length(x), "record(s)\n")
  for (i in seq_along(x)) {
    cat(sprintf("  [%d] %s: %d bp, %d feature(s)\n", i,
                x[[i]]$locus$name, nchar(x[[i]]$sequence),
                length(x[[i]]$features)))
  }
  invisible(x)
}

#' Print a GenBank Record Summary
#'
#' Print a compact summary containing the record's locus name, sequence length,
#' and feature count.
#'
#' @param x A `gbk_record` object, usually an element of a `gbk_file` returned
#'   by [read_gbk()].
#' @param ... Additional arguments. Currently ignored.
#'
#' @return `x`, invisibly.
#'
#' @seealso [read_gbk()], [print.gbk_file()]
#' @export
#'
#' @examples
#' record <- structure(
#'   list(
#'     locus = list(name = "example"),
#'     features = list(list(type = "source", location = "1..4")),
#'     sequence = "acgt"
#'   ),
#'   class = c("gbk_record", "list")
#' )
#' print(record)
print.gbk_record <- function(x, ...) {
  cat(sprintf("<gbk_record> %s: %d bp, %d feature(s)\n",
              x$locus$name, nchar(x$sequence), length(x$features)))
  invisible(x)
}



#' Print a GenBank Feature Data Frame
#'
#' Prints a compact representation of a `gbk_feature_df`. Repeated qualifier
#' values are displayed on one line separated by `" | "`. Qualifier values
#' longer than 60 characters are truncated for display.
#'
#' @param x A `gbk_feature_df` object, typically returned by
#'   [features_to_df()] or [add_gbk_feature()].
#' @param ... Additional arguments passed to [print.data.frame()].
#'
#' @return `x`, invisibly. The original object is not modified.
#'
#' @seealso [features_to_df()], [features_from_df()], [add_gbk_feature()]
#' @export
#'
#' @examples
print.gbk_feature_df <- function(x, ...) {
  display <- x
  class(display) <- "data.frame"
  qualifier_names <- setdiff(names(display), c("type", "location"))
  repeated_sep <- attr(x, "gbk_repeated_sep", exact = TRUE)
  if (is.null(repeated_sep)) repeated_sep <- "\n"
  for (qualifier_name in qualifier_names) {
    display[[qualifier_name]] <- vapply(display[[qualifier_name]], function(cell) {
      if (is.na(cell)) return(NA_character_)
      value <- gsub(repeated_sep, " | ", as.character(cell), fixed = TRUE)
      if (nchar(value) > 60L) paste0(substr(value, 1L, 57L), "...") else value
    }, character(1L))
  }
  print.data.frame(display, ...)
  invisible(x)
}
