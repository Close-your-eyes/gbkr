#' Make Consecutive Features from Named DNA Sequences
#'
#' Treat a named character vector as consecutive DNA parts in vector order.
#' Each part becomes one feature labelled with its name. The resulting table
#' can be passed directly to [splice_gbk()] or [features_from_df()].
#'
#' @param sequences A non-empty named character vector of non-empty IUPAC DNA
#'   strings. Names must be unique and non-blank. Upper- and lower-case bases
#'   are accepted. Whitespace and alignment gaps in sequences are rejected.
#' @param feature_type One feature key for all parts (default `"misc_feature"`),
#'   or a character vector named for every part. Feature types are not inferred
#'   from part names.
#' @param strand `"+"` (default) or `"-"` for every part, or a character vector
#'   named for every part. Negative strands wrap the location in `complement()`;
#'   they do not reverse-complement the supplied DNA.
#' @param start One-based start of the first part. Defaults to one. Leave this
#'   at one when supplying the features to [splice_gbk()], which offsets them
#'   to the insertion position automatically.
#'
#' @details
#' Parts are concatenated without separators or gaps. Use
#' `paste0(sequences, collapse = "")` for the corresponding replacement DNA.
#' This function does not search an existing sequence for matches. It assigns
#' coordinates from part lengths and vector order. No translations or other
#' coding qualifiers are inferred.
#'
#' @return A data frame with `type`, `location`, and `label` character columns,
#'   with one row per input part in the same order.
#' @seealso [splice_gbk()], [features_from_df()], [wide_to_gbk()]
#' @export
#' @examples
#' parts <- c(NotI = "gcggccgc", Kozak = "gccacc", insert = "ATGAGC")
#' new_features <- features_from_sequences(parts)
#' new_features
#' replacement <- paste0(parts, collapse = "")
#'
#' gbk <- wide_to_gbk(data.frame(subject = c("A", "C", "G", "T"),
#'                               position = 1:4))
#' gbk <- splice_gbk(gbk, start = 3, replacement = replacement,
#'                   features = new_features)
features_from_sequences <- function(sequences, feature_type = "misc_feature",
                                    strand = "+", start = 1L) {
  if (!is.character(sequences) || !is.null(dim(sequences)) ||
      !length(sequences) || anyNA(sequences) ||
      any(!grepl("^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+$", sequences))) {
    stop("sequences must be a non-empty character vector of non-empty IUPAC DNA strings.",
         call. = FALSE)
  }
  labels <- names(sequences)
  if (is.null(labels) || anyNA(labels) || any(!nzchar(trimws(labels))) ||
      anyDuplicated(labels) || any(grepl("[\r\n]", labels))) {
    stop("sequences must have unique, non-blank names without line breaks.", call. = FALSE)
  }
  if (!is.numeric(start) || length(start) != 1L || is.na(start) ||
      !is.finite(start) || start < 1 || start != floor(start)) {
    stop("start must be one positive whole number.", call. = FALSE)
  }
  expand <- function(value, argument) {
    if (!is.character(value) || !is.null(dim(value)) || anyNA(value) ||
        any(!nzchar(value))) {
      stop(argument, " must contain non-empty character values.", call. = FALSE)
    }
    if (length(value) == 1L && is.null(names(value))) {
      return(rep(value, length(sequences)))
    }
    if (is.null(names(value)) || anyNA(names(value)) ||
        anyDuplicated(names(value)) || !setequal(names(value), labels)) {
      stop(argument, " must be one unnamed value or a vector named for every sequence.",
           call. = FALSE)
    }
    unname(value[labels])
  }
  types <- expand(feature_type, "feature_type")
  strands <- expand(strand, "strand")
  if (any(!grepl("^[A-Za-z][A-Za-z0-9_'*-]{0,14}$", types))) {
    stop("feature_type must contain feature keys of at most 15 characters.", call. = FALSE)
  }
  if (!all(strands %in% c("+", "-"))) {
    stop("strand must contain only '+' or '-'.", call. = FALSE)
  }
  widths <- as.double(nchar(sequences))
  ends <- start - 1 + cumsum(widths)
  if (any(ends > .Machine$integer.max)) {
    stop("Feature coordinates exceed the supported integer range.", call. = FALSE)
  }
  locations <- paste0(format(ends - widths + 1, scientific = FALSE, trim = TRUE),
                       "..", format(ends, scientific = FALSE, trim = TRUE))
  negative <- strands == "-"
  locations[negative] <- paste0("complement(", locations[negative], ")")
  data.frame(type = types, location = locations, label = labels,
             stringsAsFactors = FALSE, row.names = NULL)
}
