#' Extract the Sequence of One GenBank Feature
#'
#' Extract the bases described by a GenBank location expression. Complemented
#' locations are returned as reverse complements, and joined locations are
#' concatenated in the order given by the expression.
#'
#' @details
#' Supported location syntax includes single positions, ranges, fuzzy bounds
#' marked with `<` or `>`, `join()`, `order()`, `complement()`, and descending
#' ranges that cross the origin of a circular record.
#'
#' @param record A `gbk_record` containing a non-empty `sequence` component and
#'   locus topology information when circular locations are used.
#' @param feature A parsed feature list, a one-row feature data frame containing
#'   a `location` column, or a single GenBank location string.
#'
#' @return A single character string in the orientation specified by the
#'   feature location.
#'
#' @seealso [extract_feature_sequences()], [features_to_df()]
#' @export
#'
#' @examples
#' record <- structure(
#'   list(
#'     locus = list(topology = "linear"),
#'     features = list(),
#'     sequence = "aaccttgg"
#'   ),
#'   class = c("gbk_record", "list")
#' )
#'
#' extract_feature_sequence(record, "2..5")
#' extract_feature_sequence(record, "complement(2..5)")
extract_feature_sequence <- function(record, feature) {
  if (is.null(record$sequence)) {
    stop("record must contain a sequence.", call. = FALSE)
  }
  if (is.data.frame(feature)) {
    if (nrow(feature) != 1L || !"location" %in% names(feature)) {
      stop("A feature data frame must have one row and a location column.",
           call. = FALSE)
    }
    location <- as.character(feature$location[[1L]])
  } else if (is.character(feature) && length(feature) == 1L) {
    location <- feature
  } else if (is.list(feature) && !is.null(feature$location)) {
    location <- as.character(feature$location[[1L]])
  } else {
    stop("feature must provide one GenBank location.", call. = FALSE)
  }

  sequence <- gsub("[^A-Za-z]", "", record$sequence)
  circular <- !is.null(record$locus$topology) &&
    identical(tolower(record$locus$topology), "circular")
  .gbk_extract_location(sequence, location, circular)
}

#' Extract Sequences for Several GenBank Features
#'
#' Supports ranges, fuzzy range bounds (`<` and `>`), `join()`, `order()`,
#' `complement()`, and descending ranges on circular records. Returned names use
#' the first available `label`, `gene`, or `locus_tag`, with a type/index
#' fallback.
#'
#' @param record A `gbk_record` containing a non-empty `sequence` component.
#' @param features A parsed feature list or a feature data frame. Defaults to all
#'   features in `record`.
#' @param name_by A character vector of qualifier names to try, in order, when
#'   naming the returned sequences.
#' @param use_names A logical value indicating whether to attach unique names to
#'   the result.
#'
#' @return A character vector with one sequence per input feature, in the same
#'   order. If `use_names` is `TRUE`, the vector has unique names.
#'
#' @seealso [extract_feature_sequence()], [features_to_df()]
#' @export
#'
#' @examples
#' features <- list(
#'   list(
#'     type = "gene",
#'     location = "1..4",
#'     qualifiers = list(list(name = "gene", value = "alpha", quoted = TRUE))
#'   ),
#'   list(
#'     type = "CDS",
#'     location = "complement(5..8)",
#'     qualifiers = list(list(name = "locus_tag", value = "cds_1", quoted = TRUE))
#'   )
#' )
#' record <- structure(
#'   list(
#'     locus = list(topology = "linear"),
#'     features = features,
#'     sequence = "aaccttgg"
#'   ),
#'   class = c("gbk_record", "list")
#' )
#'
#' extract_feature_sequences(record)
extract_feature_sequences <- function(
    record,
    features = record$features,
    name_by = c("label", "gene", "locus_tag"),
    use_names = TRUE) {
  if (is.data.frame(features)) {
    if (!"location" %in% names(features)) {
      stop("Feature data frame must contain a location column.", call. = FALSE)
    }
    locations <- as.character(features$location)
    types <- if ("type" %in% names(features)) as.character(features$type) else
      rep("feature", nrow(features))
    result_names <- rep("", nrow(features))
    for (qualifier_name in name_by) {
      if (qualifier_name %in% names(features)) {
        candidate <- as.character(features[[qualifier_name]])
        candidate[is.na(candidate)] <- ""
        candidate <- sub("\n.*$", "", candidate)
        take <- !nzchar(result_names) & nzchar(candidate)
        result_names[take] <- candidate[take]
      }
    }
  } else if (is.list(features)) {
    locations <- vapply(features, function(feature) {
      as.character(feature$location[[1L]])
    }, character(1L))
    types <- vapply(features, function(feature) {
      if (is.null(feature$type)) "feature" else as.character(feature$type[[1L]])
    }, character(1L))
    result_names <- vapply(features, function(feature) {
      for (qualifier_name in name_by) {
        matches <- Filter(function(q) {
          identical(q$name, qualifier_name) && !is.null(q$value)
        }, feature$qualifiers)
        if (length(matches)) return(as.character(matches[[1L]]$value))
      }
      ""
    }, character(1L))
  } else {
    stop("features must be a feature list or data frame.", call. = FALSE)
  }

  sequences <- vapply(seq_along(locations), function(i) {
    tryCatch(
      extract_feature_sequence(record, locations[[i]]),
      error = function(e) stop(
        "Could not extract feature ", i, " (", locations[[i]], "): ",
        conditionMessage(e), call. = FALSE
      )
    )
  }, character(1L))

  if (isTRUE(use_names)) {
    fallback <- paste0(types, "_", seq_along(types))
    result_names[!nzchar(result_names)] <- fallback[!nzchar(result_names)]
    names(sequences) <- make.unique(result_names, sep = "_")
  }
  sequences
}
