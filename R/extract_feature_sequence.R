#' Extract the sequence for one GenBank feature.
#'
#' @param record A `gbk_record` containing the source sequence.
#' @param feature A feature list, a one-row feature data frame, or a location
#'   string.
#' @return One character string in the orientation specified by the location.
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

#' Extract sequences for several GenBank features.
#'
#' Supports ranges, fuzzy range bounds (`<` and `>`), `join()`, `order()`,
#' `complement()`, and descending ranges on circular records. Returned names use
#' the first available `label`, `gene`, or `locus_tag`, with a type/index
#' fallback.
#'
#' @param record A `gbk_record` containing the source sequence.
#' @param features A feature list or feature data frame.
#' @param name_by Qualifiers to try, in order, when naming the result.
#' @param use_names Whether to attach useful, unique names.
#' @return A character vector aligned with `features`.
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
