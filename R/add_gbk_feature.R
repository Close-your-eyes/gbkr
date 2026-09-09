#' Convert GenBank Features to a Data Frame
#'
#' There is one row per feature. `type`, `location`, and every qualifier are
#' ordinary character columns. Missing qualifiers are `NA`. A flag qualifier
#' such as `/pseudo` is represented by `flag_value`. If a qualifier occurs more
#' than once on one feature, its values are joined with `repeated_sep`.
#'
#' Quotation information is retained as a data-frame attribute, so an unchanged
#' table can be converted back without altering quoted versus unquoted values.
#'
#' @param x A `gbk_record`, a one-record `gbk_file`, or a parsed feature list.
#' @param repeated_sep A non-empty character string used to separate repeated
#'   qualifier values in a single cell. The default newline keeps repeated
#'   values distinguishable while retaining an ordinary character column.
#' @param flag_value A single character string used to represent a qualifier
#'   that has no value, such as `/pseudo`.
#'
#' @return A `gbk_feature_df`, inheriting from `data.frame`, with character
#'   columns only. Attributes retain repeated-value, flag, and quotation
#'   metadata for [features_from_df()].
#'
#' @seealso [features_from_df()], [add_gbk_feature()]
#' @export
#'
#' @examples
#' features <- list(
#'   list(
#'     type = "CDS",
#'     location = "1..6",
#'     qualifiers = list(
#'       list(name = "gene", value = "abc", quoted = TRUE),
#'       list(name = "note", value = "first", quoted = TRUE),
#'       list(name = "note", value = "second", quoted = TRUE),
#'       list(name = "pseudo", value = NULL, quoted = FALSE)
#'     )
#'   )
#' )
#'
#' feature_df <- features_to_df(features)
#' feature_df
features_to_df <- function(x, repeated_sep = "\n", flag_value = "<flag>") {
  if (!is.character(repeated_sep) || length(repeated_sep) != 1L ||
      !nzchar(repeated_sep)) {
    stop("repeated_sep must be one non-empty character string.", call. = FALSE)
  }
  if (!is.character(flag_value) || length(flag_value) != 1L ||
      is.na(flag_value)) {
    stop("flag_value must be one non-missing character string.", call. = FALSE)
  }
  if (inherits(x, "gbk_file")) {
    if (length(x) != 1L) {
      stop("Pass one gbk_record at a time when converting features.",
           call. = FALSE)
    }
    x <- x[[1L]]
  }
  features <- if (inherits(x, "gbk_record")) x$features else x
  if (!is.list(features)) stop("x must contain a feature list.", call. = FALSE)

  if (!length(features)) {
    out <- data.frame(type = character(), location = character(),
                      stringsAsFactors = FALSE)
    attr(out, "gbk_repeated_sep") <- repeated_sep
    attr(out, "gbk_flag_value") <- flag_value
    attr(out, "gbk_quoted") <- list()
    return(structure(out, class = c("gbk_feature_df", "data.frame")))
  }

  get_field <- function(feature, field) {
    value <- feature[[field]]
    if (is.null(value) || length(value) != 1L) {
      stop("Every feature must have one ", field, ".", call. = FALSE)
    }
    as.character(value)
  }
  out <- data.frame(
    type = vapply(features, get_field, character(1L), field = "type"),
    location = vapply(features, get_field, character(1L), field = "location"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  qualifier_names <- unique(unlist(lapply(features, function(feature) {
    vapply(feature$qualifiers, function(q) q$name, character(1L))
  }), use.names = FALSE))

  quoted_metadata <- stats::setNames(vector("list", length(qualifier_names)),
                              qualifier_names)
  for (qualifier_name in qualifier_names) {
    quoted_cells <- vector("list", length(features))
    cells <- vapply(seq_along(features), function(i) {
      feature <- features[[i]]
      matches <- Filter(function(q) identical(q$name, qualifier_name),
                        feature$qualifiers)
      if (!length(matches)) {
        quoted_cells[[i]] <<- logical()
        return(NA_character_)
      }
      values <- vapply(matches, function(q) {
        if (is.null(q$value)) flag_value else as.character(q$value)[[1L]]
      }, character(1L))
      quoted_cells[[i]] <<- vapply(matches, function(q) isTRUE(q$quoted),
                                   logical(1L))
      paste(values, collapse = repeated_sep)
    }, character(1L))
    out[[qualifier_name]] <- cells
    quoted_metadata[[qualifier_name]] <- quoted_cells
  }
  attr(out, "gbk_repeated_sep") <- repeated_sep
  attr(out, "gbk_flag_value") <- flag_value
  attr(out, "gbk_quoted") <- quoted_metadata
  structure(out, class = c("gbk_feature_df", "data.frame"))
}

#' Convert a Feature Data Frame to GenBank Features
#'
#' Besides tables returned by `features_to_df()`, ordinary data frames work too.
#' All qualifier columns are coerced to character. `NA` means absent,
#' `flag_value` means a flag, and `repeated_sep` separates repeated values.
#'
#' @details
#' When `x` was created by [features_to_df()], stored metadata preserves whether
#' each qualifier value was quoted. For an ordinary data frame, conventional
#' GenBank quotation rules are applied automatically.
#'
#' @param x A data frame containing character-coercible `type` and `location`
#'   columns followed by zero or more qualifier columns.
#' @param repeated_sep A non-empty character string separating repeated
#'   qualifier values. `NULL` uses the value recorded by [features_to_df()], or
#'   a newline if no value was recorded.
#' @param flag_value A single character string representing a flag qualifier.
#'   `NULL` uses the value recorded by [features_to_df()], or `"<flag>"` if no
#'   value was recorded.
#'
#' @return A parsed feature list suitable for assignment to
#'   `record$features`.
#'
#' @seealso [features_to_df()], [add_gbk_feature()]
#' @export
#'
#' @examples
#' feature_df <- data.frame(
#'   type = c("gene", "CDS"),
#'   location = c("1..6", "1..6"),
#'   gene = c("abc", "abc"),
#'   pseudo = c(NA, "<flag>"),
#'   stringsAsFactors = FALSE
#' )
#'
#' features <- features_from_df(feature_df)
#' features[[2]]$qualifiers
features_from_df <- function(x, repeated_sep = NULL, flag_value = NULL) {
  if (!is.data.frame(x)) stop("x must be a data frame.", call. = FALSE)
  if (is.null(repeated_sep)) {
    repeated_sep <- attr(x, "gbk_repeated_sep", exact = TRUE)
    if (is.null(repeated_sep)) repeated_sep <- "\n"
  }
  if (is.null(flag_value)) {
    flag_value <- attr(x, "gbk_flag_value", exact = TRUE)
    if (is.null(flag_value)) flag_value <- "<flag>"
  }
  if (length(repeated_sep) != 1L || is.na(repeated_sep) ||
      !nzchar(repeated_sep)) {
    stop("repeated_sep must be one non-empty character string.", call. = FALSE)
  }
  core <- c("type", "location")
  missing <- setdiff(core, names(x))
  if (length(missing)) {
    stop("Feature data frame is missing: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  qualifier_names <- setdiff(names(x), core)
  quoted_metadata <- attr(x, "gbk_quoted", exact = TRUE)

  lapply(seq_len(nrow(x)), function(i) {
    qualifiers <- list()
    for (qualifier_name in qualifier_names) {
      column <- x[[qualifier_name]]
      cell <- column[[i]]
      if (is.null(cell) || !length(cell) || is.na(cell)) next
      values <- strsplit(as.character(cell), repeated_sep, fixed = TRUE)[[1L]]

      quoted <- NULL
      if (!is.null(quoted_metadata) &&
          qualifier_name %in% names(quoted_metadata) &&
          length(quoted_metadata[[qualifier_name]]) >= i) {
        quoted <- quoted_metadata[[qualifier_name]][[i]]
      }
      if (is.null(quoted) || length(quoted) != length(values)) {
        quoted <- vapply(values, function(value) {
          .gbk_default_quoted(
            qualifier_name,
            if (identical(value, flag_value)) NULL else value
          )
        }, logical(1L))
      }

      for (j in seq_along(values)) {
        is_flag <- identical(values[[j]], flag_value)
        qualifiers[[length(qualifiers) + 1L]] <- list(
          name = qualifier_name,
          value = if (is_flag) NULL else values[[j]],
          quoted = if (is_flag) FALSE else isTRUE(quoted[[j]])
        )
      }
    }
    list(
      type = as.character(x$type[[i]]),
      location = as.character(x$location[[i]]),
      qualifiers = qualifiers
    )
  })
}

#' Append a GenBank Feature to a Data Frame
#'
#' Qualifiers are supplied through `...`. Use a vector to add a repeated
#' qualifier, and use `NULL` for a flag qualifier. Override automatic quotation
#' with a named `.quoted` vector, for example `.quoted = c(label = FALSE)`.
#'
#' @param x A feature data frame containing `type` and `location` columns,
#'   typically returned by [features_to_df()].
#' @param type A single feature key, such as `"CDS"`, `"promoter"`, or
#'   `"misc_feature"`.
#' @param location A single GenBank location expression, such as `"10..50"` or
#'   `"complement(10..50)"`.
#' @param ... Named qualifier values. A vector creates a repeated qualifier,
#'   `NULL` creates a flag qualifier, and omitted qualifiers are recorded as
#'   missing.
#' @param .quoted An optional logical value or named logical vector controlling
#'   whether qualifier values are quoted when written. An unnamed scalar applies
#'   to every supplied qualifier; names target individual qualifiers.
#'
#' @return The updated data frame as a `gbk_feature_df`. Existing rows and
#'   quotation metadata are retained.
#'
#' @seealso [features_to_df()], [features_from_df()]
#' @export
#'
#' @examples
#' feature_df <- data.frame(
#'   type = "source",
#'   location = "1..100",
#'   organism = "synthetic construct",
#'   stringsAsFactors = FALSE
#' )
#'
#' feature_df <- add_gbk_feature(
#'   feature_df,
#'   type = "CDS",
#'   location = "10..60",
#'   gene = "exampleA",
#'   note = c("predicted", "reviewed"),
#'   pseudo = NULL
#' )
#' feature_df
#' features_from_df(feature_df)
add_gbk_feature <- function(x, type, location, ..., .quoted = NULL) {
  if (!is.data.frame(x)) stop("x must be a data frame.", call. = FALSE)
  if (!all(c("type", "location") %in% names(x))) {
    stop("x must contain type and location columns.", call. = FALSE)
  }
  if (length(type) != 1L || length(location) != 1L) {
    stop("type and location must each have length one.", call. = FALSE)
  }
  supplied <- list(...)
  if (length(supplied) &&
      (is.null(names(supplied)) || any(!nzchar(names(supplied))))) {
    stop("All qualifiers in ... must be named.", call. = FALSE)
  }

  repeated_sep <- attr(x, "gbk_repeated_sep", exact = TRUE)
  if (is.null(repeated_sep)) repeated_sep <- "\n"
  flag_value <- attr(x, "gbk_flag_value", exact = TRUE)
  if (is.null(flag_value)) flag_value <- "<flag>"
  quoted_metadata <- attr(x, "gbk_quoted", exact = TRUE)
  if (is.null(quoted_metadata)) quoted_metadata <- list()

  core <- c("type", "location")
  existing_qualifiers <- setdiff(names(x), core)
  x$type <- as.character(x$type)
  x$location <- as.character(x$location)
  for (qualifier_name in existing_qualifiers) {
    x[[qualifier_name]] <- as.character(x[[qualifier_name]])
  }

  new_qualifiers <- setdiff(names(supplied), names(x))
  for (qualifier_name in new_qualifiers) {
    x[[qualifier_name]] <- rep(NA_character_, nrow(x))
    quoted_metadata[[qualifier_name]] <- rep(list(logical()), nrow(x))
  }

  new_row <- nrow(x) + 1L
  x[new_row, ] <- NA_character_
  x$type[[new_row]] <- as.character(type)
  x$location[[new_row]] <- as.character(location)
  all_qualifiers <- setdiff(names(x), core)
  for (qualifier_name in all_qualifiers) {
    x[[qualifier_name]][[new_row]] <- NA_character_
    if (is.null(quoted_metadata[[qualifier_name]])) {
      quoted_metadata[[qualifier_name]] <- rep(list(logical()), new_row)
    } else {
      quoted_metadata[[qualifier_name]][[new_row]] <- logical()
    }
  }

  for (qualifier_name in names(supplied)) {
    value <- supplied[[qualifier_name]]
    if (is.null(value)) {
      values <- flag_value
      quote_value <- FALSE
    } else {
      values <- as.character(value)
      quote_value <- vapply(values, function(v) {
        .gbk_default_quoted(qualifier_name, v)
      }, logical(1L))
    }
    if (!is.null(.quoted)) {
      if (is.null(names(.quoted)) && length(.quoted) == 1L) {
        quote_value <- rep_len(as.logical(.quoted), length(values))
      } else if (qualifier_name %in% names(.quoted)) {
        quote_value <- rep_len(as.logical(.quoted[[qualifier_name]]),
                               length(values))
      }
    }
    x[[qualifier_name]][[new_row]] <- paste(values, collapse = repeated_sep)
    quoted_metadata[[qualifier_name]][[new_row]] <- quote_value
  }

  row.names(x) <- seq_len(nrow(x))
  attr(x, "gbk_repeated_sep") <- repeated_sep
  attr(x, "gbk_flag_value") <- flag_value
  attr(x, "gbk_quoted") <- quoted_metadata
  structure(x, class = c("gbk_feature_df", "data.frame"))
}
