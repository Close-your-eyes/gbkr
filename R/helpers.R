# Dependency-free reading and writing of GenBank flat files.
#
# Public functions:
#   read_gbk(file)
#   write_gbk(x, file, update_length = TRUE, validate = TRUE)
#   features_to_df(x)
#   features_from_df(x)
#   add_gbk_feature(x, type, location, ...)
#   extract_feature_sequences(record, features = record$features)
#
# A parsed file is a "gbk_file": a list of "gbk_record" objects.  Each record
# has four editable components:
#   $locus      list(name, length, unit, molecule_type, topology, division, date)
#   $headers    ordered list of list(key, prefix, value)
#   $features   ordered list of list(type, location, qualifiers)
#   $sequence   one lower-case character string
#
# Qualifiers are stored as an ordered list (rather than a named list) because a
# GenBank feature may legally contain a qualifier more than once.  Each entry is
# list(name, value, quoted); value is NULL for flag qualifiers such as /pseudo.
#
# Example:
#   source("genbank_io.R")
#   gb <- read_gbk("input.gbk")
#   gb[[1]]$locus$name <- "edited_record"
#   feature_table <- features_to_df(gb[[1]])
#   feature_table <- add_gbk_feature(
#     feature_table, "promoter", "100..200",
#     label = "new promoter", note = "added in R"
#   )
#   gb[[1]]$features <- features_from_df(feature_table)
#   feature_sequences <- extract_feature_sequences(gb[[1]])
#   write_gbk(gb, "edited.gbk")

.gbk_substr <- function(x, first, last = nchar(x)) {
  if (nchar(x) < first) return("")
  substr(x, first, last)
}

.parse_locus <- function(line) {
  raw_value <- trimws(.gbk_substr(line, 13L))
  fields <- strsplit(raw_value, "[[:space:]]+")[[1L]]

  ans <- list(
    name = if (length(fields) >= 1L) fields[[1L]] else "",
    length = if (length(fields) >= 2L && grepl("^[0-9]+$", fields[[2L]]))
      as.integer(fields[[2L]]) else NA_integer_,
    unit = if (length(fields) >= 3L) fields[[3L]] else "bp",
    molecule_type = "",
    topology = "",
    division = "",
    date = "",
    raw = line
  )

  rest <- if (length(fields) > 3L) fields[-seq_len(3L)] else character()
  if (length(rest) && grepl("^[0-9]{2}-[A-Za-z]{3}-[0-9]{4}$", tail(rest, 1L))) {
    ans$date <- tail(rest, 1L)
    rest <- head(rest, -1L)
  }
  if (length(rest) && grepl("^[A-Za-z]{3}$", tail(rest, 1L))) {
    ans$division <- tail(rest, 1L)
    rest <- head(rest, -1L)
  }
  topo_at <- which(tolower(rest) %in% c("linear", "circular"))
  if (length(topo_at)) {
    i <- topo_at[[1L]]
    ans$topology <- rest[[i]]
    rest <- rest[-i]
  }
  ans$molecule_type <- paste(rest, collapse = " ")
  ans
}

.parse_header <- function(lines) {
  entries <- list()
  current <- NULL

  flush <- function() {
    if (is.null(current)) return()
    pieces <- trimws(current$pieces)
    current$value <<- paste(pieces[nzchar(pieces)], collapse = " ")
    current$pieces <<- NULL
    entries[[length(entries) + 1L]] <<- current
  }

  for (line in lines) {
    padded <- paste0(line, strrep(" ", max(0L, 12L - nchar(line))))
    prefix <- substr(padded, 1L, 12L)
    key <- trimws(prefix)
    value <- .gbk_substr(line, 13L)
    if (nzchar(key)) {
      flush()
      current <- list(key = key, prefix = prefix, pieces = value)
    } else if (!is.null(current)) {
      current$pieces <- c(current$pieces, value)
    }
  }
  flush()
  entries
}

.parse_qualifier <- function(parts) {
  first <- trimws(parts[[1L]])
  body <- substring(first, 2L)
  eq <- regexpr("=", body, fixed = TRUE)[[1L]]

  if (eq < 0L) {
    return(list(name = trimws(body), value = NULL, quoted = FALSE))
  }

  name <- trimws(substr(body, 1L, eq - 1L))
  value_parts <- c(substr(body, eq + 1L, nchar(body)),
                   if (length(parts) > 1L) trimws(parts[-1L]) else character())
  quoted <- length(value_parts) > 0L && startsWith(value_parts[[1L]], '"')

  if (quoted) {
    value_parts[[1L]] <- substring(value_parts[[1L]], 2L)
    last <- length(value_parts)
    if (endsWith(value_parts[[last]], '"')) {
      value_parts[[last]] <- substr(value_parts[[last]], 1L,
                                    max(0L, nchar(value_parts[[last]]) - 1L))
    }
  }

  if (identical(name, "translation")) {
    value <- paste0(gsub("[[:space:]]+", "", value_parts), collapse = "")
  } else {
    value <- paste(trimws(value_parts), collapse = " ")
  }
  list(name = name, value = value, quoted = quoted)
}

.parse_features <- function(lines) {
  features <- list()
  current <- NULL
  qparts <- NULL

  flush_qualifier <- function() {
    if (is.null(qparts) || is.null(current)) return()
    current$qualifiers[[length(current$qualifiers) + 1L]] <<-
      .parse_qualifier(qparts)
    qparts <<- NULL
  }
  flush_feature <- function() {
    if (is.null(current)) return()
    flush_qualifier()
    current$location <<- paste0(current$location_parts, collapse = "")
    current$location_parts <<- NULL
    features[[length(features) + 1L]] <<- current
    current <<- NULL
  }

  for (line in lines) {
    is_feature <- grepl("^ {5}[^ ]", line) && !grepl("^ {21}", line)
    payload <- trimws(.gbk_substr(line, 22L))

    if (is_feature) {
      flush_feature()
      current <- list(
        type = trimws(.gbk_substr(line, 6L, 20L)),
        location_parts = payload,
        qualifiers = list()
      )
    } else if (!is.null(current) && startsWith(payload, "/")) {
      flush_qualifier()
      qparts <- payload
    } else if (!is.null(qparts)) {
      qparts <- c(qparts, payload)
    } else if (!is.null(current) && nzchar(payload)) {
      current$location_parts <- c(current$location_parts, payload)
    }
  }
  flush_feature()
  features
}

.parse_record <- function(lines) {
  if (!length(lines)) stop("Encountered an empty GenBank record.", call. = FALSE)
  locus_at <- which(startsWith(lines, "LOCUS"))[1L]
  if (is.na(locus_at)) stop("A GenBank record has no LOCUS line.", call. = FALSE)

  features_at <- which(startsWith(lines, "FEATURES"))[1L]
  origin_at <- which(startsWith(lines, "ORIGIN"))[1L]

  header_end <- if (!is.na(features_at)) features_at - 1L else if (!is.na(origin_at)) origin_at - 1L else length(lines)
  header_lines <- if (header_end > locus_at) lines[(locus_at + 1L):header_end] else character()

  feature_lines <- character()
  if (!is.na(features_at)) {
    feature_end <- if (!is.na(origin_at)) origin_at - 1L else length(lines)
    if (feature_end > features_at) feature_lines <- lines[(features_at + 1L):feature_end]
  }

  sequence <- ""
  if (!is.na(origin_at) && origin_at < length(lines)) {
    seq_lines <- lines[(origin_at + 1L):length(lines)]
    sequence <- paste0(gsub("[^A-Za-z]", "", seq_lines), collapse = "")
    sequence <- tolower(sequence)
  }

  structure(list(
    locus = .parse_locus(lines[[locus_at]]),
    headers = .parse_header(header_lines),
    features = .parse_features(feature_lines),
    sequence = sequence
  ), class = c("gbk_record", "list"))
}

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

# More verbose aliases for callers who prefer them.
read_genbank <- read_gbk

# Qualifiers that conventionally have unquoted values in GenBank files.
.gbk_unquoted_qualifiers <- c(
  "codon_start", "transl_table", "direction", "label", "number",
  "estimated_length"
)

.gbk_default_quoted <- function(name, value) {
  !is.null(value) && !(name %in% .gbk_unquoted_qualifiers)
}

#' Convert GenBank features to a wide, character-only data frame.
#'
#' There is one row per feature. `type`, `location`, and every qualifier are
#' ordinary character columns. Missing qualifiers are `NA`. A flag qualifier
#' such as `/pseudo` is represented by `flag_value`. If a qualifier occurs more
#' than once on one feature, its values are joined with `repeated_sep`.
#'
#' Quotation information is retained as a data-frame attribute, so an unchanged
#' table can be converted back without altering quoted versus unquoted values.
#'
#' @param x A `gbk_record`, a one-record `gbk_file`, or a feature list.
#' @param repeated_sep Separator for repeated qualifier values. A newline keeps
#'   the representation unambiguous while remaining an ordinary character cell.
#' @param flag_value Text used for a qualifier that has no value.
#' @return A `gbk_feature_df` containing only character columns.
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

  quoted_metadata <- setNames(vector("list", length(qualifier_names)),
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

#' Convert a feature data frame back to the GenBank feature-list structure.
#'
#' Besides tables returned by `features_to_df()`, ordinary data frames work too.
#' All qualifier columns are coerced to character. `NA` means absent,
#' `flag_value` means a flag, and `repeated_sep` separates repeated values.
#'
#' @param x A data frame containing `type`, `location`, and optional qualifiers.
#' @param repeated_sep Separator used for repeated qualifier values. By default,
#'   use the value recorded by `features_to_df()`, or a newline.
#' @param flag_value Text representing a flag qualifier. By default, use the
#'   value recorded by `features_to_df()`, or `<flag>`.
#' @return A list suitable for assignment to `record$features`.
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

# Clear aliases that describe the direction of conversion.
features_df_to_list <- features_from_df
features_list_to_df <- features_to_df

#' Append one feature to a wide feature data frame.
#'
#' Qualifiers are supplied through `...`. Use a vector to add a repeated
#' qualifier, and use `NULL` for a flag qualifier. Override automatic quotation
#' with a named `.quoted` vector, for example `.quoted = c(label = FALSE)`.
#'
#' @param x A feature data frame.
#' @param type Feature key such as `CDS`, `promoter`, or `misc_feature`.
#' @param location GenBank location expression such as `10..50`.
#' @param ... Named qualifier values.
#' @param .quoted Optional named logical vector controlling quotation.
#' @return The updated `gbk_feature_df`.
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

.gbk_reverse_complement <- function(sequence) {
  complemented <- chartr(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
    sequence
  )
  paste0(rev(strsplit(complemented, "", fixed = TRUE)[[1L]]), collapse = "")
}

.gbk_split_location <- function(location) {
  chars <- strsplit(location, "", fixed = TRUE)[[1L]]
  if (!length(chars)) return(character())
  depth <- 0L
  start <- 1L
  pieces <- character()
  for (i in seq_along(chars)) {
    if (chars[[i]] == "(") {
      depth <- depth + 1L
    } else if (chars[[i]] == ")") {
      depth <- depth - 1L
      if (depth < 0L) stop("Unbalanced location parentheses.", call. = FALSE)
    } else if (chars[[i]] == "," && depth == 0L) {
      pieces <- c(pieces, substr(location, start, i - 1L))
      start <- i + 1L
    }
  }
  if (depth != 0L) stop("Unbalanced location parentheses.", call. = FALSE)
  c(pieces, substr(location, start, nchar(location)))
}

.gbk_extract_location <- function(sequence, location, circular = FALSE) {
  location <- gsub("[[:space:]]+", "", location)
  if (is.na(location) || !nzchar(location)) {
    stop("Feature location is empty.", call. = FALSE)
  }

  for (operator in c("complement", "join", "order")) {
    opener <- paste0(operator, "(")
    if (startsWith(location, opener) && endsWith(location, ")")) {
      inner <- substr(location, nchar(opener) + 1L, nchar(location) - 1L)
      if (operator == "complement") {
        return(.gbk_reverse_complement(
          .gbk_extract_location(sequence, inner, circular)
        ))
      }
      parts <- .gbk_split_location(inner)
      return(paste0(vapply(parts, function(part) {
        .gbk_extract_location(sequence, part, circular)
      }, character(1L)), collapse = ""))
    }
  }

  if (grepl(":", location, fixed = TRUE)) {
    stop("Remote-accession locations cannot be extracted from the local record.",
         call. = FALSE)
  }

  clean <- gsub("[<>]", "", location)
  sequence_length <- nchar(sequence)

  if (grepl("^[0-9]+\\.\\.[0-9]+$", clean)) {
    bounds <- as.integer(strsplit(clean, "..", fixed = TRUE)[[1L]])
    if (any(bounds < 1L) || any(bounds > sequence_length)) {
      stop("Feature location is outside the record sequence.", call. = FALSE)
    }
    if (bounds[[1L]] <= bounds[[2L]]) {
      return(substr(sequence, bounds[[1L]], bounds[[2L]]))
    }
    if (!isTRUE(circular)) {
      stop("A descending range requires a circular record.", call. = FALSE)
    }
    return(paste0(substr(sequence, bounds[[1L]], sequence_length),
                  substr(sequence, 1L, bounds[[2L]])))
  }

  if (grepl("^[0-9]+\\^[0-9]+$", clean)) {
    return("") # A between-base site contains no bases.
  }

  if (grepl("^[0-9]+$", clean)) {
    position <- as.integer(clean)
    if (position < 1L || position > sequence_length) {
      stop("Feature location is outside the record sequence.", call. = FALSE)
    }
    return(substr(sequence, position, position))
  }

  stop("Unsupported GenBank location expression: ", location, call. = FALSE)
}

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

.wrap_words <- function(x, width) {
  if (!nzchar(x)) return("")
  strwrap(x, width = max(1L, width), simplify = TRUE, exdent = 0L,
          initial = "", prefix = "")
}

.format_locus <- function(locus, sequence_length = NA_integer_) {
  len <- locus$length
  if (!is.na(sequence_length)) len <- sequence_length
  if (is.null(len) || is.na(len)) len <- 0L

  tail_fields <- c(locus$molecule_type, locus$topology, locus$division, locus$date)
  tail_fields <- tail_fields[!vapply(tail_fields, is.null, logical(1L))]
  tail_fields <- tail_fields[nzchar(tail_fields)]
  sprintf("LOCUS       %-16s %11d %-5s %s",
          locus$name %||% "", as.integer(len), locus$unit %||% "bp",
          paste(tail_fields, collapse = " "))
}

`%||%` <- function(x, y) if (is.null(x)) y else x

.format_headers <- function(headers) {
  out <- character()
  for (entry in headers) {
    prefix <- entry$prefix %||% sprintf("%-12s", entry$key %||% "")
    prefix <- sprintf("%-12s", substr(prefix, 1L, 12L))
    wrapped <- .wrap_words(entry$value %||% "", 68L)
    if (!length(wrapped)) wrapped <- ""
    out <- c(out, paste0(prefix, wrapped[[1L]]))
    if (length(wrapped) > 1L) {
      out <- c(out, paste0(strrep(" ", 12L), wrapped[-1L]))
    }
  }
  out
}

.split_fixed <- function(x, width) {
  if (!nzchar(x)) return("")
  starts <- seq.int(1L, nchar(x), by = width)
  substring(x, starts, pmin(starts + width - 1L, nchar(x)))
}

.format_location <- function(type, location) {
  pieces <- .split_fixed(gsub("[[:space:]]+", "", location), 59L)
  c(sprintf("     %-15s %s", type, pieces[[1L]]),
    if (length(pieces) > 1L) paste0(strrep(" ", 21L), pieces[-1L]))
}

.format_qualifier <- function(q) {
  indent <- strrep(" ", 21L)
  name <- q$name %||% ""
  if (is.null(q$value)) return(paste0(indent, "/", name))

  quoted <- isTRUE(q$quoted)
  opener <- paste0("/", name, "=", if (quoted) '"' else "")
  closer <- if (quoted) '"' else ""
  value <- as.character(q$value)

  # A conservative common width keeps every line <= 80 characters, including
  # the qualifier name on the first line and the closing quote on the last.
  value_width <- max(1L, 58L - nchar(opener) - nchar(closer))
  if (identical(name, "translation")) {
    parts <- .split_fixed(gsub("[[:space:]]+", "", value), value_width)
  } else {
    parts <- .wrap_words(value, value_width)
  }
  if (!length(parts)) parts <- ""

  first <- paste0(indent, opener, parts[[1L]])
  rest <- if (length(parts) > 1L) paste0(indent, parts[-1L]) else character()
  result <- c(first, rest)
  result[[length(result)]] <- paste0(result[[length(result)]], closer)
  result
}

.format_sequence <- function(sequence) {
  seq <- tolower(gsub("[[:space:][:digit:]]+", "", sequence))
  if (!nzchar(seq)) return(character())
  starts <- seq.int(1L, nchar(seq), by = 60L)
  vapply(starts, function(i) {
    block <- substr(seq, i, min(i + 59L, nchar(seq)))
    groups <- .split_fixed(block, 10L)
    sprintf("%9d %s", i, paste(groups, collapse = " "))
  }, character(1L), USE.NAMES = FALSE)
}

.validate_record <- function(record, update_length) {
  required <- c("locus", "headers", "features", "sequence")
  missing <- setdiff(required, names(record))
  if (length(missing)) {
    stop("Record is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  seq_clean <- gsub("[[:space:][:digit:]]+", "", record$sequence)
  if (nzchar(seq_clean) && grepl("[^A-Za-z]", seq_clean)) {
    stop("Sequence contains characters other than letters.", call. = FALSE)
  }
  if (!update_length && nzchar(seq_clean) && !is.na(record$locus$length) &&
      nchar(seq_clean) != record$locus$length) {
    stop("LOCUS length does not equal sequence length. Use update_length = TRUE ",
         "or correct the record.", call. = FALSE)
  }
  invisible(TRUE)
}

.format_record <- function(record, update_length, validate) {
  if (validate) .validate_record(record, update_length)
  seq_clean <- gsub("[[:space:][:digit:]]+", "", record$sequence)
  locus_length <- if (update_length && nzchar(seq_clean)) nchar(seq_clean) else NA_integer_

  out <- c(.format_locus(record$locus, locus_length),
           .format_headers(record$headers))
  if (length(record$features)) {
    out <- c(out, "FEATURES             Location/Qualifiers")
    for (feature in record$features) {
      out <- c(out, .format_location(feature$type, feature$location))
      if (length(feature$qualifiers)) {
        for (q in feature$qualifiers) out <- c(out, .format_qualifier(q))
      }
    }
  }
  if (nzchar(seq_clean)) out <- c(out, "ORIGIN", .format_sequence(seq_clean))
  c(out, "//")
}

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

write_genbank <- write_gbk

print.gbk_file <- function(x, ...) {
  cat("<gbk_file>", length(x), "record(s)\n")
  for (i in seq_along(x)) {
    cat(sprintf("  [%d] %s: %d bp, %d feature(s)\n", i,
                x[[i]]$locus$name, nchar(x[[i]]$sequence),
                length(x[[i]]$features)))
  }
  invisible(x)
}

print.gbk_record <- function(x, ...) {
  cat(sprintf("<gbk_record> %s: %d bp, %d feature(s)\n",
              x$locus$name, nchar(x$sequence), length(x$features)))
  invisible(x)
}
