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
  if (length(rest) && grepl("^[0-9]{2}-[A-Za-z]{3}-[0-9]{4}$", utils::tail(rest, 1L))) {
    ans$date <- utils::tail(rest, 1L)
    rest <- utils::head(rest, -1L)
  }
  if (length(rest) && grepl("^[A-Za-z]{3}$", utils::tail(rest, 1L))) {
    ans$division <- utils::tail(rest, 1L)
    rest <- utils::head(rest, -1L)
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

# Qualifiers that conventionally have unquoted values in GenBank files.
.gbk_unquoted_qualifiers <- c(
  "codon_start", "transl_table", "direction", "label", "number",
  "estimated_length"
)

.gbk_default_quoted <- function(name, value) {
  !is.null(value) && !(name %in% .gbk_unquoted_qualifiers)
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

