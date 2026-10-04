#' Edit a GenBank Sequence and Remap Its Features
#'
#' Insert, delete, or replace DNA while updating feature coordinates. New
#' annotations use coordinates relative to the replacement sequence.
#'
#' @param x A `gbk_record`, `gbk_file`, non-empty list of records, GenBank path,
#'   or readable connection. The input is not modified.
#' @param start One-based first base to replace. For insertion, the new DNA is
#'   inserted immediately before this position; use the sequence length plus
#'   one to append.
#' @param end Inclusive last base to replace. `NULL` means insertion only.
#'   An explicit interval must be ascending, even for circular records.
#' @param replacement One IUPAC DNA string, or `""` for deletion (default).
#' @param features Optional feature data frame accepted by [features_from_df()]
#'   or parsed feature list. Locations and supported coordinate qualifiers are
#'   relative to `replacement`, starting at one. Appended after retained features.
#' @param overlap Policy for partially affected existing features: `"error"`
#'   (default), `"drop"` to remove the whole feature, or `"trim"` to retain only
#'   surviving original bases. Fully deleted features are always removed.
#' @param record One-based record index to edit when `x` contains multiple
#'   records. Defaults to one. Other records are preserved.
#'
#' @details
#' Features after the edit shift automatically. An insertion affects a range
#' only when it falls inside that range, not immediately before or after it.
#' Under `"trim"`, inserted DNA is excluded from old features, using `join()`
#' where necessary. Joined and complemented locations retain segment order;
#' an edit in the gap between joined segments does not affect their bases.
#' A full-length `source` feature is resized automatically.
#'
#' Supports local positions, ranges, fuzzy bounds, between-base sites,
#' `join()`, `order()`, `complement()`, and circular descending ranges.
#' Circular ranges are expanded into joined ascending segments. Fuzzy markers
#' on retained endpoints are preserved; cut endpoints become exact positions.
#' Between-base sites are removed if an anchor is deleted; insertion directly
#' at a site follows the overlap policy and removes that site under `"trim"`.
#' Remote-accession and unsupported locations are rejected.
#'
#' Coordinates in `transl_except`, `anticodon`, and `rpt_unit_range` qualifiers
#' are updated. A qualifier touched by the edit causes an error under
#' `"error"`; otherwise it is removed with a warning. When a feature is trimmed,
#' `translation`, `protein_id`, and `codon_start` are removed with a warning:
#' coding annotations must be reviewed before reuse. No translation or reading
#' frame is inferred. Newly supplied coding annotations are the caller's
#' responsibility.
#'
#' Sequence length and locus date are updated. Other metadata, including
#' accession/version, headers, free-text notes, and unsupported coordinate-like
#' qualifiers, is retained verbatim and is not biologically revalidated.
#' Empty output sequences and non-DNA records are rejected. Use the coordinates
#' of the current returned object when applying several edits in succession.
#'
#' @return An updated object with the same container class as `x`. A path or
#'   connection returns a `gbk_file`. Use [write_gbk()] to save the result.
#' @seealso [read_gbk()], [write_gbk()], [features_from_df()], [wide_to_gbk()]
#' @export
#' @examples
#' x <- wide_to_gbk(data.frame(subject = strsplit("AACCGGTT", "")[[1]],
#'                             position = 1:8), seq_id = "construct")
#' # Delete bases 3 through 4.
#' x <- splice_gbk(x, start = 3, end = 4)
#' # Insert before base 3, with a feature relative to the inserted DNA.
#' new_features <- data.frame(type = "misc_feature", location = "1..4",
#'                            label = "insert")
#' x <- splice_gbk(x, start = 3, replacement = "ATGC", features = new_features)
#' features_to_df(x)
splice_gbk <- function(x, start, end = NULL, replacement = "",
                       features = NULL, overlap = c("error", "drop", "trim"),
                       record = 1L) {
  overlap <- match.arg(overlap)
  if (is.character(x) || inherits(x, "connection")) x <- read_gbk(x)
  single <- inherits(x, "gbk_record")
  records <- if (single) list(x) else x
  if (!is.list(records) || !length(records) ||
      !all(vapply(records, inherits, logical(1L), "gbk_record"))) {
    stop("x must contain gbk_record objects.", call. = FALSE)
  }
  integer_scalar <- function(value) {
    is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value) && value == floor(value)
  }
  if (!integer_scalar(record) || record < 1 || record > length(records)) {
    stop("record must index one record in x.", call. = FALSE)
  }
  current <- records[[record]]
  sequence <- current$sequence
  valid_dna <- function(s) {
    is.character(s) && length(s) == 1L && !is.na(s) &&
      !grepl("[^ACGTRYSWKMBDHVNacgtryswkmbdhvn]", s)
  }
  if (!valid_dna(sequence) || !nzchar(sequence) || !valid_dna(replacement)) {
    stop("Sequence and replacement must be IUPAC DNA strings; the original sequence cannot be empty.",
         call. = FALSE)
  }
  if (identical(tolower(current$locus$unit %||% "bp"), "aa") ||
      grepl("RNA|protein", current$locus$molecule_type %||% "", ignore.case = TRUE)) {
    stop("Splicing requires a DNA record.", call. = FALSE)
  }
  size <- nchar(sequence)
  inserted <- nchar(replacement)
  if (!integer_scalar(start) || start < 1 || start > size + 1) {
    stop("start must be a position from 1 through sequence length plus one.", call. = FALSE)
  }
  if (is.null(end)) {
    end <- start - 1
  } else if (!integer_scalar(end) || end < start || end > size) {
    stop("end must be at or after start and within the sequence.", call. = FALSE)
  }
  removed <- end - start + 1
  new_size <- size - removed + inserted
  if (new_size < 1) stop("The edited sequence cannot be empty.", call. = FALSE)
  circular <- identical(tolower(current$locus$topology %||% ""), "circular")
  edit <- list(start = start, end = end, inserted = inserted,
               delta = inserted - removed, size = size, circular = circular)
  new_sequence <- paste0(substr(sequence, 1, start - 1), replacement,
                          if (end < size) substring(sequence, end + 1) else "")
  new_features <- if (is.null(features)) list() else if (is.data.frame(features)) {
    features_from_df(features)
  } else features
  if (!is.list(new_features) ||
      (length(new_features) && inserted == 0)) {
    stop("features must be a feature list or data frame and require replacement DNA.", call. = FALSE)
  }
  # Validate and offset new annotations before changing any existing features.
  new_features <- lapply(new_features, function(feature) {
    if (!is.list(feature) || !is.character(feature$type) ||
        length(feature$type) != 1L || is.na(feature$type) || !nzchar(feature$type)) {
      stop("Each new feature must have one non-empty type.", call. = FALSE)
    }
    offset <- list(start = 1, end = 0, inserted = start - 1,
                   delta = start - 1, size = inserted, circular = FALSE)
    feature$location <- .gbk_splice_location(feature$location, offset)$location
    feature$qualifiers <- .gbk_splice_qualifiers(feature$qualifiers, offset, "error")$qualifiers
    feature
  })
  retained <- list()
  invalidated <- integer()
  for (i in seq_along(current$features)) {
    feature <- current$features[[i]]
    mapped <- tryCatch(.gbk_splice_location(feature$location, edit),
      error = function(e) stop("Feature ", i, ": ", conditionMessage(e), call. = FALSE))
    full_source <- identical(feature$type, "source") &&
      identical(gsub("[[:space:]]", "", feature$location), paste0("1..", size))
    if (full_source) {
      mapped <- list(location = paste0("1..", new_size), touched = FALSE)
    } else {
      if (is.null(mapped$location) && !isTRUE(mapped$lost_site)) next
      if (mapped$touched) {
        if (overlap == "error") {
          stop("Feature ", i, " (", feature$type, " at ", feature$location,
               ") overlaps the edit; choose overlap = 'drop' or 'trim'.", call. = FALSE)
        }
        if (overlap == "drop" || is.null(mapped$location)) next
        invalidated <- c(invalidated, i)
        feature$qualifiers <- Filter(function(q) {
          !q$name %in% c("translation", "protein_id", "codon_start")
        }, feature$qualifiers)
      }
    }
    feature$location <- mapped$location
    qualifiers <- .gbk_splice_qualifiers(feature$qualifiers, edit, overlap)
    if (qualifiers$changed) invalidated <- c(invalidated, i)
    feature$qualifiers <- qualifiers$qualifiers
    retained[[length(retained) + 1L]] <- feature
  }
  current$sequence <- new_sequence
  current$features <- c(retained, new_features)
  current$locus$length <- nchar(new_sequence)
  current$locus$date <- toupper(format(Sys.Date(), "%d-%b-%Y"))
  records[[record]] <- current
  if (length(invalidated)) {
    warning("Review edited features ", paste(unique(invalidated), collapse = ", "),
            ": trimmed annotations and/or affected qualifiers; coding qualifiers on trimmed features were removed.",
            call. = FALSE)
  }
  if (single) records[[1L]] else records
}

.gbk_splice_location <- function(location, edit) {
  # Reuse the existing strict syntax/bounds validator without reversing.
  .gbk_reverse_location(location, edit$size, edit$circular, flip = FALSE)
  location <- gsub("[[:space:]]+", "", location)
  map_position <- function(p) ifelse(p > edit$end, p + edit$delta, p)
  walk <- function(loc) {
    for (op in c("complement", "join", "order")) {
      opener <- paste0(op, "(")
      if (startsWith(loc, opener)) {
        parts <- .gbk_split_location(substr(loc, nchar(opener) + 1L, nchar(loc) - 1L))
        mapped <- lapply(parts, walk)
        kept <- Filter(Negate(is.null), lapply(mapped, `[[`, "location"))
        result <- if (!length(kept)) NULL else if (length(kept) == 1L && op != "complement") {
          kept[[1L]]
        } else paste0(opener, paste(unlist(kept), collapse = ","), ")")
        return(list(location = result, touched = any(vapply(mapped, `[[`, logical(1), "touched")),
                    lost_site = any(vapply(mapped, function(m) isTRUE(m$lost_site), logical(1)))))
      }
    }
    separator <- if (grepl("..", loc, fixed = TRUE)) ".." else
      if (grepl("^", loc, fixed = TRUE)) "^" else ""
    bounds <- if (nzchar(separator)) strsplit(loc, separator, fixed = TRUE)[[1]] else loc
    p <- as.numeric(gsub("[<>]", "", bounds))
    markers <- gsub("[0-9]", "", bounds)
    if (separator == "^") {
      deleted <- any(p >= edit$start & p <= edit$end)
      at_site <- edit$inserted > 0 &&
        (edit$start == p[2] || (p[1] == edit$size && p[2] == 1 && edit$start == edit$size + 1))
      if (deleted || at_site) return(list(location = NULL, touched = TRUE,
                                         lost_site = !deleted && at_site))
      return(list(location = paste(map_position(p), collapse = "^"), touched = FALSE))
    }
    a <- p[1]
    b <- p[length(p)]
    if (a > b) {
      return(walk(paste0("join(", bounds[1], "..", edit$size, ",1..", bounds[2], ")")))
    }
    deleted <- edit$start <= b && edit$end >= a && edit$end >= edit$start
    inside <- edit$inserted > 0 && a < edit$start && edit$start <= b
    touched <- deleted || inside
    if (!touched) {
      return(list(location = paste(paste0(markers, map_position(p)), collapse = separator),
                  touched = FALSE))
    }
    chunks <- list()
    if (a < edit$start) chunks[[length(chunks) + 1L]] <- c(a, min(b, edit$start - 1))
    if (b > edit$end) chunks[[length(chunks) + 1L]] <- c(max(a, edit$end + 1), b)
    if (!length(chunks)) return(list(location = NULL, touched = TRUE))
    # Without inserted bases the two surviving pieces are now adjacent.
    if (length(chunks) == 2L && edit$inserted == 0) chunks <- list(c(chunks[[1]][1], chunks[[2]][2]))
    strings <- vapply(chunks, function(chunk) {
      left <- if (chunk[1] == a) markers[1] else ""
      right <- if (chunk[2] == b) markers[length(markers)] else ""
      paste0(left, map_position(chunk[1]), "..", right, map_position(chunk[2]))
    }, character(1))
    list(location = if (length(strings) == 1L) strings else
      paste0("join(", paste(strings, collapse = ","), ")"), touched = TRUE)
  }
  walk(location)
}

.gbk_splice_qualifiers <- function(qualifiers, edit, overlap) {
  changed <- FALSE
  out <- list()
  for (q in qualifiers) {
    if (!is.null(q$value) && q$name %in% c("rpt_unit_range", "transl_except", "anticodon")) {
      location <- q$value
      fields <- NULL
      if (q$name != "rpt_unit_range") {
        value <- trimws(q$value)
        if (!startsWith(value, "(") || !endsWith(value, ")")) {
          stop("Unsupported ", q$name, " qualifier.", call. = FALSE)
        }
        fields <- trimws(.gbk_split_location(substr(value, 2L, nchar(value) - 1L)))
        at <- which(startsWith(fields, "pos:"))
        if (length(at) != 1L) stop("Expected one pos: location in ", q$name, ".", call. = FALSE)
        location <- substring(fields[at], 5L)
      }
      mapped <- .gbk_splice_location(location, edit)
      if (mapped$touched) {
        if (overlap == "error") stop("Qualifier ", q$name, " overlaps the edit.", call. = FALSE)
        changed <- TRUE
        next
      }
      if (is.null(fields)) q$value <- mapped$location else {
        fields[at] <- paste0("pos:", mapped$location)
        q$value <- paste0("(", paste(fields, collapse = ","), ")")
      }
    }
    out[[length(out) + 1L]] <- q
  }
  list(qualifiers = out, changed = changed)
}
