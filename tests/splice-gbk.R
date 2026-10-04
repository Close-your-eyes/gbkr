library(gbkr)

expect_error <- function(expr, pattern) {
  message <- tryCatch({ force(expr); NULL }, error = conditionMessage)
  stopifnot(!is.null(message), grepl(pattern, message, fixed = TRUE))
}
make_splice_record <- function(locations, circular = FALSE) {
  x <- wide_to_gbk(data.frame(subject = strsplit("AACCGGTTACGT", "")[[1]],
                              position = 1:12), topology = if (circular) "circular" else "linear")[[1]]
  x$features <- c(x$features, features_from_df(data.frame(
    type = rep("misc_feature", length(locations)), location = locations,
    label = paste0("f", seq_along(locations)))))
  x
}
locations <- function(x) features_to_df(x)$location

# Fully deleted features disappear; later and reverse-strand features shift.
x <- make_splice_record(c("1..2", "3..5", "complement(8..10)", "join(1..2,8..10)"))
original <- x
y <- splice_gbk(x, 3, 5)
stopifnot(identical(x, original), inherits(y, "gbk_record"),
          y$sequence == "aagttacgt", y$locus$length == 9L,
          identical(locations(y), c("1..9", "1..2", "complement(5..7)", "join(1..2,5..7)")))

# Partial overlaps default to errors, even with same-length replacement.
x <- make_splice_record(c("2..8", "complement(join(2..4,7..10))"))
expect_error(splice_gbk(x, 4, 5), "overlaps the edit")
expect_error(splice_gbk(x, 4, 5, "AT"), "overlaps the edit")
y <- splice_gbk(x, 4, 5, overlap = "drop")
stopifnot(identical(locations(y), "1..10"))
y <- suppressWarnings(splice_gbk(x, 4, 5, "ATGC", overlap = "trim"))
stopifnot(identical(locations(y),
                    c("1..14", "join(2..3,8..10)", "complement(join(2..3,9..12))")),
          extract_feature_sequence(y, y$features[[2]]) == "acgtt")

# Insertion boundaries, append, prepend, and relative new annotations.
x <- make_splice_record(c("2..4", "6..8"))
new <- data.frame(type = "CDS", location = "complement(1..3)", label = "new",
                  anticodon = "(pos:1..3,aa:Leu,seq:tag)")
y <- splice_gbk(x, 5, replacement = "CTA", features = new)
stopifnot(identical(locations(y), c("1..15", "2..4", "9..11", "complement(5..7)")),
          extract_feature_sequence(y, y$features[[4]]) == "TAG",
          features_to_df(y)$anticodon[4] == "(pos:5..7,aa:Leu,seq:tag)")
stopifnot(identical(locations(splice_gbk(x, 1, replacement = "AA")),
                    c("1..14", "4..6", "8..10")),
          identical(locations(splice_gbk(x, 13, replacement = "AA")),
                    c("1..14", "2..4", "6..8")))
expect_error(splice_gbk(x, 3, replacement = "A"), "overlaps")
y <- suppressWarnings(splice_gbk(x, 3, replacement = "A", overlap = "trim"))
stopifnot(identical(locations(y), c("1..13", "join(2..2,4..5)", "7..9")))

# Coding and coordinate qualifiers cannot silently remain stale.
x <- make_splice_record("2..10")
x$features[[2]]$qualifiers <- features_from_df(data.frame(
  type = "CDS", location = "2..10", translation = "TEST", protein_id = "old",
  codon_start = "2", transl_except = "(pos:complement(7..9),aa:Sec)",
  rpt_unit_range = "7..9", note = "keep me"))[[1]]$qualifiers
y <- suppressWarnings(splice_gbk(x, 3, 4, overlap = "trim"))
q <- features_to_df(y)
stopifnot(!any(c("translation", "protein_id", "codon_start") %in% names(q)),
          q$transl_except[2] == "(pos:complement(5..7),aa:Sec)",
          q$rpt_unit_range[2] == "5..7", q$note[2] == "keep me")
y <- suppressWarnings(splice_gbk(x, 7, 8, overlap = "trim"))
stopifnot(!any(c("transl_except", "rpt_unit_range") %in% names(features_to_df(y))))

# Fuzzy bounds, circular ranges, mixed strand/order, and between-base sites.
x <- make_splice_record(c("<2..>9", "complement(10..3)", "order(1..2,8..9)"), TRUE)
y <- suppressWarnings(splice_gbk(x, 2, 3, overlap = "trim"))
stopifnot(identical(locations(y), c("1..10", "2..>7", "complement(join(8..10,1..1))",
                                  "order(1..1,6..7)")))
x <- make_splice_record(c("4^5", "12^1"), TRUE)
expect_error(splice_gbk(x, 5, replacement = "A"), "overlaps")
stopifnot(identical(locations(splice_gbk(x, 5, replacement = "A", overlap = "trim")),
                    c("1..13", "13^1")))
expect_error(splice_gbk(x, 1, replacement = "A"), "overlaps")
expect_error(splice_gbk(x, 13, replacement = "A"), "overlaps")
stopifnot(identical(locations(splice_gbk(x, 4, 4)), c("1..11", "11^1")))
x <- make_splice_record("complement(4^5)")
expect_error(splice_gbk(x, 5, replacement = "A"), "overlaps")
stopifnot(length(splice_gbk(x, 5, replacement = "A", overlap = "trim")$features) == 1L)

# Trimming preserves exactly the original surviving bases, irrespective of
# strand, edit size, feature boundaries, or replacement length.
set.seed(41)
for (i in 1:120) {
  bounds <- sort(sample(1:12, 2))
  edit <- sort(sample(1:12, 2))
  inserted <- sample(c("", "A", "AGT", "CCGTA"), 1)
  reverse <- sample(c(TRUE, FALSE), 1)
  loc <- paste(bounds, collapse = "..")
  if (reverse) loc <- paste0("complement(", loc, ")")
  x <- make_splice_record(loc)
  if (edit[1] == 1 && edit[2] == 12 && !nzchar(inserted)) next
  y <- suppressWarnings(splice_gbk(x, edit[1], edit[2], inserted, overlap = "trim"))
  old_positions <- seq.int(bounds[1], bounds[2])
  kept <- old_positions[old_positions < edit[1] | old_positions > edit[2]]
  if (!length(kept)) {
    stopifnot(length(y$features) == 1L)
  } else {
    expected <- paste0(strsplit(x$sequence, "")[[1]][kept], collapse = "")
    if (reverse) expected <- chartr("acgt", "tgca", paste0(rev(strsplit(expected, "")[[1]]), collapse = ""))
    stopifnot(extract_feature_sequence(y, y$features[[2]]) == expected)
  }
}

# Multi-record selection, file input, input preservation, and write/read.
x <- make_splice_record("8..10")
container <- structure(list(first = x, second = x), class = c("gbk_file", "list"))
y <- splice_gbk(container, 2, 3, record = 2)
stopifnot(identical(y[[1]], x), identical(names(y), names(container)), inherits(y, "gbk_file"))
path <- tempfile(fileext = ".gbk")
write_gbk(container, path)
before <- readLines(path)
z <- splice_gbk(path, 2, 3, record = 2)
stopifnot(identical(before, readLines(path)), identical(features_to_df(z[[2]]), features_to_df(y[[2]])))
write_gbk(y, path)
z <- read_gbk(path)
stopifnot(identical(z[[2]]$sequence, y[[2]]$sequence),
          identical(features_to_df(z[[2]]), features_to_df(y[[2]])))
unlink(path)

expect_error(splice_gbk(x, 0), "start")
expect_error(splice_gbk(x, 2.5), "start")
expect_error(splice_gbk(x, 4, 3), "end")
expect_error(splice_gbk(x, 1, 12), "empty")
expect_error(splice_gbk(x, 1, replacement = "AX"), "IUPAC")
expect_error(splice_gbk(x, 1, features = new), "require replacement")
expect_error(splice_gbk(x, 1, replacement = "AA", features = new), "outside")
expect_error(splice_gbk(x, 1, record = 2), "record")
expect_error(splice_gbk(make_splice_record("OTHER:1..2"), 3), "Remote-accession")
expect_error(splice_gbk(make_splice_record("one-of(1,2)"), 3), "Unsupported")

