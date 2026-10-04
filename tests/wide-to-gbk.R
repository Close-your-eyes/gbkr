library(gbkr)

expect_error <- function(expr, pattern) {
  message <- tryCatch({ force(expr); NULL }, error = conditionMessage)
  stopifnot(!is.null(message), grepl(pattern, message, fixed = TRUE))
}

wide <- data.frame(
  subject = strsplit("ATGCCATG", "")[[1]],
  position = 1:8,
  subject.position = 11:18,
  forward = c("A", "T", NA, NA, "C", "A", NA, NA),
  reverse = c(NA, NA, "G", "C", "C", NA, NA, NA),
  empty = rep(NA_character_, 8),
  check.names = FALSE
)
original <- wide
gbk <- wide_to_gbk(wide[c(8, 2, 6, 4, 5, 3, 7, 1), ], seq_id = "example",
                   organism = "synthetic construct", topology = "linear",
                   feature_type = c(empty = "misc_feature", reverse = "CDS",
                                    forward = "misc_feature"),
                   strand = c(reverse = "-", forward = "+", empty = "+"))
df <- features_to_df(gbk)
stopifnot(
  identical(wide, original), inherits(gbk, "gbk_file"),
  inherits(gbk[[1]], "gbk_record"), gbk[[1]]$locus$name == "example",
  gbk[[1]]$locus$topology == "linear", gbk[[1]]$sequence == "atgccatg",
  identical(df$type, c("source", "misc_feature", "CDS")),
  identical(df$location, c("1..8", "join(1..2,5..6)", "complement(3..5)")),
  identical(df$label, c(NA_character_, "forward", "reverse")),
  identical(df$organism, c("synthetic construct", NA_character_, NA_character_)),
  extract_feature_sequence(gbk[[1]], gbk[[1]]$features[[2]]) == "atca",
  extract_feature_sequence(gbk[[1]], gbk[[1]]$features[[3]]) == "ggc"
)

# Direct output and the ordinary writer both survive parsing.
path <- tempfile(fileext = ".gbk")
direct <- wide_to_gbk(wide, seq_id = "example", file = path)
parsed <- read_gbk(path)
stopifnot(identical(parsed[[1]]$sequence, direct[[1]]$sequence),
          identical(features_to_df(parsed), features_to_df(direct)),
          !"organism" %in% names(features_to_df(direct)),
          !grepl("NA", direct[[1]]$headers[[1]]$value, fixed = TRUE))
write_gbk(gbk, path)
stopifnot(identical(features_to_df(read_gbk(path)), features_to_df(gbk)))
unlink(path)

# Alternate columns, factors, whitespace, single bases, overlapping masks,
# reverse-strand joins, and explicitly selected metadata.
custom <- data.frame(base = factor(c("a", "n", "t", "g")),
                     a = c("x", "", "x", " "), b = c(NA, "x", "x", NA))
custom_gbk <- wide_to_gbk(custom, sequence_col = "base", position_col = NULL,
                          strand = "-")
stopifnot(custom_gbk[[1]]$sequence == "antg",
          identical(features_to_df(custom_gbk)$location,
                    c("1..4", "complement(join(1..1,3..3))", "complement(2..3)")))
source_only <- wide_to_gbk(wide, annotation_cols = character())
stopifnot(length(source_only[[1]]$features) == 1L,
          length(wide_to_gbk(wide, annotation_cols = "subject.position")[[1]]$features) == 2L)
markers <- data.frame(subject = c("A", "C"), position = 1:2,
                      logical = c(FALSE, NA), numeric = c(0, NA))
stopifnot(identical(features_to_df(wide_to_gbk(markers))$location,
                    c("1..2", "1..1", "1..1")))

# Reject data that would silently corrupt coordinates or sequence.
expect_error(wide_to_gbk(wide[FALSE, ]), "non-empty data frame")
expect_error(wide_to_gbk(wide, seq_id = "two names"), "without whitespace")
expect_error(wide_to_gbk(wide, organism = ""), "organism")
expect_error(wide_to_gbk(wide, sequence_col = "missing"), "sequence_col")
expect_error(wide_to_gbk(wide, position_col = "subject"), "must differ")
expect_error(wide_to_gbk(wide, annotation_cols = "subject"), "annotation_cols")
expect_error(wide_to_gbk(wide, annotation_cols = c("reverse", "reverse")), "annotation_cols")
expect_error(wide_to_gbk(wide, strand = "?"), "strand")
expect_error(wide_to_gbk(wide, strand = c(reverse = "-")), "every annotation")
expect_error(wide_to_gbk(wide, feature_type = NA_character_), "feature_type")
expect_error(wide_to_gbk(wide, feature_type = "bad key"), "feature_type")
expect_error(wide_to_gbk(wide, feature_type = c("CDS", "gene")), "every annotation")
for (bad in list(c(1:7, 7), c(2:9), c(1:7, NA), c(1:7, 8.5))) {
  broken <- wide
  broken$position <- bad
  expect_error(wide_to_gbk(broken), "each integer")
}
for (bad in c(NA_character_, "", "-", "AT", "U")) {
  broken <- wide
  broken$subject[1] <- bad
  expect_error(wide_to_gbk(broken), "IUPAC DNA")
}
broken <- wide
names(broken)[4] <- "subject"
expect_error(wide_to_gbk(broken), "unique")
broken <- wide
broken$forward <- as.list(broken$forward)
expect_error(wide_to_gbk(broken), "atomic vector")
