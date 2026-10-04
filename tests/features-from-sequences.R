library(gbkr)

parts <- c(NotI = "gcggccgc", Kozak = "gccacc", insert = "ATGAGC", EcoRI = "gaattc")
original <- parts
features <- features_from_sequences(parts)
stopifnot(identical(parts, original), is.data.frame(features),
          identical(features$label, names(parts)),
          identical(features$type, rep("misc_feature", 4)),
          identical(features$location, c("1..8", "9..14", "15..20", "21..26")))
mapped <- features_from_sequences(parts, start = 101,
  strand = c(EcoRI = "-", insert = "+", Kozak = "+", NotI = "-"),
  feature_type = c(insert = "CDS", NotI = "misc_feature", Kozak = "misc_feature", EcoRI = "misc_feature"))
stopifnot(identical(mapped$location,
                    c("complement(101..108)", "109..114", "115..120", "complement(121..126)")),
          mapped$type[3] == "CDS",
          features_from_sequences(c(one = "N"))$location == "1..1",
          features_from_sequences(c(one = "AT"), start = 1e6)$location == "1000000..1000001")

# New feature extraction after splicing must exactly recover each supplied part.
gbk <- wide_to_gbk(data.frame(subject = c("A", "C", "G", "T"), position = 1:4))
edited <- splice_gbk(gbk, 3, replacement = paste0(parts, collapse = ""), features = features)
stopifnot(identical(unname(extract_feature_sequences(edited[[1]], edited[[1]]$features[-1])),
                    unname(parts)))
path <- tempfile(fileext = ".gbk")
write_gbk(edited, path)
reread <- read_gbk(path)
stopifnot(identical(features_to_df(reread), features_to_df(edited)),
          identical(tolower(reread[[1]]$sequence), tolower(edited[[1]]$sequence)))
unlink(path)

expect_error <- function(expr, pattern) {
  message <- tryCatch({ force(expr); NULL }, error = conditionMessage)
  stopifnot(!is.null(message), grepl(pattern, message, fixed = TRUE))
}
for (bad in list(character(), c(a = ""), c(a = NA_character_), c(a = "AT G"),
                 c(a = "AT-"), c(a = "AU"), list(a = "AT"))) {
  expect_error(features_from_sequences(bad), "IUPAC DNA")
}
for (bad in list("AT", c(a = "AT", a = "GC"), c(" " = "AT"))) {
  expect_error(features_from_sequences(bad), "names")
}
expect_error(features_from_sequences(parts, start = 0), "start")
expect_error(features_from_sequences(parts, start = 1.5), "start")
expect_error(features_from_sequences(parts, start = .Machine$integer.max), "integer range")
expect_error(features_from_sequences(parts, strand = "?"), "strand")
expect_error(features_from_sequences(parts, strand = c(NotI = "-")), "every sequence")
expect_error(features_from_sequences(parts, feature_type = "bad key"), "feature_type")
