library(gbkr)

make_record <- function(locations = character(), circular = FALSE,
                        sequence = "aacgttagccat") {
  structure(list(
    locus = list(name = "test", length = nchar(sequence), unit = "bp",
                 molecule_type = "DNA",
                 topology = if (circular) "circular" else "linear"),
    headers = list(list(key = "DEFINITION", prefix = "DEFINITION  ",
                        value = "Test record.")),
    features = lapply(locations, function(location) {
      list(type = "misc_feature", location = location,
           qualifiers = list(list(name = "label", value = "test", quoted = TRUE)))
    }),
    sequence = sequence
  ), class = c("gbk_record", "list"))
}

expect_error <- function(expr, pattern) {
  message <- tryCatch({ force(expr); NULL }, error = conditionMessage)
  stopifnot(!is.null(message), grepl(pattern, message, fixed = TRUE))
}

# Explicit coordinates detect off-by-one errors independently of extraction.
locations <- c("2..5", "complement(2..5)", "1", "12", "<2..>5", "<2",
               "5^6", "join(2..4,8..10)", "complement(join(2..4,8..10))",
               "order(2..4,8..10)", "join(complement(2..4),8..10)")
expected <- c("complement(8..11)", "8..11", "complement(12)", "complement(1)",
              "complement(<8..>11)", "complement(>11)", "complement(7^8)",
              "complement(join(3..5,9..11))", "join(3..5,9..11)",
              "complement(order(3..5,9..11))",
              "complement(join(3..5,complement(9..11)))")
record <- make_record(locations)
original <- record
reversed <- reverse_complement_gbk(record)
stopifnot(
  identical(record, original),
  identical(reversed$sequence, "atggctaacgtt"),
  identical(vapply(reversed$features, `[[`, "", "location"), expected),
  identical(extract_feature_sequences(record), extract_feature_sequences(reversed)),
  identical(reverse_complement_gbk(reversed), record),
  identical(reversed$headers, record$headers)
)

# DNA ambiguity codes and case must survive a double reversal.
ambiguous <- make_record(sequence = "ACGTRYSWKMBDHVNacgtryswkmbdhvn")
stopifnot(identical(reverse_complement_gbk(ambiguous)$sequence,
                    "nbdhvkmwsryacgtNBDHVKMWSRYACGT"),
          identical(reverse_complement_gbk(reverse_complement_gbk(ambiguous)),
                    ambiguous))

circular <- make_record(c("10..3", "complement(10..3)", "12^1",
                          "join(10..12,1..3)"), circular = TRUE)
flipped <- reverse_complement_gbk(circular)
stopifnot(identical(flipped$features[[1L]]$location, "complement(10..3)"),
          identical(flipped$features[[3L]]$location, "complement(12^1)"),
          identical(extract_feature_sequences(circular),
                    extract_feature_sequences(flipped)),
          identical(reverse_complement_gbk(flipped), circular))

# Standard qualifiers with coordinate or direction semantics must also change.
qualified <- make_record("2..10")
values <- list(direction = "LEFT", direction = "RIGHT", direction = "BOTH",
               transl_except = "(pos:complement(2..4),aa:Sec)",
               anticodon = "(pos:join(2..3,5),aa:Leu,seq:taa)",
               rpt_unit_range = "2..4", translation = "MLK",
               pseudo = NULL, note = "unchanged")
qualified$features[[1L]]$qualifiers <- Map(function(name, value) {
  list(name = name, value = value, quoted = TRUE)
}, names(values), unname(values))
q <- reverse_complement_gbk(qualified)$features[[1L]]$qualifiers
stopifnot(identical(unname(lapply(q, `[[`, "value")),
                    list("RIGHT", "LEFT", "BOTH", "(pos:9..11,aa:Sec)",
                         "(pos:complement(join(8,10..11)),aa:Leu,seq:taa)",
                         "9..11", "MLK", NULL, "unchanged")),
          identical(reverse_complement_gbk(reverse_complement_gbk(qualified)),
                    qualified))

# All entry points, including multi-record files and connections.
records <- structure(list(first = record, second = circular),
                     class = c("gbk_file", "list"))
stopifnot(identical(reverse_complement_gbk(reverse_complement_gbk(records)), records),
          identical(reverse_complement_gbk(list(record))[[1L]], reversed))
input <- tempfile(fileext = ".gbk")
output <- tempfile(fileext = ".gbk")
write_gbk(records, input)
before <- readLines(input)
from_path <- reverse_complement_gbk(input)
connection <- file(input, "r")
from_connection <- reverse_complement_gbk(connection)
close(connection)
stopifnot(identical(from_path, from_connection), identical(readLines(input), before))
write_gbk(from_path, output)
roundtrip <- read_gbk(output)
for (i in seq_along(records)) {
  stopifnot(identical(extract_feature_sequences(records[[i]]),
                      extract_feature_sequences(roundtrip[[i]])),
            identical(roundtrip[[i]]$features, from_path[[i]]$features))
}
unlink(c(input, output))

for (location in c("0..3", "2..13", "999999999999999999999999999999")) {
  expect_error(reverse_complement_gbk(make_record(location)), "outside")
}
expect_error(reverse_complement_gbk(make_record("10..3")), "circular")
expect_error(reverse_complement_gbk(make_record("12^1")), "adjacent")
expect_error(reverse_complement_gbk(make_record("2^4")), "adjacent")
expect_error(reverse_complement_gbk(make_record("OTHER.1:2..4")), "Remote-accession")
expect_error(reverse_complement_gbk(make_record("one-of(2,4)")), "Unsupported")
expect_error(reverse_complement_gbk(make_record("join(2..4,8..10")), "parentheses")
expect_error(reverse_complement_gbk(make_record("complement(2,4)")), "Invalid")
expect_error(reverse_complement_gbk(make_record("join(2..4,)")), "locations require")
expect_error(reverse_complement_gbk(make_record(NA_character_)), "non-empty string")
expect_error(reverse_complement_gbk(make_record(sequence = "")), "IUPAC DNA")
expect_error(reverse_complement_gbk(make_record(sequence = "acgu")), "IUPAC DNA")
expect_error(reverse_complement_gbk(make_record(sequence = "acgt-")), "IUPAC DNA")
expect_error(reverse_complement_gbk(list()), "x must be")
protein <- make_record()
protein$locus$unit <- "aa"
expect_error(reverse_complement_gbk(protein), "DNA record")

# A stale declared length must not affect mapping.
stale <- make_record("2..5")
stale$locus$length <- 100L
stopifnot(identical(reverse_complement_gbk(stale)$features[[1L]]$location,
                    "complement(8..11)"),
          identical(reverse_complement_gbk(stale)$locus$length, 12L))
