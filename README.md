# gbkr

`gbkr` is a dependency-free R package for reading, editing, and writing GenBank
flat files (`.gb` and `.gbk`). It represents records as ordinary R lists and
provides a data-frame interface for working with feature annotations.

## Installation

Install the development version from GitHub:

```r
pak::pak("Close-your-eyes/gbkr")
```

To install a local checkout instead:

```r
install.packages("/path/to/gbkr", repos = NULL, type = "source")
```

## Quick start

Read a GenBank file, edit its first record, and write the result:

```r
library(gbkr)

gbk <- read_gbk("input.gbk")
gbk

record <- gbk[[1]]
record$locus$name <- "edited_record"

feature_df <- features_to_df(record)
feature_df <- add_gbk_feature(
  feature_df,
  type = "promoter",
  location = "100..200",
  label = "new_promoter",
  note = "added in R"
)

record$features <- features_from_df(feature_df)
gbk[[1]] <- record

write_gbk(gbk, "edited.gbk")
```

`read_gbk()` accepts either a file path or a readable connection. Files may
contain one or more records separated by the standard `//` terminator.

## Record structure

`read_gbk()` returns a `gbk_file`, which is a list of `gbk_record` objects. Each
record has four editable components:

```text
record
├── locus       Parsed fields from the LOCUS line
├── headers     Ordered GenBank header entries
├── features    Ordered feature definitions and qualifiers
└── sequence    Lower-case sequence string
```

For example:

```r
record <- gbk[[1]]

record$locus
record$headers
record$features
record$sequence
```

Each feature contains a `type`, a GenBank `location`, and an ordered list of
qualifiers. Qualifiers are stored as separate list entries so repeated
qualifiers can be represented without losing their order. A flag qualifier such
as `/pseudo` has a `NULL` value.

## Work with features as a data frame

Convert a record's feature list to a wide, character-only data frame:

```r
feature_df <- features_to_df(record)
feature_df
```

The table has one row per feature. `type` and `location` are followed by one
column for each qualifier.

- `NA` means that a qualifier is absent.
- `"<flag>"` represents a flag qualifier with no value.
- Repeated qualifier values are stored in one cell, separated by a newline.
- Printing a `gbk_feature_df` displays repeated values separated by ` | `.

You can edit the table with ordinary data-frame operations:

```r
feature_df$note[1] <- "reviewed annotation"
feature_df$gene[feature_df$gene == "old_name"] <- "new_name"
```

Add a feature with `add_gbk_feature()`:

```r
feature_df <- add_gbk_feature(
  feature_df,
  type = "CDS",
  location = "complement(250..900)",
  gene = "exampleA",
  note = c("predicted", "manually reviewed"),
  pseudo = NULL
)
```

A vector creates a repeated qualifier, while `NULL` creates a flag qualifier.
Use `.quoted` to override automatic qualifier quoting:

```r
feature_df <- add_gbk_feature(
  feature_df,
  type = "misc_feature",
  location = "950..1000",
  label = "site_1",
  .quoted = c(label = FALSE)
)
```

Convert the edited table back to the feature-list representation before
writing the record:

```r
record$features <- features_from_df(feature_df)
```

Tables created by `features_to_df()` retain attributes describing repeated
values, flags, and whether each qualifier was quoted. If those attributes are
removed, `features_from_df()` applies conventional GenBank quoting rules.

## Extract feature sequences

Extract a sequence using a feature object, a one-row feature data frame, or a
location string:

```r
extract_feature_sequence(record, "10..50")
extract_feature_sequence(record, "complement(10..50)")
extract_feature_sequence(record, record$features[[1]])
```

Extract several features at once:

```r
sequences <- extract_feature_sequences(record)
sequences
```

By default, sequences are named using the first available `label`, `gene`, or
`locus_tag` qualifier. You can change that order or omit names:

```r
extract_feature_sequences(
  record,
  name_by = c("locus_tag", "gene"),
  use_names = TRUE
)
```

Supported location syntax includes:

- single positions and ranges, such as `10` and `10..50`;
- fuzzy bounds marked with `<` or `>`;
- `join()`, `order()`, and `complement()`;
- descending ranges across the origin of circular records; and
- between-base sites such as `10^11`, which return an empty sequence.

Remote-accession locations, such as `J00194.1:100..200`, cannot be extracted
from the local record sequence.

## Reverse-complement a GenBank file

Reverse the entire DNA sequence and update feature coordinates and strands:

```r
reversed <- reverse_complement_gbk("input.gbk")
write_gbk(reversed, "reverse_complement.gbk")
```

You can also pass an object returned by `read_gbk()`, a single record, or a list
of records. The input is left unchanged. Each feature still describes the same
biological sequence, including joined features and features across a circular
origin. Feature order and translations are retained.

The function supports IUPAC DNA bases and the local location syntax listed
above. It also updates LEFT/RIGHT direction qualifiers and coordinates in
`transl_except`, `anticodon`, and `rpt_unit_range` qualifiers. Unsupported or
remote-accession locations produce an error. Coordinates embedded in free-text
notes or headers are not rewritten.

## Main functions

| Function | Purpose |
| --- | --- |
| `read_gbk()` | Read one or more GenBank records |
| `write_gbk()` | Write records in GenBank flat-file format |
| `reverse_complement_gbk()` | Reverse-complement DNA and remap feature locations |
| `features_to_df()` | Convert parsed features to a wide data frame |
| `features_from_df()` | Convert a feature data frame back to parsed features |
| `add_gbk_feature()` | Append a feature and its qualifiers |
| `extract_feature_sequence()` | Extract the sequence for one feature |
| `extract_feature_sequences()` | Extract sequences for several features |

## Notes

- `write_gbk()` validates records by default and updates the `LOCUS` sequence
  length from the sequence.
- Output formatting is normalized, so a read-write round trip is not intended
  to reproduce the original file byte for byte.
- The package focuses on common local GenBank location expressions and does not
  download sequence data referenced by remote accessions.
