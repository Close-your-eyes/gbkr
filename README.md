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

## Create a GenBank file from a per-base wide table

`wide_to_gbk()` accepts one row per base, with the DNA in `subject`, positions
in `position`, and one column per annotation. Non-missing, non-empty annotation
cells mark feature positions. By default, `subject.position` is ignored.

```r
# wide_df is your per-base data frame.
gbk <- wide_to_gbk(
  wide_df,
  seq_id = "my_construct",
  organism = "synthetic construct",
  topology = "circular"
)
write_gbk(gbk, "my_construct.gbk")

# Or write directly while still returning the gbk_file object:
gbk <- wide_to_gbk(wide_df, seq_id = "my_construct", file = "my_construct.gbk")
```

Each annotation becomes a `misc_feature` labelled with its column name.
Separated runs within one column are joined; empty columns are skipped.
Choose `annotation_cols` to exclude extra metadata. `feature_type` and `strand`
accept either one value for every annotation or a vector named for every
selected annotation column:

```r
gbk <- wide_to_gbk(
  wide_df,
  seq_id = "my_construct",
  annotation_cols = c("NotI", "EcoRI"),
  feature_type = "misc_feature",
  strand = c(NotI = "+", EcoRI = "-")
)
```

Use `sequence_col` and `position_col` for different column names, or
`position_col = NULL` to use row order. Otherwise, positions must contain every
integer from one through the sequence length exactly once; rows are sorted
before conversion. Sequence cells must each contain one IUPAC DNA base.

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

## Edit sequence and annotations together

`splice_gbk()` updates DNA, feature coordinates, full-length source features,
and record length in one operation. It accepts an object from `read_gbk()` or
a single record and returns an edited copy.

```r
gbk <- read_gbk("input.gbk")

# Delete bases 100 through 200, including any wholly deleted features.
gbk <- splice_gbk(gbk, start = 100, end = 200)

# Insert before position 100 in the CURRENT edited sequence.
new_features <- data.frame(
  type = "misc_feature",
  location = "1..6",
  label = "new_insert"
)
gbk <- splice_gbk(
  gbk, start = 100, replacement = "ATGCGT", features = new_features
)
write_gbk(gbk, "edited.gbk")
```

Build the new features directly from a named vector of DNA parts:

```r
parts <- c(NotI = "gcggccgc", Kozak = "gccacc", insert = "ATGAGC", EcoRI = "gaattc")
new_features <- features_from_sequences(parts)
gbk <- splice_gbk(
  gbk,
  start = 100,
  replacement = paste0(parts, collapse = ""),
  features = new_features
)
```

`features_from_sequences()` returns a feature data frame with consecutive
locations and labels taken from the vector names. It defaults to `misc_feature`
on the forward strand. Set `feature_type` or `strand` to one value for all
parts, or a vector named for every part. Negative strands change annotation
orientation only; supplied DNA is concatenated as written. Keep the helper's
`start = 1` when passing features to `splice_gbk()`.

To replace an interval in one step, supply `start`, `end`, and `replacement`
in the same call. New feature coordinates refer to the replacement DNA,
starting at one. To append, use the current sequence length plus one as
`start` and omit `end`. For multi-record files, choose the record with
`record = 2`, for example.

Features after the edit shift automatically. Features that are only partly
affected cause an error by default. Set `overlap = "drop"` to remove those
features entirely, or `overlap = "trim"` to keep their surviving original
bases. Trimming excludes inserted DNA from old annotations and can produce
joined locations. It removes `translation`, `protein_id`, and `codon_start`
qualifiers from trimmed features and warns that annotations need review.
Translations and reading frames are not recalculated.

Joined and reverse-strand locations are supported, as are coordinate qualifiers
`transl_except`, `anticodon`, and `rpt_unit_range`. The edit interval itself must
not wrap around a circular origin. Other metadata, including accession/version
and free-text coordinates, is preserved verbatim; review it for the edited
construct. Call `write_gbk()` only after completing your edits.

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
| `wide_to_gbk()` | Build a GenBank record from per-base wide annotations |
| `splice_gbk()` | Insert, delete, or replace DNA and remap annotations |
| `reverse_complement_gbk()` | Reverse-complement DNA and remap feature locations |
| `features_to_df()` | Convert parsed features to a wide data frame |
| `features_from_sequences()` | Make consecutive features from named DNA parts |
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
