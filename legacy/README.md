# legacy

The original single-pass BacDive scraper, kept for reference and provenance. It
is superseded by the pipeline in `R/` and is not maintained.

| File | What it was |
|---|---|
| `bacive_trait_extractor.R` | queried BacDive one species at a time and extracted six traits from a single record per species |
| `helper_functions.R` | the trait accessors it used |
| `README_original.md` | the repository README as it stood before the rewrite |

## Why it was replaced

- **One record per species.** It kept the first type strain `retrieve()` returned
  and discarded the rest, although the API had already paginated through every
  strain. That makes the result depend on API result ordering, and it cannot
  distinguish "tested and negative" from "never tested".
- **No verification of what came back.** BacDive's taxon endpoint truncates a
  query to its first three word components and does not promise the result
  matches the query; the species name was taken from the query string rather than
  from the record.
- **Synonyms not reconciled.** `D-glucose` and `glucose`, or `esterase lipase
  (C 8)` and `esterase Lipase (C 8)`, became separate outputs.
- **Assay context discarded.** Growth on a substrate, acid production from it and
  gas production from it were merged into one trait.
- **A length-fragile type-strain test.** `!is.null(val) && tolower(val) == "yes"`
  is a hard error on R >= 4.3 when the field has length > 1, and it sat outside
  the loop's `tryCatch`.

The current pipeline addresses each of these; see the root `README.md`.
