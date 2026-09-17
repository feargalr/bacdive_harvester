## 07_export_review.R -- export the built sets as TSVs for manual curation.
##
## Three files, in increasing detail:
##
##   output/review_sets_summary.tsv   one row per set. The triage sheet: decide
##                                    keep / drop / merge / rename here.
##   output/review_analytes.tsv       one row per analyte, aggregating its
##                                    positive and negative sets side by side.
##                                    Use this to spot synonyms still unmerged.
##   output/review_sets_long.tsv.gz   the full membership table, one row per
##                                    (set, taxon), gzipped for storage.
##   output/review_sets_long.tsv      the same table uncompressed, so it opens
##                                    directly in Excel/Numbers without a
##                                    decompression step.
##
## TSV rather than CSV because analyte names contain commas and parentheses;
## nothing here contains a tab.

source("R/00_config.R")

b      <- readRDS(PIPE$taxon_sets)
sets   <- b$sets
traits <- b$traits
univ   <- b$tested_universe

write_tsv <- function(x, path, gz = FALSE) {
  con <- if (gz) gzfile(path, "w") else file(path, "w")
  on.exit(close(con))
  utils::write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE,
                     na = "")
  invisible(path)
}

## Analyte names must not contain a tab or newline or the TSV breaks.
scrub <- function(x) gsub("[\t\r\n]+", " ", x)

## ---- 1. one row per set -----------------------------------------------------

t2 <- traits[!is.na(traits$ncbi_taxid) & nzchar(traits$ncbi_taxid), , drop = FALSE]
by_set <- split(seq_len(nrow(t2)), t2$set_name)

## A few example species per set makes the triage sheet readable without
## cross-referencing the membership table.
example_species <- function(i) {
  sp <- sort(unique(t2$species[i]))
  paste(utils::head(sp, 3), collapse = "; ")
}

rows <- lapply(names(by_set), function(nm) {
  i <- by_set[[nm]]
  n_tested_universe <- length(univ[[paste(t2$set_infix[i][[1]],
                                          t2$analyte[i][[1]],
                                          t2$context[i][[1]], sep = "|")]] %||% character(0))
  data.frame(
    set_name        = nm,
    family          = paste(sort(unique(t2$family[i])), collapse = "|"),
    set_infix       = t2$set_infix[i][[1]],
    analyte         = scrub(t2$analyte[i][[1]]),
    context         = t2$context[i][[1]],
    context_raw     = scrub(paste(sort(unique(unlist(strsplit(
                        stats::na.omit(t2$context_raw[i]), "|", fixed = TRUE)))),
                        collapse = "|")),
    polarity        = if (grepl("_negative$", nm)) "negative" else "positive",
    evidence        = if (all(t2$predicted[i])) "predicted" else "measured",
    is_api_code     = grepl("^BacDive_API", nm),
    n_taxa          = length(sets[[nm]] %||% character(0)),
    n_species       = length(unique(t2$species[i])),
    n_species_tested= n_tested_universe,
    median_n_tested = stats::median(t2$n_tested[i], na.rm = TRUE),
    max_n_tested    = suppressWarnings(max(t2$n_tested[i], na.rm = TRUE)),
    example_species = scrub(example_species(i)),
    keep            = "",   # <- for you: keep / drop / merge
    merge_into      = "",   # <- for you: target set_name if merging
    notes           = "",
    stringsAsFactors = FALSE
  )
})
summ <- do.call(rbind, rows)
summ <- summ[order(-summ$n_taxa, summ$set_name), , drop = FALSE]
write_tsv(summ, file.path(PIPE$output_dir, "review_sets_summary.tsv"))

## ---- 2. one row per analyte, polarities side by side ------------------------
## The synonym-hunting view: two rows with near-identical n_pos over the same
## infix are very likely the same assay under two spellings.

key <- paste(t2$set_infix, t2$analyte, t2$context, sep = "\r")
ag <- lapply(split(seq_len(nrow(t2)), key), function(i) {
  data.frame(
    set_infix   = t2$set_infix[i][[1]],
    analyte     = scrub(t2$analyte[i][[1]]),
    context     = t2$context[i][[1]],
    is_api_code = grepl("^API", t2$set_infix[i][[1]]),
    families    = paste(sort(unique(t2$family[i])), collapse = "|"),
    n_species_pos = sum(t2$call[i] == "pos"),
    n_species_neg = sum(t2$call[i] == "neg"),
    n_species_tested = length(unique(t2$species[i])),
    pct_positive = round(100 * sum(t2$call[i] == "pos") / length(i), 1),
    canonical_name = "",  # <- for you: the name this analyte should collapse to
    notes          = "",
    stringsAsFactors = FALSE
  )
})
an <- do.call(rbind, ag)
an <- an[order(an$set_infix, -an$n_species_tested, an$analyte), , drop = FALSE]
write_tsv(an, file.path(PIPE$output_dir, "review_analytes.tsv"))

## ---- 3. full membership ------------------------------------------------------

## context_raw preserves the original BacDive assay contexts that were collapsed
## into each feature, so the granularity is recoverable without re-harvesting.
memb <- t2[, c("set_name", "ncbi_taxid", "species", "family", "analyte", "context",
               "context_raw", "call", "n_pos", "n_neg", "n_tested")]
memb$analyte <- scrub(memb$analyte)
memb <- memb[order(memb$set_name, memb$species), , drop = FALSE]
write_tsv(memb, file.path(PIPE$output_dir, "review_sets_long.tsv.gz"), gz = TRUE)
## Uncompressed twin: same content, opens without a decompression step.
write_tsv(memb, file.path(PIPE$output_dir, "review_sets_long.tsv"))

log_msg(sprintf("stage 07 done: %d sets, %d analytes, %d memberships exported",
                nrow(summ), nrow(an), nrow(memb)))
log_msg("  output/review_sets_summary.tsv   -- triage sheet (keep/drop/merge columns)")
log_msg("  output/review_analytes.tsv       -- synonym hunting (canonical_name column)")
log_msg("  output/review_sets_long.tsv.gz   -- full membership (gzipped)")
log_msg("  output/review_sets_long.tsv      -- same, uncompressed")
