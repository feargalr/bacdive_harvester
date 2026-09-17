## 05_validate.R -- assertions and reports over the built sets.
## Exits non-zero if any gold-standard assertion actually fails, so this can gate
## a release. Species absent from the harvest are reported, not failed.

source("R/00_config.R")
source("R/lib/normalise.R")
source("R/lib/validate.R")

b      <- readRDS(PIPE$taxon_sets)
sets   <- b$sets
traits <- b$traits
prov   <- utils::read.csv(PIPE$set_provenance, stringsAsFactors = FALSE)
gold   <- load_curation("gold_standard.csv", c("species", "trait_set", "expectation"))

gs <- check_gold_standard(traits, sets, gold)
utils::write.csv(gs, PIPE$validation, row.names = FALSE)

cov <- coverage_report(traits, sets, prov)
utils::write.csv(cov, PIPE$coverage, row.names = FALSE)

depth <- annotation_depth(traits)
utils::write.csv(depth, file.path(PIPE$output_dir, "05_annotation_depth.csv"), row.names = FALSE)

red <- set_redundancy(sets, threshold = 0.9)
if (!is.null(red)) {
  utils::write.csv(red, file.path(PIPE$output_dir, "05_redundant_sets.csv"), row.names = FALSE)
}

cat("\n=== coverage ===\n");  print(cov, row.names = FALSE)
cat("\n=== gold standard ===\n"); print(table(gs$status))

nd <- gs[gs$status == "no_data", , drop = FALSE]
if (nrow(nd)) {
  cat("\n--- BacDive records nothing in this trait family (coverage ceiling, not a defect) ---\n")
  print(nd[, c("species", "trait_set")], row.names = FALSE)
}

fails <- gs[gs$status == "fail", , drop = FALSE]
if (nrow(fails)) {
  cat("\n--- FAILED assertions (BacDive has a call and it contradicts) ---\n")
  print(fails[, c("species", "trait_set", "expectation", "observed")], row.names = FALSE)
}
nc <- gs[gs$status == "not_covered", , drop = FALSE]
if (nrow(nc)) {
  cat(sprintf("\n%d assertion(s) skipped: species not present in this harvest (%s)\n",
              nrow(nc), paste(unique(nc$species), collapse = ", ")))
}

if (!is.null(red)) {
  cat(sprintf("\n=== %d set pair(s) with Jaccard >= 0.9 (candidate unmerged synonyms) ===\n",
              nrow(red)))
  print(utils::head(red, 20), row.names = FALSE)
} else {
  cat("\n=== no set pairs above Jaccard 0.9 ===\n")
}

cat("\n=== annotation depth: most and least annotated species ===\n")
print(utils::head(depth, 5), row.names = FALSE)
print(utils::tail(depth, 5), row.names = FALSE)

if (nrow(fails)) {
  quit(status = 1L, save = "no")
}
