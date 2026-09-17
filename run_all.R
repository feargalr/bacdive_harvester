## run_all.R -- drive the whole pipeline in order.
##
##   Rscript run_all.R            # stages 02-07 from the existing cache
##   Rscript run_all.R --harvest  # also run 01a/01 (needs BacDive credentials)
##
## R/01b_resolve_synonyms.R is deliberately not in this sequence. It is a second
## pass: run it after a first harvest to find target species that were missed
## because they have been reclassified into a genus you did not query, then merge
## data/target_genera_extra.csv into data/target_genera.csv and re-run the
## harvest (which skips genera already cached).
##
## Each stage reads the previous stage's cached artefact, so re-running any
## subset is safe and produces the same result.

args <- commandArgs(trailingOnly = TRUE)
do_harvest <- "--harvest" %in% args
## 05_validate.R exits non-zero on a real gold-standard failure, which gates the
## merge by default. Pass --ignore-validation to build the .rda anyway.
ignore_validation <- "--ignore-validation" %in% args

stages <- c(
  if (do_harvest) c("R/01a_target_list.R", "R/01_harvest.R"),
  "R/02_flatten.R",
  "R/03_aggregate.R",
  "R/04_build_sets.R",
  "R/05_validate.R",
  "R/07_export_review.R",
  ## Last, and only if TaxSEA is installed: everything above is standalone.
  if (requireNamespace("TaxSEA", quietly = TRUE)) "R/06_merge_taxsea_db.R"
)

for (s in stages) {
  cat("\n", strrep("=", 72), "\n", s, "\n", strrep("=", 72), "\n", sep = "")
  st <- system2("Rscript", shQuote(s))
  if (!identical(as.integer(st), 0L)) {
    if (grepl("05_validate", s) && ignore_validation) {
      cat("\n!! validation reported failures; continuing because --ignore-validation\n")
      next
    }
    stop(sprintf("stage %s exited with status %d", s, st),
         if (grepl("05_validate", s))
           "\n  fix the failing assertions, or re-run with --ignore-validation"
         else "",
         call. = FALSE)
  }
}
cat("\nall stages complete\n")
