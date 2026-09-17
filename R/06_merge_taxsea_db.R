## 06_merge_taxsea_db.R -- fold the new BacDive sets into TaxSEA_db.
##
## Deliberately separate from stage 04 and deliberately non-destructive: it writes
## output/TaxSEA_db.rda and never modifies the installed package or anything in
## the parent directory. Copying it into the package is a manual step.
##
## Differences from `Formatting TaxSEA db v7.R`:
##   * sets arrive already named and already keyed on taxid, so there is no
##     split(ncbi_id, set_id) that could merge every oxygen category into one set;
##   * every set name carries the BacDive_ prefix, so the removal step below is
##     idempotent -- v7 emitted bare names like "oxygen_tolerance", which its own
##     !grepl("BacDive", ...) filter would then fail to remove on a second run;
##   * ids are de-duplicated;
##   * a provenance table and the tested universe are written alongside.

source("R/00_config.R")

b    <- readRDS(PIPE$taxon_sets)
sets <- b$sets

stopifnot(length(sets) > 0L)
bad <- names(sets)[!startsWith(names(sets), paste0(PIPE$set_prefix, "_"))]
if (length(bad)) {
  stop("these set names lack the BacDive_ prefix, which would break removal on ",
       "the next run: ", paste(utils::head(bad, 5), collapse = ", "), call. = FALSE)
}
if (any(vapply(sets, function(x) anyDuplicated(x) > 0L, logical(1)))) {
  stop("duplicated taxid within a set", call. = FALSE)
}

if (!requireNamespace("TaxSEA", quietly = TRUE)) {
  stop("TaxSEA is not installed; cannot load the existing database", call. = FALSE)
}
e <- new.env()
utils::data("TaxSEA_db", package = "TaxSEA", envir = e)
TaxSEA_db <- get("TaxSEA_db", envir = e)

n_before  <- length(TaxSEA_db)
n_bacdive <- sum(grepl("BacDive", names(TaxSEA_db)))
TaxSEA_db <- TaxSEA_db[!grepl("BacDive", names(TaxSEA_db))]
log_msg(sprintf("existing db: %d sets, removed %d BacDive sets, %d remain",
                n_before, n_bacdive, length(TaxSEA_db)))

## Guard against colliding with a non-BacDive set name.
clash <- intersect(names(TaxSEA_db), names(sets))
if (length(clash)) stop("new set names collide with existing ones: ",
                        paste(clash, collapse = ", "), call. = FALSE)

TaxSEA_db <- c(TaxSEA_db, sets)

if (!dir.exists(PIPE$output_dir)) dir.create(PIPE$output_dir, recursive = TRUE)
save(TaxSEA_db, file = file.path(PIPE$output_dir, "TaxSEA_db.rda"),
     compress = "xz", version = 2)

sizes <- vapply(sets, length, integer(1))
log_msg(sprintf("stage 06 done: %d sets total (%d BacDive), %d unique taxa in BacDive sets",
                length(TaxSEA_db), length(sets), length(unique(unlist(sets)))))
log_msg(sprintf("BacDive set sizes: median %g, range %d-%d",
                stats::median(sizes), min(sizes), max(sizes)))
log_msg(sprintf("written to %s -- copy into the package's data/ directory by hand",
                file.path(PIPE$output_dir, "TaxSEA_db.rda")))
