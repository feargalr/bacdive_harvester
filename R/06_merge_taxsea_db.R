## 06_merge_taxsea_db.R -- build TaxSEA's data files from the harvested sets.
##
## Reads an existing TaxSEA_db.rda and NCBI_ids.rda and writes updated copies to
## output/. Never modifies a package in place; copying the files into TaxSEA's
## data/ directory is a deliberate manual step.
##
## Source of the existing files, in order of preference:
##   TAXSEA_DATA_DIR  a directory holding TaxSEA_db.rda and NCBI_ids.rda, e.g. a
##                    TaxSEA source checkout's data/. Needs no TaxSEA install.
##   installed TaxSEA otherwise.
##
## What it does:
##   1. ships only sets with at least PIPE$ship_min_members members. TaxSEA
##      intersects a set with the observed taxa before size-filtering, so a set
##      smaller than this can never be tested at that threshold;
##   2. replaces every existing BacDive_* set and keeps every other source;
##   3. extends NCBI_ids ADD-ONLY with species name -> taxid pairs from the
##      harvest, plus former names resolved by 01b, in both "Genus species" and
##      "Genus_species" forms. Existing taxids are never overwritten, so a name
##      that maps today maps the same way afterwards; disagreements are reported
##      instead. Names present with no taxid (dead ends) are filled;
##   4. writes a per-set provenance table and a build record.

source("R/00_config.R")

min_members <- as.integer(Sys.getenv("TAXSEA_SHIP_MIN_MEMBERS",
                                     PIPE$ship_min_members))

## ---- inputs ------------------------------------------------------------------

b <- readRDS(PIPE$taxon_sets)
meta <- readRDS(PIPE$strain_meta)

data_dir <- Sys.getenv("TAXSEA_DATA_DIR")
e <- new.env()
if (nzchar(data_dir)) {
  for (f in c("TaxSEA_db.rda", "NCBI_ids.rda")) {
    p <- file.path(data_dir, f)
    if (!file.exists(p)) stop("missing ", p, call. = FALSE)
    load(p, envir = e)
  }
  log_msg("existing data read from ", data_dir)
} else if (requireNamespace("TaxSEA", quietly = TRUE)) {
  utils::data("TaxSEA_db", "NCBI_ids", package = "TaxSEA", envir = e)
  log_msg("existing data read from the installed TaxSEA package")
} else {
  stop("set TAXSEA_DATA_DIR to a directory containing TaxSEA_db.rda and ",
       "NCBI_ids.rda, or install TaxSEA", call. = FALSE)
}
TaxSEA_db <- e$TaxSEA_db
NCBI_ids <- e$NCBI_ids

## ---- 1. sets to ship -----------------------------------------------------------

sets <- b$sets
ship <- sets[lengths(sets) >= min_members]
log_msg(sprintf("%d BacDive sets built, %d have >= %d members and will ship",
                length(sets), length(ship), min_members))

prefix <- paste0(PIPE$set_prefix, "_")
if (!length(ship)) stop("no sets to ship", call. = FALSE)
if (!all(startsWith(names(ship), prefix))) {
  stop("set names must start with ", prefix, call. = FALSE)
}
if (any(grepl("[[:space:]]", names(ship)))) stop("set names contain whitespace", call. = FALSE)
if (any(vapply(ship, anyDuplicated, integer(1)) > 0L)) {
  stop("duplicated taxid within a set", call. = FALSE)
}

## ---- 2. merge ------------------------------------------------------------------

old <- grepl("BacDive", names(TaxSEA_db))
kept <- TaxSEA_db[!old]
clash <- intersect(names(kept), names(ship))
if (length(clash)) stop("new set names collide with existing ones: ",
                        paste(clash, collapse = ", "), call. = FALSE)
TaxSEA_db <- c(kept, ship)
log_msg(sprintf("TaxSEA_db: removed %d old BacDive sets, kept %d other sets, now %d",
                sum(old), length(kept), length(TaxSEA_db)))

## ---- 3. NCBI_ids, add-only -----------------------------------------------------

shipped_tx <- unique(unlist(ship, use.names = FALSE))

pairs <- unique(meta[!is.na(meta$ncbi_taxid) & meta$ncbi_taxid %in% shipped_tx,
                     c("species", "ncbi_taxid")])
names(pairs) <- c("name", "taxid")
pairs$source <- rep("harvest", nrow(pairs))

syn_file <- file.path("data", "synonym_resolution.csv")
if (file.exists(syn_file)) {
  syn <- utils::read.csv(syn_file, stringsAsFactors = FALSE, colClasses = "character")
  syn <- syn[!is.na(syn$taxid) & syn$taxid %in% shipped_tx, , drop = FALSE]
  if (nrow(syn)) {
    pairs <- rbind(pairs, data.frame(name = syn$missing_species, taxid = syn$taxid,
                                     source = "former name", stringsAsFactors = FALSE))
  }
}

## Both spellings TaxSEA users pass: "Genus species" and "Genus_species".
pairs$name <- trimws(gsub("_", " ", pairs$name))
## Binomials only. BacDive species fields also hold "Bacillus sp." and
## "X y subsp. z"; the first is not a species and the second would map a
## subspecies name onto its parent species' taxid.
pairs <- pairs[grepl("^[A-Z][a-z]+ [a-z][a-z-]+$", pairs$name), , drop = FALSE]
underscored <- pairs
underscored$name <- gsub(" ", "_", underscored$name)
pairs <- unique(rbind(pairs, underscored))

## A name the harvest itself maps to more than one taxid is ambiguous: skip it.
n_ids <- tapply(pairs$taxid, pairs$name, function(z) length(unique(z)))
ambiguous <- names(n_ids)[n_ids > 1L]
pairs <- pairs[!pairs$name %in% ambiguous, , drop = FALSE]
pairs <- pairs[!duplicated(pairs$name), , drop = FALSE]

in_lookup <- pairs$name %in% names(NCBI_ids)
existing <- lapply(pairs$name, function(n)
  if (n %in% names(NCBI_ids)) as.character(unlist(NCBI_ids[[n]])) else character(0))
## A name already present with no taxid is a dead end in TaxSEA: get_NCBI_sets()
## will not query NCBI for a name it holds, so it silently maps to nothing.
## Filling it adds a mapping and changes no existing one.
empty <- in_lookup & lengths(existing) == 0L
agrees <- vapply(seq_len(nrow(pairs)), function(i) pairs$taxid[i] %in% existing[[i]],
                 logical(1))
conflicts <- pairs[in_lookup & !empty & !agrees, , drop = FALSE]
conflicts$existing_taxid <- vapply(existing[in_lookup & !empty & !agrees],
                                   paste, character(1), collapse = "|")
additions <- pairs[!in_lookup | empty, , drop = FALSE]
additions$action <- ifelse(additions$name %in% names(NCBI_ids), "filled empty entry", "added")

reach_before <- mean(shipped_tx %in% as.character(unlist(NCBI_ids)))
fill <- additions[additions$action != "added", , drop = FALSE]
for (i in seq_len(nrow(fill))) NCBI_ids[[fill$name[i]]] <- fill$taxid[i]
add <- additions[additions$action == "added", , drop = FALSE]
NCBI_ids <- c(NCBI_ids, stats::setNames(as.list(add$taxid), add$name))
reach_after <- mean(shipped_tx %in% as.character(unlist(NCBI_ids)))
log_msg(sprintf("NCBI_ids: +%d names, %d empty entries filled (%d from former names); %d conflicts left untouched; %d ambiguous skipped",
                nrow(add), nrow(fill), sum(additions$source == "former name"),
                nrow(conflicts), length(ambiguous)))
log_msg(sprintf("shipped BacDive taxa reachable by name: %.0f%% -> %.0f%%",
                100 * reach_before, 100 * reach_after))

## ---- 4. provenance -------------------------------------------------------------

prov <- utils::read.csv(PIPE$set_provenance, stringsAsFactors = FALSE)
prov <- prov[prov$set_name %in% names(ship), , drop = FALSE]
prov$n_taxa <- unname(lengths(ship)[prov$set_name])

man <- if (file.exists(PIPE$manifest_file)) {
  utils::read.csv(PIPE$manifest_file, stringsAsFactors = FALSE)
} else NULL
harvest_dates <- c(NA_character_, NA_character_)
if (!is.null(man)) {
  d <- sort(as.Date(substr(man$harvested_at, 1, 10)))
  if (length(d)) harvest_dates <- as.character(range(d))
}

commit <- Sys.getenv("TAXSEA_HARVESTER_COMMIT")
if (!nzchar(commit)) {
  commit <- tryCatch(system2("git", c("rev-parse", "--short", "HEAD"),
                             stdout = TRUE, stderr = FALSE)[1],
                     error = function(err) "unknown",
                     warning = function(w) "unknown")
}

build <- data.frame(
  key = c("harvester", "harvester_commit", "bacdive_r_package",
          "harvest_first_date", "harvest_last_date", "genera_queried",
          "strain_records", "species_with_calls", "sets_built",
          "sets_shipped", "ship_min_members", "taxa_in_shipped_sets",
          "ncbi_ids_added", "built_on"),
  value = c("https://github.com/feargalr/bacdive_harvester", commit,
            if (!is.null(man)) paste(sort(unique(man$bacdive_pkg)), collapse = "|") else NA,
            harvest_dates[1], harvest_dates[2],
            if (!is.null(man)) nrow(man) else NA,
            if (!is.null(man)) sum(man$n_records, na.rm = TRUE) else NA,
            length(unique(meta$species)),
            length(sets), length(ship), min_members, length(shipped_tx),
            nrow(additions), as.character(Sys.Date())),
  stringsAsFactors = FALSE)

## ---- write -----------------------------------------------------------------------

out <- PIPE$output_dir
save(TaxSEA_db, file = file.path(out, "TaxSEA_db.rda"), compress = "xz", version = 3)
save(NCBI_ids, file = file.path(out, "NCBI_ids.rda"), compress = "xz", version = 3)
tsv <- function(x, f) utils::write.table(x, file.path(out, f), sep = "\t",
                                         quote = FALSE, row.names = FALSE, na = "")
tsv(prov, "BacDive_set_provenance.tsv")
tsv(build, "BacDive_build_info.tsv")
tsv(additions, "NCBI_ids_additions.tsv")
tsv(conflicts, "NCBI_ids_conflicts.tsv")

log_msg(sprintf("stage 06 done: TaxSEA_db.rda (%d KB), NCBI_ids.rda (%d KB) in %s",
                round(file.size(file.path(out, "TaxSEA_db.rda")) / 1024),
                round(file.size(file.path(out, "NCBI_ids.rda")) / 1024), out))
log_msg("copy both into TaxSEA's data/ directory by hand")
