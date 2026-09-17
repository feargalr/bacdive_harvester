## 01a_target_list.R -- build the list of genera to harvest.
##
## Harvesting by GENUS rather than by species is the key change. BacDive's taxon
## endpoint accepts a genus, retrieve() paginates through the whole result, and
## one call returns every strain of every species in that genus. So:
##   * ~350 calls instead of 1,500+;
##   * species we never thought to ask for arrive for free;
##   * reclassified names resolve themselves, because we are not matching on the
##     species epithet at query time.
##
## Writes data/target_genera.csv. Run this once; 01_harvest.R reads the file.

source("R/00_config.R")

## Datasets to draw species names from. v7 used HMP + iHMP + LifeLines only;
## widening this is Phase 3.2 of the plan. Set TAXSEA_CMD_ALL=1 to sweep every
## relative_abundance dataset in curatedMetagenomicData instead.
CMD_PATTERNS <- c(
  "2021-03-31.HMP_2012.relative_abundance",
  "2021-03-31.HMP_2019_ibdmdb.relative_abundance",
  "2021-03-31.LifeLinesDeep_2016.relative_abundance",
  "2021-03-31.OhJ_2014.relative_abundance",
  "2021-03-31.QinJ_2012.relative_abundance",
  "2021-03-31.QinN_2014.relative_abundance",
  "2021-03-31.ZellerG_2014.relative_abundance",
  "2021-03-31.FengQ_2015.relative_abundance",
  "2021-03-31.AsnicarF_2017.relative_abundance",
  "2021-03-31.VogtmannE_2016.relative_abundance"
)

## One quality filter, shared by every source. Applied late so it cannot be
## bypassed: TaxSEA::NCBI_ids in particular carries `_sp_`/CAG placeholders,
## literal ".1" suffixes and non-bacterial entries, none of which BacDive holds.
## Vectorised and alignment-preserving: returns the canonical binomial, or NA for
## anything that is not one. Keeping the length lets callers filter a parallel
## vector of ids without the two drifting apart.
canon_species_or_na <- function(x) {
  x <- gsub("_", " ", as.character(x))
  x <- trimws(gsub("\\s+", " ", x))
  ## Strip the ".1"/".2" de-duplication artefacts present in TaxSEA's own lookup
  ## (e.g. "Slackia equolifaciens.1", "Eubacterium sp..1").
  x <- sub("[.][0-9]+$", "", x)
  x <- sub("[.]$", "", x)
  ## Placeholder and MAG names: BacDive only holds cultured, named strains, so
  ## these can never match and would only burn API calls.
  ## Junk words must be WHOLE words. An unanchored "bacterium" matches inside
  ## Faecali-bacterium, Eu-bacterium, Coryne-bacterium, Fuso-bacterium and
  ## Propioni-bacterium, which silently deleted most of the core gut genera.
  ## \\b does not fire inside those names because the preceding character is a
  ## word character, so the boundary form is both safe and still catches the
  ## placeholders ("Firmicutes bacterium", "gut metagenome").
  bad <- grepl("(^|\\s)sp(\\s|$|[.])", x) |
    grepl("\\b(bacterium|archaeon|metagenome|uncultured|unclassified|symbiont|endosymbiont)\\b",
          x, ignore.case = TRUE) |
    grepl("\\bCAG\\b", x) |
    grepl("oral taxon", x, ignore.case = TRUE) |
    ## Require a clean two-word binomial: Genus epithet. This also rejects the
    ## bare-numeric self-entries in NCBI_ids (id "820" appears as its own name).
    !grepl("^[A-Z][a-z]{2,} [a-z][a-z-]{2,}$", x)
  x[bad] <- NA_character_
  x
}

clean_species_names <- function(x) {
  x <- canon_species_or_na(x)
  sort(unique(x[!is.na(x)]))
}

## Parse a MetaPhlAn lineage string to a species binomial.
species_from_lineage <- function(x) {
  x <- x[grepl("s__", x, fixed = TRUE)]
  sp <- vapply(strsplit(x, "|", fixed = TRUE), function(p) {
    hit <- grep("^s__", p, value = TRUE)
    if (length(hit)) sub("^s__", "", hit[[1]]) else NA_character_
  }, character(1))
  unique(sp[!is.na(sp)])
}

collect_cmd_species <- function(patterns) {
  if (!requireNamespace("curatedMetagenomicData", quietly = TRUE)) {
    log_msg("curatedMetagenomicData not installed; skipping that source")
    return(character(0))
  }
  all_sp <- character(0)
  for (p in patterns) {
    log_msg("loading ", p)
    obj <- tryCatch(
      curatedMetagenomicData::curatedMetagenomicData(pattern = p, counts = TRUE,
                                                     dryrun = FALSE),
      error = function(e) { log_msg("  skipped: ", conditionMessage(e)); NULL })
    if (is.null(obj) || !length(obj)) next
    se <- obj[[1]]
    all_sp <- union(all_sp, species_from_lineage(rownames(se)))
  }
  all_sp
}

## The taxa already in TaxSEA_db are worth covering, because those are what users
## test against.
##
## NB: the source is the taxids that appear as SET MEMBERS in TaxSEA_db (~1,800),
## not the names in TaxSEA::NCBI_ids. NCBI_ids is a 12,772-entry name->id lookup
## covering far more than the database itself, and taking it wholesale pulled in
## placeholder names and non-bacterial entries.
species_from_taxsea <- function() {
  if (!requireNamespace("TaxSEA", quietly = TRUE)) return(character(0))
  e <- new.env()
  ok <- tryCatch({
    utils::data("TaxSEA_db", package = "TaxSEA", envir = e)
    utils::data("NCBI_ids",  package = "TaxSEA", envir = e)
    TRUE
  }, error = function(err) FALSE)
  if (!ok) return(character(0))

  db  <- get("TaxSEA_db", envir = e)
  ids <- get("NCBI_ids",  envir = e)
  member_ids <- unique(as.character(unlist(db, use.names = FALSE)))

  ## Invert the lookup to id -> name.
  ##
  ## NCBI_ids is NOT a 1:1 named vector: of its 12,772 elements, 360 are NULL and
  ## 38 hold 2-3 ids. So names(ids) is longer than unlist(ids) and pairing them
  ## directly silently scrambles the mapping -- it had id 820 pointing at
  ## "Dietzia maris". rep(names, lengths) is the correct expansion.
  nm  <- rep(names(ids), lengths(ids))
  val <- as.character(unlist(ids, use.names = FALSE))
  stopifnot(length(nm) == length(val))

  ## Reduce to names that are genuine binomials BEFORE picking one per id. Each
  ## id carries several spellings ("Bacteroides_uniformis", "Bacteroides
  ## uniformis", "820", "Bacteroides uniformis.1"); choosing by string length
  ## would pick the bare id.
  canon <- canon_species_or_na(nm)
  keep <- !is.na(canon) & !is.na(val) & nzchar(val)
  canon <- canon[keep]; val <- val[keep]

  first <- !duplicated(val)
  hit <- canon[first][match(member_ids, val[first])]
  unique(hit[!is.na(hit)])
}

## A user-supplied species list short-circuits everything. This is the route that
## needs no Bioconductor packages at all: drop a one-column CSV at
## data/target_species.csv (header `species`, values like "Bacteroides uniformis")
## and the genus list is derived from it directly.
USER_LIST <- file.path("data", "target_species.csv")
if (file.exists(USER_LIST) && nzchar(Sys.getenv("TAXSEA_USE_USER_LIST", ""))) {
  sp <- clean_species_names(utils::read.csv(USER_LIST, stringsAsFactors = FALSE)$species)
  genera <- sort(unique(sub(" .*$", "", sp)))
  utils::write.csv(data.frame(genus = genera), file.path("data", "target_genera.csv"),
                   row.names = FALSE)
  log_msg(sprintf("using supplied %s: %d species -> %d genera", USER_LIST,
                  length(sp), length(genera)))
  quit(status = 0L, save = "no")
}

## TAXSEA_CMD_ALL=1 sweeps every relative_abundance dataset instead of the list
## above. dryrun = TRUE returns the matching resource names rather than the data.
patterns <- if (nzchar(Sys.getenv("TAXSEA_CMD_ALL"))) {
  avail <- curatedMetagenomicData::curatedMetagenomicData("relative_abundance",
                                                          dryrun = TRUE)
  avail <- as.character(unlist(avail, use.names = FALSE))
  avail <- unique(avail[grepl("relative_abundance$", avail)])
  log_msg(sprintf("TAXSEA_CMD_ALL set: sweeping %d datasets", length(avail)))
  avail
} else CMD_PATTERNS

raw_cmd    <- collect_cmd_species(patterns)
raw_taxsea <- species_from_taxsea()
sp_cmd     <- clean_species_names(raw_cmd)
sp_taxsea  <- clean_species_names(raw_taxsea)
species    <- sort(unique(c(sp_cmd, sp_taxsea)))
if (!length(species)) {
  stop("no species found. Install curatedMetagenomicData and/or TaxSEA, or supply ",
       "data/target_species.csv and set TAXSEA_USE_USER_LIST=1.", call. = FALSE)
}
genera     <- sort(unique(sub(" .*$", "", species)))

log_msg(sprintf("curatedMetagenomicData: %d raw -> %d clean species",
                length(raw_cmd), length(sp_cmd)))
log_msg(sprintf("TaxSEA_db members:      %d raw -> %d clean species",
                length(raw_taxsea), length(sp_taxsea)))
log_msg(sprintf("union: %d species (%d only in TaxSEA_db, %d only in cMD)",
                length(species), length(setdiff(sp_taxsea, sp_cmd)),
                length(setdiff(sp_cmd, sp_taxsea))))

utils::write.csv(data.frame(genus = genera, stringsAsFactors = FALSE),
                 file.path("data", "target_genera.csv"), row.names = FALSE)
utils::write.csv(data.frame(species = sort(species), stringsAsFactors = FALSE),
                 file.path("data", "target_species.csv"), row.names = FALSE)

log_msg(sprintf("stage 01a done: %d species -> %d genera to harvest",
                length(species), length(genera)))
