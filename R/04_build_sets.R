## 04_build_sets.R -- species-level trait calls -> named taxon sets.
## Writes the set list, a per-set provenance table and the tested universe.
## Does NOT touch the installed TaxSEA_db; merging is a separate explicit step.

source("R/00_config.R")
source("R/lib/vocab.R")
source("R/lib/normalise.R")
source("R/lib/sets.R")

traits <- readRDS(PIPE$species_traits)
spec   <- load_trait_spec()

before <- length(unique(paste(traits$set_infix, traits$analyte)))
traits <- collapse_analytes(traits)
after  <- length(unique(paste(traits$set_infix, traits$analyte)))
log_msg(sprintf("cross-family analyte collapse: %d -> %d distinct analytes", before, after))

## emit_negative must be known before the information filter runs.
em <- stats::setNames(tolower(spec$emit_negative) %in% c("yes","true","1"), spec$family)
traits$emit_negative <- unname(em[traits$family])
traits$emit_negative[is.na(traits$emit_negative)] <- FALSE
traits <- drop_uninformative(traits, spec)

traits <- name_sets(traits, spec)
if (!nrow(traits)) stop("no sets after naming", call. = FALSE)

## Unions are added AFTER naming, because they carry their own set_name.
traits <- add_derived_unions(traits)

sets <- build_taxon_sets(traits)
prov <- set_provenance(traits, sets)
univ <- build_tested_universe(traits)

no_taxid <- attr(sets, "species_without_taxid")
if (length(no_taxid)) {
  log_msg(sprintf("%d species dropped for lack of a species-level NCBI taxid", length(no_taxid)))
  writeLines(sort(no_taxid), file.path(PIPE$output_dir, "04_species_without_taxid.txt"))
}

saveRDS(list(sets = sets, tested_universe = univ, traits = traits), PIPE$taxon_sets)
utils::write.csv(prov, PIPE$set_provenance, row.names = FALSE)

log_msg(sprintf("stage 04 done: %d sets, %d unique taxa, %d measured / %d predicted sets",
                length(sets), length(unique(unlist(sets))),
                sum(prov$evidence == "measured"), sum(prov$evidence == "predicted")))
