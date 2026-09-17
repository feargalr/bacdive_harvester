## 03_aggregate.R -- strain-level long table -> species-level trait calls.

source("R/00_config.R")
source("R/lib/vocab.R")
source("R/lib/normalise.R")
source("R/lib/aggregate.R")

long <- readRDS(PIPE$strain_long)
meta <- readRDS(PIPE$strain_meta)
spec <- load_trait_spec()
log_msg(sprintf("aggregating %d trait families over %d strains / %d species",
                nrow(spec), nrow(meta), length(unique(meta$species))))
traits <- aggregate_all(long, meta, spec)
if (is.null(traits)) stop("no species-level trait calls produced", call. = FALSE)

## Attach the species-level NCBI taxid, taken from the records themselves.
tax <- unique(meta[!is.na(meta$ncbi_taxid), c("species", "ncbi_taxid")])
tax <- tax[!duplicated(tax$species), , drop = FALSE]
traits <- merge(traits, tax, by = "species", all.x = TRUE)

saveRDS(traits, PIPE$species_traits)
log_msg(sprintf("stage 03 done: %d species-trait calls, %d species, %d without a taxid",
                nrow(traits), length(unique(traits$species)),
                length(unique(traits$species[is.na(traits$ncbi_taxid)]))))
