## 02_flatten.R -- cached BacDive records -> one tidy long table + strain metadata.
## Reads only the local cache; makes no network calls.

source("R/00_config.R")
source("R/lib/flatten.R")

res <- flatten_cache()
saveRDS(res$long, PIPE$strain_long)
saveRDS(res$meta, PIPE$strain_meta)
log_msg(sprintf("stage 02 done: %d long rows, %d strains, %d species",
                nrow(res$long), nrow(res$meta), length(unique(res$meta$species))))
