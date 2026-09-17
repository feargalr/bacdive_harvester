## 01_harvest.R -- pull raw BacDive records into a local cache, one .rds per genus.
##
## This is the ONLY stage that touches the network. Everything downstream reads
## the cache, which is what makes the pipeline deterministic: the v7 scripts
## re-queried on every run and two runs a day apart disagreed on the oxygen call
## for 300 of 730 species.
##
## Credentials come from the environment, never from source. Copy
## .Renviron.example to .Renviron (which is gitignored) and fill it in.
##
## Resumable: a genus whose .rds already exists is skipped, so an interrupted
## harvest can be restarted without re-querying.

source("R/00_config.R")

if (!requireNamespace("BacDive", quietly = TRUE)) {
  stop("the BacDive package is required: install.packages('BacDive')", call. = FALSE)
}

user <- Sys.getenv("BACDIVE_USER")
pass <- Sys.getenv("BACDIVE_PASSWORD")
if (!nzchar(user) || !nzchar(pass)) {
  stop("BACDIVE_USER / BACDIVE_PASSWORD are not set.\n",
       "  cp .Renviron.example .Renviron   # then fill it in and restart R\n",
       "Register for API access at https://api.bacdive.dsmz.de", call. = FALSE)
}

genus_file <- file.path("data", "target_genera.csv")
if (!file.exists(genus_file)) {
  stop("missing ", genus_file, " -- run R/01a_target_list.R first", call. = FALSE)
}
genera <- utils::read.csv(genus_file, stringsAsFactors = FALSE)$genus
genera <- sort(unique(genera[nzchar(genera)]))

## TAXSEA_MAX_GENERA=10 harvests only the first N genera. Use it to smoke-test
## against the live API before committing to the full sweep.
lim <- suppressWarnings(as.integer(Sys.getenv("TAXSEA_MAX_GENERA", "")))
if (!is.na(lim) && lim > 0L) {
  genera <- utils::head(genera, lim)
  log_msg(sprintf("TAXSEA_MAX_GENERA=%d: limiting to %d genera", lim, length(genera)))
}
log_msg(sprintf("%d genera to harvest", length(genera)))

bd <- BacDive::open_bacdive(username = user, password = pass)

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

## One retry with backoff. retrieve() paginates internally; a mid-pagination
## failure loses the whole genus, so it is worth retrying once before giving up.
harvest_genus <- function(genus, attempts = 3L) {
  for (k in seq_len(attempts)) {
    res <- tryCatch(
      BacDive::retrieve(bd, query = genus, search = "taxon", sleep = PIPE$api_sleep),
      error = function(e) structure(list(msg = conditionMessage(e)), class = "harvest_error"))
    if (!inherits(res, "harvest_error")) return(res)
    if (k < attempts) {
      wait <- 2^k
      log_msg(sprintf("  %s failed (%s); retrying in %ds", genus, res$msg, wait))
      Sys.sleep(wait)
    } else {
      log_msg(sprintf("  %s FAILED after %d attempts: %s", genus, attempts, res$msg))
      return(NULL)
    }
  }
  NULL
}

manifest <- list()
pkg_ver  <- as.character(utils::packageVersion("BacDive"))

for (g in genera) {
  out <- file.path(PIPE$cache_dir, paste0(safe_name(g), ".rds"))
  if (file.exists(out)) {
    blob <- readRDS(out)
    manifest[[length(manifest) + 1L]] <- data.frame(
      genus = g, n_records = blob$n_records %||% NA_integer_,
      harvested_at = as.character(blob$harvested_at %||% NA),
      bacdive_pkg = blob$bacdive_pkg %||% NA_character_,
      status = "cached", stringsAsFactors = FALSE)
    next
  }

  recs <- harvest_genus(g)
  n <- if (is.null(recs)) 0L else length(recs)

  ## Record what actually came back. BacDive's taxon endpoint truncates a query
  ## to its first three word components and does not guarantee the result matches
  ## what you asked for, so we log the genera present rather than assume.
  genera_seen <- if (n) {
    sort(unique(vapply(recs, function(r) {
      v <- r[["Name and taxonomic classification"]]$genus
      if (is.null(v) || !length(v)) NA_character_ else as.character(v)[[1]]
    }, character(1))))
  } else character(0)

  saveRDS(list(genus = g, harvested_at = Sys.time(), n_records = n,
               bacdive_pkg = pkg_ver, genera_seen = genera_seen,
               records = if (is.null(recs)) list() else recs), out)

  manifest[[length(manifest) + 1L]] <- data.frame(
    genus = g, n_records = n, harvested_at = as.character(Sys.time()),
    bacdive_pkg = pkg_ver,
    status = if (is.null(recs)) "failed" else if (!n) "empty" else "ok",
    stringsAsFactors = FALSE)

  log_msg(sprintf("%-28s %4d records  [%s]", g, n, paste(genera_seen, collapse = ", ")))
  Sys.sleep(PIPE$api_sleep)
}

man <- do.call(rbind, manifest)
utils::write.csv(man, PIPE$manifest_file, row.names = FALSE)

log_msg(sprintf("stage 01 done: %d genera, %d records, %d failed, %d empty",
                nrow(man), sum(man$n_records, na.rm = TRUE),
                sum(man$status == "failed"), sum(man$status == "empty")))
if (any(man$status == "failed")) {
  log_msg("re-run this script to retry the failed genera (cached ones are skipped)")
}
