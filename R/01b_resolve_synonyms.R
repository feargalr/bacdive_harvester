## 01b_resolve_synonyms.R -- find target species the genus harvest missed because
## they have been reclassified into a genus we never queried.
##
## The problem: harvesting by genus resolves synonyms in one direction only.
## Querying "Absiella" returns records BacDive files under *Eubacterium*, so that
## works. But *Eubacterium rectale* is now *Agathobacter rectalis*: we queried
## genus "Eubacterium", BacDive files the organism under "Agathobacter", and a
## taxon search on the old genus does not return the new one. 289 of 1,851 target
## species (15.6%) were missed this way.
##
## The fix: for every missing species, ask NCBI Taxonomy for its CURRENT
## scientific name plus all its synonyms, take the genera implied by those names,
## and add any we have not harvested to a supplementary genus list. 01_harvest.R
## then picks them up on the next run (it skips genera already cached).
##
## Uses only the public NCBI E-utilities; no credentials, no BacDive calls.
## Writes data/target_genera_extra.csv and data/synonym_resolution.csv.

source("R/00_config.R")

EUTILS <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
## NCBI asks for <=3 requests/second without an API key. Set NCBI_API_KEY in
## .Renviron to raise that to 10/s.
NCBI_KEY <- Sys.getenv("NCBI_API_KEY")
NCBI_SLEEP <- if (nzchar(NCBI_KEY)) 0.12 else 0.35

.key_param <- function() if (nzchar(NCBI_KEY)) paste0("&api_key=", NCBI_KEY) else ""

.get <- function(url, attempts = 3L) {
  for (k in seq_len(attempts)) {
    res <- tryCatch(suppressWarnings(readLines(url, warn = FALSE)),
                    error = function(e) NULL)
    if (!is.null(res)) return(paste(res, collapse = "\n"))
    Sys.sleep(2^k)
  }
  NULL
}

## Resolve a species name to an NCBI taxid. NCBI keeps former names searchable,
## which is exactly the property we need.
##
## XML rather than JSON so no JSON parser is needed: esearch returns
## <IdList><Id>1382</Id></IdList>, and the efetch parsing below is already
## regex-based, so this keeps the whole script on base R.
ncbi_taxid_for_name <- function(name) {
  url <- sprintf("%s/esearch.fcgi?db=taxonomy&term=%s&retmode=xml%s",
                 EUTILS, utils::URLencode(name, reserved = TRUE), .key_param())
  txt <- .get(url)
  if (is.null(txt)) return(NA_character_)
  m <- regmatches(txt, regexpr("<Id>[0-9]+</Id>", txt))
  if (!length(m)) return(NA_character_)
  gsub("[^0-9]", "", m[[1]])
}

## Current scientific name + every synonym / equivalent name, for a batch of
## taxids. efetch accepts a comma-separated list, so this is a handful of calls.
ncbi_names_for_taxids <- function(taxids, batch = 150L) {
  taxids <- unique(taxids[!is.na(taxids) & nzchar(taxids)])
  if (!length(taxids)) return(NULL)
  out <- list()
  chunks <- split(taxids, ceiling(seq_along(taxids) / batch))
  for (ci in seq_along(chunks)) {
    ids <- chunks[[ci]]
    url <- sprintf("%s/efetch.fcgi?db=taxonomy&id=%s&retmode=xml%s",
                   EUTILS, paste(ids, collapse = ","), .key_param())
    txt <- .get(url)
    Sys.sleep(NCBI_SLEEP)
    if (is.null(txt)) {
      log_msg(sprintf("  efetch batch %d failed", ci))
      next
    }
    ## One <Taxon> block per record. Split on the top-level opening tag.
    blocks <- strsplit(txt, "<Taxon>", fixed = TRUE)[[1]][-1]
    for (b in blocks) {
      tid <- sub(".*?<TaxId>([0-9]+)</TaxId>.*", "\\1", b)
      sn  <- if (grepl("<ScientificName>", b)) {
        sub(".*?<ScientificName>(.*?)</ScientificName>.*", "\\1", b)
      } else NA_character_
      syn <- unlist(regmatches(b, gregexpr("<Synonym>.*?</Synonym>", b)))
      eqv <- unlist(regmatches(b, gregexpr("<EquivalentName>.*?</EquivalentName>", b)))
      alt <- gsub("<[^>]+>", "", c(syn, eqv))
      ## Division/Lineage let us drop non-prokaryotes. The target list contains
      ## fungi (Aspergillus kawachii, Saccharomyces) whose synonyms would
      ## otherwise add fungal genera that BacDive cannot hold.
      dv <- if (grepl("<Division>", b)) {
        sub(".*?<Division>(.*?)</Division>.*", "\\1", b)
      } else NA_character_
      lin <- if (grepl("<Lineage>", b)) {
        sub(".*?<Lineage>(.*?)</Lineage>.*", "\\1", b)
      } else NA_character_
      out[[length(out) + 1L]] <- data.frame(
        taxid = tid, current_name = sn,
        alt_names = paste(unique(alt), collapse = " | "),
        division = dv, lineage = lin,
        stringsAsFactors = FALSE)
    }
    log_msg(sprintf("  efetch batch %d/%d (%d taxids)", ci, length(chunks), length(ids)))
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

## ---- which target species did we miss? --------------------------------------

target <- utils::read.csv(file.path("data", "target_species.csv"),
                          stringsAsFactors = FALSE)$species
meta <- readRDS(PIPE$strain_meta)
found <- unique(meta$species)
harvested_genera <- unique(sub(" .*$", "", found))
missing <- sort(setdiff(target, found))

log_msg(sprintf("%d target species, %d found, %d missing", length(target),
                length(target) - length(missing), length(missing)))
if (!length(missing)) {
  log_msg("nothing to resolve")
  quit(status = 0L, save = "no")
}

## ---- taxids: prefer TaxSEA's lookup, fall back to esearch -------------------

## TaxSEA is only a shortcut: it saves an NCBI round-trip for species whose taxid
## is already in its lookup. Everything it provides is obtainable from esearch, so
## the script works without it -- just more slowly.
taxid_of <- stats::setNames(rep(NA_character_, length(missing)), missing)
if (requireNamespace("TaxSEA", quietly = TRUE)) {
  e <- new.env()
  ok <- tryCatch({ utils::data("NCBI_ids", package = "TaxSEA", envir = e); TRUE },
                 error = function(err) FALSE)
  if (ok) {
    ids <- get("NCBI_ids", envir = e)
    nm  <- rep(names(ids), lengths(ids))
    val <- as.character(unlist(ids, use.names = FALSE))
    nm_clean <- trimws(gsub("_", " ", nm))
    hit <- match(missing, nm_clean)
    taxid_of[!is.na(hit)] <- val[hit[!is.na(hit)]]
  }
}
log_msg(sprintf("%d taxids from the TaxSEA lookup, %d need an NCBI search",
                sum(!is.na(taxid_of)), sum(is.na(taxid_of))))

need <- names(taxid_of)[is.na(taxid_of)]
for (i in seq_along(need)) {
  taxid_of[[need[[i]]]] <- ncbi_taxid_for_name(need[[i]])
  Sys.sleep(NCBI_SLEEP)
  if (i %% 25 == 0) log_msg(sprintf("  esearch %d/%d", i, length(need)))
}
log_msg(sprintf("%d of %d missing species now have a taxid",
                sum(!is.na(taxid_of)), length(taxid_of)))

## ---- current names and synonyms ---------------------------------------------

info <- ncbi_names_for_taxids(taxid_of)
if (is.null(info)) stop("no names returned from NCBI", call. = FALSE)

res <- data.frame(missing_species = names(taxid_of),
                  taxid = unname(taxid_of), stringsAsFactors = FALSE)
res <- merge(res, info, by = "taxid", all.x = TRUE)

## Every genus implied by the current name or any synonym.
## "Candidatus" is a nomenclatural status prefix, not a genus: the genus of
## "Candidatus Arthromitus" is Arthromitus. Querying BacDive for "Candidatus" or
## "Bacterium" would match an enormous, meaningless result set.
GENUS_BLOCKLIST <- c("Bacterium", "Bacillus sensu", "Bovine", "Einheimischer",
                     "Candidatus", "Uncultured", "Unidentified", "Endosymbiont",
                     "Symbiont", "Pseudobacterium", "Lactobacterium")

genus_of <- function(x) {
  x <- trimws(unlist(strsplit(x %||% "", "\\|")))
  x <- sub("^Candidatus\\s+", "", x)
  x <- x[grepl("^[A-Z][a-z]{2,} ", x)]
  g <- unique(sub(" .*$", "", x))
  setdiff(g, GENUS_BLOCKLIST)
}
res$candidate_genera <- vapply(seq_len(nrow(res)), function(i) {
  g <- unique(c(genus_of(res$current_name[i]), genus_of(res$alt_names[i])))
  paste(g, collapse = " | ")
}, character(1))

## BacDive holds bacteria and archaea only, so a fungal or metazoan synonym can
## never yield a record; including them just burns API calls.
is_prok <- (!is.na(res$division) & res$division %in% c("Bacteria", "Archaea")) |
  grepl("^(Bacteria|Archaea)\\b", res$lineage %||% "")
log_msg(sprintf("%d of %d resolved taxa are bacteria/archaea (%d dropped as non-prokaryote)",
                sum(is_prok), nrow(res), sum(!is_prok)))

prok <- res[is_prok, , drop = FALSE]
extra <- unique(unlist(lapply(seq_len(nrow(prok)), function(i) {
  unique(c(genus_of(prok$current_name[i]), genus_of(prok$alt_names[i])))
})))
extra <- sort(setdiff(extra, harvested_genera))

utils::write.csv(res[order(res$missing_species), ],
                 file.path("data", "synonym_resolution.csv"), row.names = FALSE)
utils::write.csv(data.frame(genus = extra, stringsAsFactors = FALSE),
                 file.path("data", "target_genera_extra.csv"), row.names = FALSE)

log_msg(sprintf("stage 01b done: %d new genera to harvest", length(extra)))
if (length(extra)) {
  log_msg(paste("  ", paste(utils::head(extra, 20), collapse = ", ")))
  log_msg("  merge into data/target_genera.csv, then re-run R/01_harvest.R")
}
