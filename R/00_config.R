## 00_config.R -- shared configuration for the TaxSEA BacDive pipeline.
## Sourced by every stage. Contains no credentials and makes no network calls.

## Directories can be redirected with environment variables so the test suite can
## run the real stages against tests/fixtures/cache without touching real data.
.cache_dir  <- Sys.getenv("TAXSEA_CACHE_DIR",  unset = "cache")
.output_dir <- Sys.getenv("TAXSEA_OUTPUT_DIR", unset = "output")

PIPE <- list(
  cache_dir     = .cache_dir,        # raw BacDive JSON, one .rds per genus
  output_dir    = .output_dir,       # derived tables and the final set list
  log_dir       = "logs",
  curation_dir  = file.path("data", "curation"),

  ## Stage artefacts. Each stage reads the previous stage's file and never the API.
  manifest_file = file.path(.cache_dir,  "harvest_manifest.csv"),
  strain_long   = file.path(.output_dir, "02_strain_long.rds"),
  strain_meta   = file.path(.output_dir, "02_strain_meta.rds"),
  species_traits= file.path(.output_dir, "03_species_traits.rds"),
  taxon_sets    = file.path(.output_dir, "04_taxon_sets.rds"),
  set_provenance= file.path(.output_dir, "04_set_provenance.csv"),
  validation    = file.path(.output_dir, "05_validation.csv"),
  coverage      = file.path(.output_dir, "05_coverage.csv"),

  ## Aggregation thresholds (strain -> species).
  ## A species is called positive for a binary assay when at least this fraction
  ## of the strains that were *tested* came back positive.
  min_pos_fraction = 0.5,
  ## ...and at least this many strains were tested.
  min_strains_tested = 1L,

  ## A `dominant` categorical value (see oxygen_values.csv) overrides majority
  ## voting, but only with real support behind it: at least this fraction of the
  ## strains tested, and at least `dominant_min_n` of them. Without the floor, one
  ## stray annotation flips a species -- Pseudomonas aeruginosa became a
  ## facultative anaerobe on 1 strain out of 142.
  dominant_min_frac = 0.10,
  dominant_min_n    = 2L,

  ## Information-content filter. A polarity set covering more than this share of
  ## its own tested universe is not a biological hypothesis -- it is "almost
  ## everything that was tested". Uses_L-xylose_acid_negative held 98% of the
  ## species ever tested for L-xylose, and 38 such near-universal negative sets
  ## reached significance on the IBD data while spanning only 70 distinct taxa
  ## between them. Applies to assay families that emit both polarities.
  max_universe_fraction = 0.90,
  ## ...but only once enough species were tested for the share to mean anything.
  ## With a universe of one, every set is trivially 100% of it, which would
  ## delete small well-defined traits rather than near-universal ones.
  min_universe_for_filter = 10L,

  ## Set naming. The BacDive_ prefix must be retained: TaxSEA identifies and
  ## replaces this family with grepl("BacDive", names(TaxSEA_db)).
  set_prefix        = "BacDive",
  predicted_infix   = "Pred",   # genome-based / confidence-scored calls

  ## Politeness. retrieve() enforces a floor of 0.1s between pages.
  api_sleep = 0.35
)

for (d in c(PIPE$cache_dir, PIPE$output_dir, PIPE$log_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

## ---- small shared utilities -------------------------------------------------

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

log_msg <- function(...) {
  message(format(Sys.time(), "[%H:%M:%S] "), ...)
}

## Trim, collapse internal whitespace, drop empties. Used everywhere before
## comparing or keying on a free-text value.
clean_str <- function(x) {
  x <- gsub("\\s+", " ", trimws(as.character(x)))
  x[!nzchar(x) | x %in% c("NA", "na")] <- NA_character_
  x
}

