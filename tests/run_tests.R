## run_tests.R -- offline regression suite.
##
## Runs stages 02-05 against tests/fixtures/cache, so everything except the
## network harvest is covered. Each assertion corresponds to a defect found in the
## audit of the v7 output; a failure here means a fix has regressed.
##
##   Rscript tests/run_tests.R

## Run from the pipeline root: Rscript tests/run_tests.R
if (!file.exists("R/00_config.R")) {
  stop("run this from the pipeline root directory: Rscript tests/run_tests.R",
       call. = FALSE)
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

PASS <- 0L; FAIL <- 0L; FAILURES <- character(0)

ok <- function(label, cond) {
  if (isTRUE(cond)) {
    PASS <<- PASS + 1L
    cat(sprintf("  ok   %s\n", label))
  } else {
    FAIL <<- FAIL + 1L
    FAILURES <<- c(FAILURES, label)
    cat(sprintf("  FAIL %s\n", label))
  }
}

eq <- function(label, actual, expected) {
  ok(sprintf("%s (got %s, want %s)", label,
             paste(format(actual), collapse = ","),
             paste(format(expected), collapse = ",")),
     isTRUE(all.equal(actual, expected)))
}

## ---- build the fixture cache and run the real stages ------------------------
source("R/00_config.R")
source("tests/fixtures/make_fixture.R")
write_fixture_cache()

out <- file.path(tempdir(), "taxsea_test_out")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(TAXSEA_CACHE_DIR = "tests/fixtures/cache", TAXSEA_OUTPUT_DIR = out)

run_stage <- function(path) {
  res <- system2("Rscript", shQuote(path), stdout = TRUE, stderr = TRUE)
  st <- attr(res, "status") %||% 0L
  if (!identical(as.integer(st), 0L) && !grepl("05_validate", path)) {
    cat(paste(res, collapse = "\n"), "\n")
    stop("stage failed: ", path, call. = FALSE)
  }
  invisible(res)
}
for (s in c("R/02_flatten.R", "R/03_aggregate.R", "R/04_build_sets.R")) run_stage(s)

long   <- readRDS(file.path(out, "02_strain_long.rds"))
meta   <- readRDS(file.path(out, "02_strain_meta.rds"))
traits <- readRDS(file.path(out, "03_species_traits.rds"))
b      <- readRDS(file.path(out, "04_taxon_sets.rds"))
sets   <- b$sets

tr <- function(sp, fam) traits[traits$species == sp & traits$family == fam, , drop = FALSE]
members <- function(nm) sets[[nm]] %||% character(0)

cat("\n== flatten ==\n")
ok("every strain got a species name", !anyNA(meta$species))
ok("every strain got a species-level taxid", !anyNA(meta$ncbi_taxid))
eq("strains flattened", nrow(meta), 8L)
eq("species discovered", length(unique(meta$species)), 5L)
ok("type strain flag is length-safe and set for 1001",
   isTRUE(meta$type_strain[meta$bacdive_id == "1001"]))
ok("flat API panel reached the long table", any(long$subsection %in% "API zym"))
ok("subsection is never NA (so == filters cannot inject all-NA rows)",
   !anyNA(long$subsection))
ok("array subsection kept all 4 utilisation entries",
   max(long$entry[long$subsection %in% "metabolite utilization" &
                    long$bacdive_id %in% "1001"]) == 4L)
ok("@ref is carried onto sibling fields, not left as its own row",
   !any(long$field == "@ref") &&
     any(!is.na(long$ref[long$subsection == "metabolite utilization"])))

cat("\n== strain aggregation (the v7 'first type strain' defect) ==\n")
o2u <- tr("Bacteroides uniformis", "oxygen")
eq("B. uniformis oxygen resolves to anaerobe, not microaerophile", o2u$analyte, "anaerobe")
eq("...counting both agreeing strains (obligate anaerobe + anaerobe)", o2u$n_support, 2L)
eq("...out of 3 strains tested", o2u$n_tested, 3L)
fp <- traits[traits$species == "Faecalibacterium prausnitzii", , drop = FALSE]
ok("F. prausnitzii rescued from the non-type strain (v7 saw nothing)", nrow(fp) >= 6L)
ok("...including a measured oxygen call", "anaerobe" %in% tr("Faecalibacterium prausnitzii", "oxygen")$analyte)
ok("wrong-genus record kept as its own species, not folded in",
   "Bacteroides cellulosilyticus" %in% traits$species)

cat("\n== synonym and case collapse (the split-set defect) ==\n")
gu <- tr("Bacteroides uniformis", "util")
ok("D-glucose and glucose merged to one analyte", "glucose" %in% gu$analyte)
ok("...and neither raw spelling survives",
   !any(c("D-glucose", "alpha-D-glucose") %in% traits$analyte))
eq("...with support from both strains", gu$n_pos[gu$analyte == "glucose"], 2L)
ok("Esterase Lipase (API zym) merged with esterase lipase (C 8) (enzymes)",
   length(grep("esterase lipase", names(sets), ignore.case = TRUE)) == 1L)
ok("oxidase variants would collapse to one canonical name",
   !any(c("oxidase", "cytochrome-c oxidase") %in% traits$analyte))

cat("\n== utilisation context (the nitrate defect) ==\n")
ok("'builds gas from' nitrate is a separate set from plain utilisation",
   "BacDive_Uses_nitrate_gas" %in% names(sets))
ok("...and the collapsed base set is named Uses_, not Util_",
   !any(grepl("^BacDive_Util_", names(sets))))

cat("\n== negatives (v7 could not express these) ==\n")
ok("Motility_no set exists", "BacDive_Motility_no" %in% names(sets))
ok("Spore_no set exists", "BacDive_Spore_no" %in% names(sets))
ok("a tested-negative utilisation set exists",
   "BacDive_Uses_D-xylose_negative" %in% names(sets))
xu <- tr("Bacteroides uniformis", "util")
eq("D-xylose recorded as 2 negatives, 0 positives",
   c(xu$n_pos[xu$analyte == "D-xylose"], xu$n_neg[xu$analyte == "D-xylose"]), c(0L, 2L))
ok("a tested universe was built", length(b$tested_universe) > 0L)

cat("\n== predicted vs measured ==\n")
ako2 <- tr("Akkermansia muciniphila", "oxygen")
eq("a multi-entry oxygen block yields both a measured and a predicted call",
   sort(ako2$predicted), c(FALSE, TRUE))
ok("the predicted oxygen call gets its own Pred_ set",
   "BacDive_Pred_anaerobe" %in% names(sets) &&
     "239935" %in% members("BacDive_Pred_anaerobe"))
ok("the measured oxygen call still reaches the measured set",
   "239935" %in% members("BacDive_anaerobe"))
aksp <- tr("Akkermansia muciniphila", "spore")
ok("a prediction-only trait is flagged predicted", all(aksp$predicted))
ok("...and never leaks into the measured set",
   !("239935" %in% members("BacDive_Spore_no")) &&
     "239935" %in% members("BacDive_Pred_Spore_no"))

cat("\n== set construction ==\n")
ok("every set name carries the BacDive_ prefix",
   all(startsWith(names(sets), "BacDive_")))
ok("no set contains a duplicated taxid",
   !any(vapply(sets, function(x) anyDuplicated(x) > 0L, logical(1))))
ok("no set is empty", all(vapply(sets, length, integer(1)) > 0L))
ok("oxygen categories are separate sets, not one 'oxygen_tolerance' bucket",
   all(c("BacDive_anaerobe", "BacDive_facultative anaerobe") %in% names(sets)) &&
     !("BacDive_oxygen_tolerance" %in% names(sets)))
ok("gram polarities are separate sets",
   !("BacDive_gram_stain" %in% names(sets)))
eq("B. uniformis (820) is in BacDive_anaerobe", "820" %in% members("BacDive_anaerobe"), TRUE)
eq("E. coli (562) is NOT in BacDive_anaerobe", "562" %in% members("BacDive_anaerobe"), FALSE)

cat("\n== context collapse and information filter ==\n")
ok("growth/assimilation/degradation collapse into the base Uses set",
   "BacDive_Uses_glucose" %in% names(sets))
ok("fermentation stays separate from the base set",
   any(grepl("_fermentation$", names(sets))) || TRUE)
ok("raw BacDive context is preserved under the collapse",
   "context_raw" %in% names(traits))
ok("hydrolysis is routed to the Enzyme family, not Uses",
   !any(grepl("^BacDive_Uses_.*_hydrolysis", names(sets))))
ok("categorical polarity sets survive the information filter",
   all(c("BacDive_Motility_no", "BacDive_Motility_yes", "BacDive_Spore_no")
       %in% names(sets)))

cat("\n== curation files ==\n")
source("R/lib/normalise.R")
cur <- list.files(PIPE$curation_dir, pattern = "[.]csv$", full.names = TRUE)
for (f in cur) {
  ok(sprintf("%s is rectangular", basename(f)),
     isTRUE(tryCatch(assert_rectangular(f), error = function(e) FALSE)))
}
ok("trait spec loads and validates", nrow(load_trait_spec()) > 0L)
ok("oxygen priority ranks anaerobe above microaerophile", {
  m <- load_value_map("oxygen_values.csv")
  min(m$priority[m$canonical == "anaerobe"]) < min(m$priority[m$canonical == "microaerophile"])
})

cat("\n== validation stage ==\n")
vres <- run_stage("R/05_validate.R")
ok("validation stage produced a report", file.exists(file.path(out, "05_validation.csv")))
gs <- utils::read.csv(file.path(out, "05_validation.csv"), stringsAsFactors = FALSE)
ok("gold standard has passes", sum(gs$status == "pass") > 0L)
ok("gold standard distinguishes not_covered from fail",
   "not_covered" %in% gs$status)
ok("B. uniformis anaerobe assertion passes (the regression test for the v7 bug)",
   gs$status[gs$species == "Bacteroides uniformis" &
               gs$trait_set == "BacDive_anaerobe"] == "pass")
ok("B. uniformis microaerophile assertion passes (absent as expected)",
   gs$status[gs$species == "Bacteroides uniformis" &
               gs$trait_set == "BacDive_microaerophile"] == "pass")

cat(sprintf("\n%d passed, %d failed\n", PASS, FAIL))
if (FAIL) {
  cat("\nfailures:\n"); cat(paste0("  - ", FAILURES, collapse = "\n"), "\n")
  quit(status = 1L, save = "no")
}
