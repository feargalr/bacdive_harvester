## make_fixture.R -- synthetic BacDive records that reproduce the real record
## shapes observed in https://api.bacdive.dsmz.de/v2/example/fetch/24493 plus the
## specific failure modes found in the audit of the v7 output.
##
## Deliberately exercises:
##   * a subsection arriving as a single named list  (cell morphology)
##   * a subsection arriving as an unnamed list of entries (metabolite utilization)
##   * a flat API panel block  (API zym) -- unused by v7, carries explicit negatives
##   * strains of one species disagreeing on oxygen tolerance
##     (the Bacteroides uniformis anaerobe/microaerophile case)
##   * a species whose type strain carries no physiology at all, while a
##     non-type strain does (the Faecalibacterium prausnitzii case)
##   * D-glucose vs glucose and "Esterase Lipase" vs "esterase lipase"
##     (the synonym / case-variant set splits)
##   * "kind of utilization tested" = "builds gas from" on nitrate, which v7
##     folded into BacDive_Utilizes_nitrate
##   * a confidence-scored (predicted) oxygen call that must not be merged with
##     lab-measured calls
##   * a record from the wrong genus, to prove the species guard works

.tax <- function(genus, epithet, type = "no", sub = NULL) {
  out <- list(genus = genus, species = paste(genus, epithet),
              `species epithet` = epithet, `type strain` = type)
  if (!is.null(sub)) out[["subspecies epithet"]] <- sub
  out
}

.general <- function(id, taxid, keywords = NULL) {
  out <- list(`@ref` = 20729L, `BacDive-ID` = id,
              `NCBI tax id` = list(
                list(`NCBI tax id` = taxid, `Matching level` = "species"),
                list(`NCBI tax id` = taxid + 900000L, `Matching level` = "strain")
              ))
  if (!is.null(keywords)) out$keywords <- keywords
  out
}

.util <- function(...) {
  ## each argument: c(metabolite, activity, kind, chebi)
  lapply(list(...), function(v) {
    list(`@ref` = 119508L, `Chebi-ID` = v[[4]], metabolite = v[[1]],
         `utilization activity` = v[[2]], `kind of utilization tested` = v[[3]])
  })
}

.enz <- function(...) {
  lapply(list(...), function(v) {
    list(`@ref` = 20729L, value = v[[1]], activity = v[[2]], ec = v[[3]])
  })
}

fixture_records <- function() {
  recs <- list()

  ## ---- Bacteroides uniformis: three strains, two call it anaerobe, one
  ## ---- microaerophile. Majority + priority must return anaerobe.
  recs[["1001"]] <- list(
    General = .general(1001L, 820L),
    `Name and taxonomic classification` = .tax("Bacteroides", "uniformis", "yes"),
    Morphology = list(`cell morphology` = list(
      `@ref` = 119508L, `gram stain` = "negative",
      `cell shape` = "rod-shaped", motility = "no")),
    `Physiology and metabolism` = list(
      `oxygen tolerance` = list(`@ref` = 68369L, `oxygen tolerance` = "anaerobe"),
      `spore formation` = list(`@ref` = 68369L, `spore formation` = "no"),
      `metabolite utilization` = .util(
        c("D-glucose", "+", "carbon source", "17634"),
        c("maltose",   "+", "carbon source", "17306"),
        c("nitrate",   "+", "builds gas from", "17632"),
        c("D-xylose",  "-", "carbon source", "53455")),
      enzymes = .enz(c("beta-galactosidase", "+", "3.2.1.23"),
                     c("catalase", "-", "1.11.1.6")),
      `API zym` = list(`@ref` = 119508L, Control = "-",
                       `Alkaline phosphatase` = "+", `Esterase Lipase` = "+",
                       `Leucine arylamidase` = "+", Trypsin = "-")
    ),
    `Isolation, sampling and environmental information` = list(
      `isolation source categories` = list(
        list(Cat1 = "#Host", Cat2 = "#Human", Cat3 = "#Intestine"))),
    `Interaction and safety` = list(
      `risk assessment` = list(`@ref` = 20729L, `biosafety level` = "1"))
  )

  recs[["1002"]] <- list(
    General = .general(1002L, 820L),
    `Name and taxonomic classification` = .tax("Bacteroides", "uniformis"),
    `Physiology and metabolism` = list(
      ## Same species, conflicting call -- and the one v7 let win.
      `oxygen tolerance` = list(`@ref` = 68371L, `oxygen tolerance` = "microaerophile"),
      `metabolite utilization` = .util(
        c("glucose", "+", "carbon source", "17234"),
        c("maltose", "+", "carbon source", "17306")),
      enzymes = .enz(c("esterase lipase (C 8)", "+", NA))
    )
  )

  recs[["1003"]] <- list(
    General = .general(1003L, 820L),
    `Name and taxonomic classification` = .tax("Bacteroides", "uniformis"),
    `Physiology and metabolism` = list(
      `oxygen tolerance` = list(`@ref` = 68372L, `oxygen tolerance` = "obligate anaerobe"),
      `metabolite utilization` = .util(c("D-xylose", "-", "carbon source", "53455"))
    )
  )

  ## ---- Faecalibacterium prausnitzii: the type strain has no physiology,
  ## ---- a non-type strain has plenty. v7 saw nothing at all.
  recs[["2001"]] <- list(
    General = .general(2001L, 853L),
    `Name and taxonomic classification` = .tax("Faecalibacterium", "prausnitzii", "yes")
  )
  recs[["2002"]] <- list(
    General = .general(2002L, 853L),
    `Name and taxonomic classification` = .tax("Faecalibacterium", "prausnitzii"),
    Morphology = list(`cell morphology` = list(
      `@ref` = 68369L, `gram stain` = "negative",
      `cell shape` = "rod-shaped", motility = "no")),
    `Physiology and metabolism` = list(
      `oxygen tolerance` = list(`@ref` = 68369L, `oxygen tolerance` = "obligate anaerobe"),
      `spore formation` = list(`@ref` = 68369L, `spore formation` = "no"),
      `metabolite production` = list(
        list(`@ref` = 68369L, `Chebi-ID` = 17968L, metabolite = "butyrate", production = "yes"),
        list(`@ref` = 68369L, `Chebi-ID` = 35581L, metabolite = "indole", production = "no")),
      `metabolite utilization` = .util(c("D-glucose", "+", "carbon source", "17634"))
    )
  )

  ## ---- Escherichia coli: facultative anaerobe, motile, rich panel.
  recs[["3001"]] <- list(
    General = .general(3001L, 562L, keywords = c("Gram-negative", "motile", "rod-shaped")),
    `Name and taxonomic classification` = .tax("Escherichia", "coli", "yes"),
    Morphology = list(`cell morphology` = list(
      `@ref` = 20729L, `gram stain` = "negative",
      `cell shape` = "rod-shaped", motility = "yes")),
    `Physiology and metabolism` = list(
      `oxygen tolerance` = list(`@ref` = 20729L,
                                `oxygen tolerance` = "facultative anaerobe"),
      `spore formation` = list(`@ref` = 20729L, `spore formation` = "no"),
      `metabolite utilization` = .util(
        c("D-glucose", "+", "carbon source", "17634"),
        c("lactose",   "+", "carbon source", "17716"),
        c("citrate",   "-", "carbon source", "16947")),
      `metabolite tests` = list(`@ref` = 68369L, `Chebi-ID` = 35581L,
                                metabolite = "indole", `indole test` = "+"),
      `metabolite production` = list(
        list(`@ref` = 68369L, `Chebi-ID` = 35581L, metabolite = "indole", production = "yes")),
      enzymes = .enz(c("beta-galactosidase", "+", "3.2.1.23"),
                     c("catalase", "+", "1.11.1.6"),
                     c("cytochrome-c oxidase", "-", "1.9.3.1")),
      `antibiotic resistance` = list(`@ref` = 119508L, metabolite = "ampicillin",
                                     `is sensitive` = "no", `is resistant` = "yes")
    ),
    `Interaction and safety` = list(
      `risk assessment` = list(`@ref` = 20729L, `pathogenicity human` = "yes",
                               `biosafety level` = "2"))
  )

  ## ---- Akkermansia muciniphila: the oxygen tolerance block is a LIST of two
  ## ---- entries -- one laboratory observation and one genome-based prediction
  ## ---- carrying a confidence score. Both must be kept, in separate families.
  ## ---- Spore formation is prediction-ONLY, which proves a prediction never
  ## ---- leaks into a measured set.
  recs[["4001"]] <- list(
    General = .general(4001L, 239935L),
    `Name and taxonomic classification` = .tax("Akkermansia", "muciniphila", "yes"),
    Morphology = list(`cell morphology` = list(
      `@ref` = 68369L, `gram stain` = "negative", `cell shape` = "oval-shaped",
      motility = "no")),
    `Physiology and metabolism` = list(
      `oxygen tolerance` = list(
        list(`@ref` = 68369L, `oxygen tolerance` = "anaerobe"),
        list(`@ref` = 69480L, `oxygen tolerance` = "anaerobe", confidence = "94.2")),
      `spore formation` = list(`@ref` = 69480L, `spore formation` = "no",
                               confidence = "88.1"),
      observation = list(`@ref` = 68369L, observation = "mucin degradation")
    )
  )

  ## ---- A record returned by a genus query that is not the species we asked
  ## ---- for. Proves the species guard: it must be kept under its own species,
  ## ---- never folded into the queried one.
  recs[["5001"]] <- list(
    General = .general(5001L, 246787L),
    `Name and taxonomic classification` = .tax("Bacteroides", "cellulosilyticus", "yes"),
    `Physiology and metabolism` = list(
      `oxygen tolerance` = list(`@ref` = 68369L, `oxygen tolerance` = "anaerobe"))
  )

  recs
}

## Write the fixture as if it were a harvested genus cache, so stage 02 can be
## run against it unchanged.
write_fixture_cache <- function(dir = file.path("tests", "fixtures", "cache")) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  recs <- fixture_records()
  by_genus <- split(recs, vapply(recs, function(r) {
    r[["Name and taxonomic classification"]]$genus
  }, character(1)))
  for (g in names(by_genus)) {
    saveRDS(list(genus = g, harvested_at = as.POSIXct("2026-09-16 00:00:00", tz = "UTC"),
                 n_records = length(by_genus[[g]]), records = by_genus[[g]]),
            file.path(dir, paste0(gsub("[^A-Za-z0-9]+", "_", g), ".rds")))
  }
  invisible(names(by_genus))
}
