## lib/sets.R -- species-level trait calls -> named taxon sets for TaxSEA_db.

## Collapse analyte spellings across every family that shares a set_infix.
## `enzymes` and the `API zym` panel both emit BacDive_Enzyme_* names, so if this
## ran per family (as it did at first) "Esterase Lipase" and "esterase lipase"
## would still become two sets -- the exact split being fixed.
collapse_analytes <- function(traits) {
  key <- ifelse(is.na(traits$set_infix) | !nzchar(traits$set_infix),
                paste0("#", traits$family), traits$set_infix)
  for (k in unique(key)) {
    i <- key == k
    traits$analyte[i] <- collapse_case(traits$analyte[i], traits$curated[i])
  }
  traits
}

## Drop polarity sets that cover almost the whole tested universe for their
## analyte. Such a set cannot discriminate: "does not ferment L-xylose" is true of
## 98% of the species ever tested for it, so it carries no information and is
## strongly correlated with every other rare-sugar negative set.
##
## Deliberately one-sided: a SMALL set is a specific hypothesis and is kept, no
## matter how small. Only near-universal sets are removed.
drop_uninformative <- function(traits, spec) {
  key <- paste(traits$set_infix, traits$analyte, traits$context, sep = "\r")
  universe <- table(key)
  n_pos <- table(key[traits$call == "pos"])
  n_neg <- table(key[traits$call == "neg"])

  u <- as.numeric(universe[key])
  share <- ifelse(traits$call == "pos",
                  as.numeric(n_pos[key]) / u,
                  as.numeric(n_neg[key]) / u)
  share[is.na(share)] <- 0

  ## Only assay/panel families, where pos and neg are two outcomes of the SAME
  ## analyte and so genuinely partition a universe.
  ##
  ## It must NOT touch categorical families: there the polarity is encoded in the
  ## analyte itself (Motility_yes vs Motility_no), so each set is 100% of its own
  ## universe by construction and gating on emit_negative alone deleted the
  ## entire Motility and Spore families.
  assay_fams <- spec$family[spec$kind %in% c("assay", "panel")]
  ## A share only means something once a reasonable number of species were
  ## tested: with a universe of 1 every set is 100% of it by definition, so an
  ## unguarded filter deletes small specific traits instead of near-universal ones.
  eligible <- traits$family %in% assay_fams & u >= PIPE$min_universe_for_filter
  drop <- eligible & share > PIPE$max_universe_fraction
  if (any(drop)) {
    lost <- unique(paste(traits$set_infix[drop], traits$analyte[drop],
                         traits$context[drop], traits$call[drop]))
    log_msg(sprintf("information filter: dropped %d species-rows across %d near-universal polarity sets (>%.0f%% of their tested universe)",
                    sum(drop), length(lost), 100 * PIPE$max_universe_fraction))
  }
  traits[!drop, , drop = FALSE]
}

## Name every set. Negative-polarity sets are only emitted for families that
## declare emit_negative, because for most traits "not reported" and "negative"
## are not the same thing and only the assay families can tell them apart.
name_sets <- function(traits, spec) {
  em <- stats::setNames(tolower(spec$emit_negative) %in% c("yes", "true", "1"), spec$family)
  traits$emit_negative <- unname(em[traits$family])
  traits$emit_negative[is.na(traits$emit_negative)] <- FALSE

  ## A categorical family encodes its polarity in the value itself
  ## (Motility_yes / Motility_no), so its rows are all "pos" and need no suffix.
  is_cat <- traits$family %in% spec$family[spec$kind == "categorical"]

  keep <- traits$call == "pos" | (traits$emit_negative & !is_cat)
  traits <- traits[keep, , drop = FALSE]
  if (!nrow(traits)) return(traits)

  polarity <- ifelse(traits$call == "neg", "negative", "")
  ctx <- ifelse(is.na(traits$context), "", traits$context)
  suffix <- ifelse(nzchar(ctx) & nzchar(polarity), paste(ctx, polarity, sep = "_"),
                   ifelse(nzchar(ctx), ctx, polarity))

  traits$set_name <- vapply(seq_len(nrow(traits)), function(i) {
    make_set_name(traits$set_infix[i], traits$analyte[i],
                  predicted = traits$predicted[i],
                  polarity = suffix[i])
  }, character(1))
  traits
}

## Higher-order sets formed by unioning canonical values within a family, defined
## in data/curation/derived_unions.csv. Additive: the constituent sets are kept.
add_derived_unions <- function(traits) {
  path <- file.path(PIPE$curation_dir, "derived_unions.csv")
  if (!file.exists(path)) return(traits)
  u <- load_curation("derived_unions.csv", c("set_name", "family", "member_analytes"))
  if (!nrow(u)) return(traits)

  extra <- list()
  for (i in seq_len(nrow(u))) {
    members <- trimws(strsplit(u$member_analytes[i], "|", fixed = TRUE)[[1]])
    src <- traits[traits$family == u$family[i] & traits$analyte %in% members &
                    traits$call == "pos", , drop = FALSE]
    if (!nrow(src)) next
    ## One row per species: a species in both constituent sets joins the union once.
    src <- src[order(src$species, -src$n_support), , drop = FALSE]
    src <- src[!duplicated(src$species), , drop = FALSE]
    src$set_name  <- u$set_name[i]
    src$family    <- paste0(u$family[i], "_union")
    src$set_infix <- NA_character_
    src$analyte   <- u$set_name[i]
    src$context   <- ""
    extra[[length(extra) + 1L]] <- src
  }
  if (!length(extra)) return(traits)
  out <- do.call(rbind, extra)
  log_msg(sprintf("derived unions: %d set(s), %d species-rows added",
                  length(unique(out$set_name)), nrow(out)))
  rbind(traits, out)
}

## Build the named list of NCBI taxon ids that TaxSEA consumes.
##
## Two things v7 got wrong here:
##  * split() was keyed on set_id alone, so with the 4-column table every oxygen
##    category collapsed into one "oxygen_tolerance" set containing every species
##    regardless of call, and gram_stain held both polarities;
##  * ids were never de-duplicated, so the 13 sets in which Lactobacillus
##    fructivorans and L. homohiochii both map to taxid 1614 double-counted it.
build_taxon_sets <- function(traits) {
  t2 <- traits[!is.na(traits$ncbi_taxid) & nzchar(traits$ncbi_taxid), , drop = FALSE]
  dropped <- setdiff(unique(traits$species), unique(t2$species))
  sets <- split(t2$ncbi_taxid, t2$set_name)
  sets <- lapply(sets, function(x) unique(x[!is.na(x)]))
  sets <- sets[order(names(sets))]
  attr(sets, "species_without_taxid") <- dropped
  sets
}

## A per-set provenance row, shipped alongside the sets. This is what lets a
## reader see that a set is lab-measured rather than predicted, how many strains
## stand behind it, and how many species were actually tested -- the denominator
## v7 could not express.
set_provenance <- function(traits, sets) {
  t2 <- traits[!is.na(traits$ncbi_taxid) & nzchar(traits$ncbi_taxid), , drop = FALSE]
  sp <- split(seq_len(nrow(t2)), t2$set_name)
  rows <- lapply(names(sp), function(nm) {
    i <- sp[[nm]]
    data.frame(
      set_name         = nm,
      family           = paste(sort(unique(t2$family[i])), collapse = "|"),
      section_source   = paste(sort(unique(t2$set_infix[i])), collapse = "|"),
      analyte          = t2$analyte[i][[1]],
      polarity         = t2$call[i][[1]],
      evidence         = if (all(t2$predicted[i])) "predicted" else "measured",
      n_species        = length(unique(t2$species[i])),
      n_taxa           = length(sets[[nm]] %||% character(0)),
      n_strains_total  = sum(t2$n_tested[i], na.rm = TRUE),
      median_n_tested  = stats::median(t2$n_tested[i], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(-out$n_taxa, out$set_name), , drop = FALSE]
}

## The "tested" universe per set family: the species BacDive actually assayed for
## that analyte, positive or negative. The correct enrichment background for
## "utilises glucose" is the species tested for glucose, not every species in the
## database -- without this, annotation depth masquerades as biology.
build_tested_universe <- function(traits) {
  t2 <- traits[!is.na(traits$ncbi_taxid) & nzchar(traits$ncbi_taxid) &
                 !is.na(traits$n_tested) & traits$n_tested > 0, , drop = FALSE]
  key <- paste(t2$set_infix, t2$analyte, t2$context, sep = "|")
  u <- split(t2$ncbi_taxid, key)
  lapply(u, function(x) unique(x[!is.na(x)]))
}
