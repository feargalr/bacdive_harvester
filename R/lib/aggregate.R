## lib/aggregate.R -- strain-level observations to species-level trait calls.
##
## The v7 pipeline collapsed this problem by taking one strain record per species
## (the first type strain retrieve() happened to return) and reading it. That
## discarded every other strain BacDive had already sent, made the answer depend
## on API record ordering, and left no way to express "tested and negative" as
## distinct from "never tested".
##
## Here every strain contributes, and each call carries its evidence:
##   n_support / n_strains  for categorical traits
##   n_pos / n_neg / n_tested for assays
## so a downstream user can filter on evidence and so the annotation bias is
## visible rather than implicit.

## Reshape the long table to one row per (strain, subsection entry) with the
## requested fields as columns. Used for assays (analyte + result + context live
## in the same entry) and for categorical traits (value + confidence likewise).
pivot_entries <- function(d, fields) {
  if (!nrow(d)) return(NULL)
  key <- paste(d$bacdive_id, d$entry, sep = "\r")
  ukey <- unique(key)
  out <- data.frame(
    bacdive_id = sub("\r.*$", "", ukey),
    entry      = as.integer(sub("^.*\r", "", ukey)),
    stringsAsFactors = FALSE
  )
  for (f in fields) {
    if (is.na(f) || !nzchar(f)) next
    sel <- d$field == f
    out[[f]] <- d$value[sel][match(ukey, key[sel])]
  }
  ## Carry the reference id through for provenance.
  sel <- !is.na(d$ref)
  out$ref <- d$ref[sel][match(ukey, key[sel])]
  out
}

## Rows of the flat long table belonging to one trait-spec family.
spec_rows <- function(strain_long, sp) {
  ## subsection is "" (never NA) for fields sitting directly on a section, so
  ## these comparisons cannot produce NA.
  keep <- strain_long$section == sp$section & strain_long$subsection == sp$subsection
  strain_long[keep, , drop = FALSE]
}

classify_result <- function(x, pos_pattern, neg_pattern) {
  x <- trimws(as.character(x))
  out <- rep(NA_character_, length(x))
  if (!is.na(pos_pattern) && nzchar(pos_pattern)) {
    out[grepl(pos_pattern, x, ignore.case = TRUE, perl = TRUE)] <- "pos"
  }
  if (!is.na(neg_pattern) && nzchar(neg_pattern)) {
    out[is.na(out) & grepl(neg_pattern, x, ignore.case = TRUE, perl = TRUE)] <- "neg"
  }
  out
}

## A confidence score on an entry means the value came from BacDive's
## genome-based prediction rather than a laboratory observation. v7 read the
## value and dropped the confidence, merging the two kinds of evidence.
is_predicted_entry <- function(tab) {
  if (is.null(tab$confidence)) return(rep(FALSE, nrow(tab)))
  !is.na(tab$confidence) & nzchar(tab$confidence)
}

## ---- categorical ------------------------------------------------------------
## One value per strain. Species call = most-supported value; ties broken by the
## curated priority (lower wins), which is how "anaerobe; microaerophile" lands
## on anaerobe instead of microaerophile.
aggregate_categorical <- function(strain_long, meta, sp) {
  d <- spec_rows(strain_long, sp)
  if (!nrow(d)) return(NULL)
  tab <- pivot_entries(d, unique(c(sp$result_field, "confidence")))
  if (is.null(tab) || is.null(tab[[sp$result_field]])) return(NULL)

  vmap <- load_value_map(sp$alias_table)
  norm <- normalise_value(tab[[sp$result_field]], vmap)
  tab$canonical <- norm$canonical
  tab$priority  <- norm$priority
  tab$dominant  <- norm$dominant
  tab$predicted <- is_predicted_entry(tab)
  tab <- tab[!is.na(tab$canonical) & nzchar(tab$canonical), , drop = FALSE]
  if (!nrow(tab)) return(NULL)

  tab <- merge(tab, meta[, c("bacdive_id", "species")], by = "bacdive_id", all.x = TRUE)
  tab <- tab[!is.na(tab$species), , drop = FALSE]
  if (!nrow(tab)) return(NULL)

  ## Measured and predicted evidence are resolved separately so a predicted call
  ## can never outvote a laboratory one.
  do.call(rbind, lapply(c(FALSE, TRUE), function(pred) {
    t2 <- tab[tab$predicted == pred, , drop = FALSE]
    if (!nrow(t2)) return(NULL)
    ## One vote per (strain, value): a strain reporting the same value twice
    ## under two references must not count twice.
    t2 <- t2[!duplicated(t2[, c("species", "bacdive_id", "canonical")]), , drop = FALSE]
    n_strains <- tapply(t2$bacdive_id, t2$species, function(z) length(unique(z)))
    ## Count support per canonical value only. Priority must NOT be part of the
    ## grouping key: "anaerobe" and "obligate anaerobe" both canonicalise to
    ## "anaerobe" but carry priorities 2 and 1, so grouping on it would split two
    ## agreeing strains into two groups of one and could flip the winner.
    agg <- aggregate(list(n_support = t2$bacdive_id),
                     by = list(species = t2$species, canonical = t2$canonical),
                     FUN = function(z) length(unique(z)))
    ## Priority of a canonical value is the strongest (lowest) seen for it.
    prio <- tapply(t2$priority, t2$canonical, min, na.rm = TRUE)
    agg$priority <- as.integer(prio[agg$canonical])

    ## A `dominant` value wins on presence alone, ahead of majority voting. This
    ## is what makes oxygen tolerance behave as a capability rather than a vote:
    ## a facultative anaerobe grown only aerobically is recorded as "aerobe", so
    ## counting strains made Enterobacteriaceae come out aerobic.
    dom <- tapply(t2$dominant, t2$canonical, any)
    agg$dominant <- unname(dom[agg$canonical])
    agg$dominant[is.na(agg$dominant)] <- FALSE
    ## Dominance needs support behind it, or a single stray annotation flips the
    ## species (P. aeruginosa: 1 facultative call among 142 strains).
    n_sp <- as.integer(n_strains[agg$species])
    agg$dominant <- agg$dominant &
      agg$n_support >= PIPE$dominant_min_n &
      (agg$n_support / pmax(n_sp, 1L)) >= PIPE$dominant_min_frac

    ## Winner per species: dominant first, then most strains, then curated
    ## priority, then name for determinism.
    ord <- order(agg$species, -agg$dominant, -agg$n_support, agg$priority,
                 agg$canonical)
    agg <- agg[ord, , drop = FALSE]
    win <- agg[!duplicated(agg$species), , drop = FALSE]
    data.frame(
      species   = win$species,
      family    = sp$family,
      set_infix = sp$set_infix,
      analyte   = win$canonical,
      context   = "",
      call      = "pos",
      n_support = win$n_support,
      n_pos     = win$n_support,
      n_neg     = NA_integer_,
      n_tested  = as.integer(n_strains[win$species]),
      predicted = pred,
      curated   = TRUE,   # categorical values always come from a curated value map
      context_raw = NA_character_,
      stringsAsFactors = FALSE
    )
  }))
}

## ---- assay and panel --------------------------------------------------------
## Binary +/- readouts. Both polarities are counted, which is what gives us real
## negatives and a defensible "tested" universe per set.
aggregate_assay <- function(strain_long, meta, sp, panel = FALSE) {
  d <- spec_rows(strain_long, sp)
  if (!nrow(d)) return(NULL)

  if (panel) {
    ## A flat panel: every field is an analyte, its value is the result.
    tab <- data.frame(bacdive_id = d$bacdive_id, entry = d$entry,
                      analyte_raw = d$field, result = d$value,
                      context_raw = NA_character_, id = NA_character_,
                      confidence = NA_character_, stringsAsFactors = FALSE)
    if (!is.na(sp$exclude_analyte_pattern) && nzchar(sp$exclude_analyte_pattern)) {
      tab <- tab[!grepl(sp$exclude_analyte_pattern, tab$analyte_raw, perl = TRUE), , drop = FALSE]
    }
  } else {
    want <- unique(c(sp$analyte_field, sp$result_field, sp$context_field,
                     "Chebi-ID", "ec", "confidence"))
    pv <- pivot_entries(d, want)
    if (is.null(pv) || is.null(pv[[sp$analyte_field]]) || is.null(pv[[sp$result_field]])) {
      return(NULL)
    }
    ctx <- if (!is.na(sp$context_field) && nzchar(sp$context_field)) pv[[sp$context_field]] else NULL
    id  <- if (!is.null(pv[["Chebi-ID"]])) pv[["Chebi-ID"]] else pv[["ec"]]
    tab <- data.frame(
      bacdive_id  = pv$bacdive_id, entry = pv$entry,
      analyte_raw = pv[[sp$analyte_field]],
      result      = pv[[sp$result_field]],
      context_raw = if (is.null(ctx)) NA_character_ else ctx,
      id          = if (is.null(id)) NA_character_ else id,
      confidence  = if (is.null(pv$confidence)) NA_character_ else pv$confidence,
      stringsAsFactors = FALSE
    )
    if (!is.na(sp$exclude_analyte_pattern) && nzchar(sp$exclude_analyte_pattern)) {
      tab <- tab[!grepl(sp$exclude_analyte_pattern, tab$analyte_raw, perl = TRUE), , drop = FALSE]
    }
  }

  ## Context include/exclude: one subsection can feed two families.
  ctx_for_filter <- ifelse(is.na(tab$context_raw), "", tab$context_raw)
  if (!is.na(sp$context_include) && nzchar(sp$context_include)) {
    tab <- tab[grepl(sp$context_include, ctx_for_filter, perl = TRUE), , drop = FALSE]
    ctx_for_filter <- ctx_for_filter[grepl(sp$context_include, ctx_for_filter, perl = TRUE)]
  }
  if (!is.na(sp$context_exclude) && nzchar(sp$context_exclude) && nrow(tab)) {
    keep_ctx <- !grepl(sp$context_exclude, ctx_for_filter, perl = TRUE)
    tab <- tab[keep_ctx, , drop = FALSE]
  }
  tab <- tab[!is.na(tab$analyte_raw) & !is.na(tab$result), , drop = FALSE]
  if (!nrow(tab)) return(NULL)

  amap <- load_alias_map(sp$alias_table)
  cmap <- load_alias_map("context_aliases.csv")

  na <- normalise_analyte(tab$analyte_raw, amap)
  tab$analyte <- na$value
  tab$curated <- na$curated
  ## Prefer BacDive's own identifier as the merge key where it exists: immune to
  ## spelling, and it is what ties "D-glucose" to "glucose" when both carry the
  ## same Chebi-ID.
  tab$analyte <- apply_id_bridge(tab$analyte, tab$id, build_id_bridge(tab$analyte, tab$id))
  ## NB: case collapse is deliberately NOT done here. Families that share a
  ## set_infix (enzymes + API zym both produce BacDive_Enzyme_*) must be
  ## collapsed together, which only stage 04 can see. See collapse_analytes().

  tab$context <- if (identical(tolower(sp$context_mode), "suffix")) {
    normalise_context(tab$context_raw, cmap)
  } else {
    rep("", nrow(tab))
  }
  ## Keep the raw BacDive assay context. Several raw contexts collapse into one
  ## TaxSEA feature (growth, assimilation, degradation, oxidation -> Uses_X), so
  ## without this the original granularity would be unrecoverable.
  tab$ctx_raw_keep <- ifelse(is.na(tab$context_raw) | !nzchar(tab$context_raw),
                             "<none>", tab$context_raw)

  tab$cls <- classify_result(tab$result, sp$positive_pattern, sp$negative_pattern)
  tab <- tab[!is.na(tab$cls), , drop = FALSE]
  if (!nrow(tab)) return(NULL)

  tab$predicted <- is_predicted_entry(tab)
  tab <- merge(tab, meta[, c("bacdive_id", "species")], by = "bacdive_id", all.x = TRUE)
  tab <- tab[!is.na(tab$species), , drop = FALSE]
  if (!nrow(tab)) return(NULL)

  ## One vote per strain per analyte: if a strain was tested twice with the same
  ## outcome that is one observation, not two. Conflicting outcomes within a
  ## strain are kept and cancel out in the fraction.
  tab <- tab[!duplicated(tab[, c("species", "bacdive_id", "analyte", "context", "cls")]), ,
             drop = FALSE]

  ## NB: `curated` must NOT be a grouping key. It describes where the label came
  ## from, not what was measured -- grouping on it split one strain reporting
  ## "D-glucose" (alias hit) from another reporting "glucose" (no hit) into two
  ## rows of n_pos=1 instead of one row of n_pos=2.
  curated_by_analyte <- tapply(tab$curated, tab$analyte, any)
  grp <- list(species = tab$species, analyte = tab$analyte, context = tab$context,
              predicted = tab$predicted)
  n_pos <- aggregate(list(v = tab$cls == "pos"), by = grp, FUN = sum)
  n_neg <- aggregate(list(v = tab$cls == "neg"), by = grp, FUN = sum)
  ## Which raw contexts contributed to each collapsed group.
  raw_ctx <- aggregate(list(v = tab$ctx_raw_keep), by = grp,
                       FUN = function(z) paste(sort(unique(z)), collapse = "|"))
  key_cols <- c("species", "analyte", "context", "predicted")
  agg <- merge(n_pos, n_neg, by = key_cols, suffixes = c("_pos", "_neg"))
  agg <- merge(agg, raw_ctx, by = key_cols)
  names(agg)[names(agg) == "v"] <- "context_raw"
  agg$n_pos    <- as.integer(agg$v_pos)
  agg$n_neg    <- as.integer(agg$v_neg)
  agg$n_tested <- agg$n_pos + agg$n_neg

  frac <- ifelse(agg$n_tested > 0, agg$n_pos / agg$n_tested, NA_real_)
  agg$call <- ifelse(is.na(frac), NA_character_,
                     ifelse(frac >= PIPE$min_pos_fraction, "pos", "neg"))
  agg <- agg[!is.na(agg$call) & agg$n_tested >= PIPE$min_strains_tested, , drop = FALSE]
  if (!nrow(agg)) return(NULL)

  data.frame(
    species   = agg$species,
    family    = sp$family,
    set_infix = sp$set_infix,
    analyte   = agg$analyte,
    context   = agg$context,
    call      = agg$call,
    n_support = ifelse(agg$call == "pos", agg$n_pos, agg$n_neg),
    n_pos     = agg$n_pos,
    n_neg     = agg$n_neg,
    n_tested  = agg$n_tested,
    predicted = agg$predicted,
    curated   = unname(curated_by_analyte[agg$analyte]),
    context_raw = agg$context_raw,
    stringsAsFactors = FALSE
  )
}

## ---- label ------------------------------------------------------------------
## Presence-only descriptors (isolation source category, observation, nutrition
## type...). There is no meaningful negative: BacDive not recording "mucin
## degradation" for a strain is absence of evidence, not evidence of absence.
aggregate_label <- function(strain_long, meta, sp) {
  d <- spec_rows(strain_long, sp)
  d <- d[d$field == sp$result_field, , drop = FALSE]
  if (!nrow(d)) return(NULL)

  d$analyte <- canon_whitespace(clean_vocab(d$value))  # case collapse in stage 04
  d <- d[!is.na(d$analyte) & nzchar(d$analyte), , drop = FALSE]

  ## Same semantic allowlist facility as the typed_label path.
  if (isTRUE(tolower(sp$restrict_to_alias %||% "no") %in% c("yes", "true", "1"))) {
    allow <- load_allowlist(sp$alias_table)
    if (is.null(allow)) stop("restrict_to_alias needs an alias_table: ", sp$family,
                             call. = FALSE)
    cls <- apply_allowlist(d$analyte, allow)
    d <- d[!is.na(cls), , drop = FALSE]
    cls <- cls[!is.na(cls)]
    if (!nrow(d)) return(NULL)
    d$analyte <- cls
  }
  d <- merge(d, meta[, c("bacdive_id", "species")], by = "bacdive_id", all.x = TRUE)
  d <- d[!is.na(d$species), , drop = FALSE]
  if (!nrow(d)) return(NULL)

  d <- d[!duplicated(d[, c("species", "bacdive_id", "analyte")]), , drop = FALSE]
  n_strains <- tapply(d$bacdive_id, d$species, function(z) length(unique(z)))
  agg <- aggregate(list(n_support = d$bacdive_id),
                   by = list(species = d$species, analyte = d$analyte),
                   FUN = function(z) length(unique(z)))
  data.frame(
    species   = agg$species,
    family    = sp$family,
    set_infix = sp$set_infix,
    analyte   = agg$analyte,
    context   = "",
    call      = "pos",
    n_support = agg$n_support,
    n_pos     = agg$n_support,
    n_neg     = NA_integer_,
    n_tested  = as.integer(n_strains[agg$species]),
    predicted = FALSE,
    curated   = FALSE,
    context_raw = NA_character_,
    stringsAsFactors = FALSE
  )
}

## ---- driver -----------------------------------------------------------------
aggregate_all <- function(strain_long, meta, spec) {
  out <- vector("list", nrow(spec))
  for (i in seq_len(nrow(spec))) {
    sp <- as.list(spec[i, ])
    res <- switch(sp$kind,
      categorical = aggregate_categorical(strain_long, meta, sp),
      assay       = aggregate_assay(strain_long, meta, sp, panel = FALSE),
      panel       = aggregate_assay(strain_long, meta, sp, panel = TRUE),
      label       = aggregate_label(strain_long, meta, sp),
      typed_label = aggregate_typed_label(strain_long, meta, sp),
      numeric     = aggregate_numeric(strain_long, meta, sp),
      stop("unknown kind: ", sp$kind)
    )
    n <- if (is.null(res)) 0L else nrow(res)
    log_msg(sprintf("  %-16s %-11s %5d species-trait calls", sp$family, sp$kind, n))
    out[[i]] <- res
  }
  res <- do.call(rbind, out)
  if (is.null(res)) return(NULL)
  rownames(res) <- NULL
  res
}

## ---- derived: numeric traits ------------------------------------------------
## Generalises what was originally a salinity-only path. Several BacDive blocks
## have the same shape: one entry per tested value, with a `type` saying whether
## the value is a growth point, an optimum, a minimum or a maximum. Treated
## categorically they manufacture sets that are not distinct traits --
## 633 distinct NaCl concentrations, 767 temperatures, 456 pH values.
##
## Instead: parse every value to a number, reduce to one statistic per species
## (the highest or lowest supporting growth, or the median for a measured
## quantity), and emit conventional thresholds and classes. Configured by
## data/curation/numeric_traits.csv and numeric_classes.csv.
aggregate_numeric <- function(strain_long, meta, sp) {
  cfg <- load_curation("numeric_traits.csv",
                       c("family", "stat", "relations", "thresholds", "threshold_fmt"))
  cfg <- cfg[cfg$family == sp$family, , drop = FALSE]
  if (!nrow(cfg)) stop("numeric_traits.csv has no row for family ", sp$family,
                       call. = FALSE)
  cfg <- as.list(cfg[1, ])

  cls <- load_curation("numeric_classes.csv", c("family", "lower", "upper", "label"))
  cls <- cls[cls$family == sp$family, , drop = FALSE]

  d <- spec_rows(strain_long, sp)
  if (!nrow(d)) return(NULL)
  want <- unique(c(sp$analyte_field, sp$result_field, sp$context_field, "salt"))
  tab <- pivot_entries(d, want[!is.na(want) & nzchar(want)])
  if (is.null(tab) || is.null(tab[[sp$analyte_field]])) return(NULL)

  ## Salt identity only applies to halophily; NaCl alone, because marine salts,
  ## MgCl2 and KCl are different physiology and far too sparse to support traits.
  if (!is.null(tab$salt)) {
    ok <- is.na(tab$salt) | grepl("^NaCl$", trimws(tab$salt), ignore.case = TRUE)
    tab <- tab[ok, , drop = FALSE]
  }

  ## Admissible `type` values. "optimum" is a different quantity from a maximum
  ## and must never be read as one.
  if (!is.na(cfg$relations) && nzchar(cfg$relations) &&
      !is.na(sp$context_field) && nzchar(sp$context_field)) {
    rel <- tolower(trimws(tab[[sp$context_field]] %||% rep(NA_character_, nrow(tab))))
    allowed <- strsplit(tolower(cfg$relations), "|", fixed = TRUE)[[1]]
    tab <- tab[is.na(rel) | rel %in% allowed, , drop = FALSE]
  }
  if (!nrow(tab)) return(NULL)

  ## A growth/ability column, where the block has one.
  if (!is.na(sp$result_field) && nzchar(sp$result_field) &&
      !is.null(tab[[sp$result_field]])) {
    grew <- grepl("^(positive|yes|\\+)$", trimws(tab[[sp$result_field]]),
                  ignore.case = TRUE)
    tab <- tab[grew, , drop = FALSE]
  }
  if (!nrow(tab)) return(NULL)

  rng <- parse_range(tab[[sp$analyte_field]])
  ## A max-statistic takes the top of a range, a min-statistic the bottom.
  tab$val <- if (identical(cfg$stat, "min")) rng$low else rng$high
  tab <- tab[!is.na(tab$val), , drop = FALSE]
  if (!nrow(tab)) return(NULL)

  tab <- merge(tab, meta[, c("bacdive_id", "species")], by = "bacdive_id", all.x = TRUE)
  tab <- tab[!is.na(tab$species), , drop = FALSE]
  if (!nrow(tab)) return(NULL)

  fn <- switch(cfg$stat, max = max, min = min, median = stats::median,
               stop("unknown stat: ", cfg$stat))
  stat_by_sp <- tapply(tab$val, tab$species, function(z) fn(z, na.rm = TRUE))
  ns <- tapply(tab$bacdive_id, tab$species, function(z) length(unique(z)))
  spn <- names(stat_by_sp)

  mk <- function(species, analyte) {
    data.frame(species = species, family = sp$family, set_infix = sp$set_infix,
               analyte = analyte, context = "", call = "pos",
               n_support = as.integer(ns[species]), n_pos = as.integer(ns[species]),
               n_neg = NA_integer_, n_tested = as.integer(ns[species]),
               predicted = FALSE, curated = TRUE, context_raw = NA_character_,
               stringsAsFactors = FALSE)
  }

  rows <- list()
  ## Cumulative thresholds: a species growing at 65 C also satisfies >=45 and >=55.
  if (!is.na(cfg$thresholds) && nzchar(cfg$thresholds)) {
    thr <- as.numeric(strsplit(cfg$thresholds, "|", fixed = TRUE)[[1]])
    for (t in thr[!is.na(t)]) {
      hit <- spn[!is.na(stat_by_sp[spn]) &
                   if (identical(cfg$stat, "min")) stat_by_sp[spn] <= t
                   else stat_by_sp[spn] >= t]
      if (length(hit)) rows[[length(rows) + 1L]] <-
        mk(hit, sprintf(cfg$threshold_fmt, t))
    }
  }
  ## One class per species, from the bounds table.
  if (nrow(cls)) {
    lo <- as.numeric(cls$lower); hi <- as.numeric(cls$upper)
    lab <- rep(NA_character_, length(spn))
    v <- as.numeric(stat_by_sp[spn])
    for (k in seq_len(nrow(cls))) {
      inb <- !is.na(v) & v >= lo[k] & v < hi[k] & is.na(lab)
      lab[inb] <- cls$label[k]
    }
    ok <- !is.na(lab)
    if (any(ok)) rows[[length(rows) + 1L]] <- mk(spn[ok], lab[ok])
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

## ---- typed labels -----------------------------------------------------------
## For fields whose values are "<concept>: <value>" or comma-separated lists.
## "quinones: MK-9(H6), MK-9(H8)" becomes two rows, MK-9(H6) and MK-9(H8), so the
## same quinone is comparable across strains instead of every combination (and
## every ordering of a combination) becoming its own set.
aggregate_typed_label <- function(strain_long, meta, sp) {
  d <- spec_rows(strain_long, sp)
  d <- d[d$field == sp$result_field, , drop = FALSE]
  if (!nrow(d)) return(NULL)

  tp <- split_typed_prefix(d$value)
  mode <- if (identical(tolower(sp$context_mode), "explode")) "explode" else "sort"

  if (mode == "explode") {
    comps <- split_components(tp$value, "explode")
    n <- vapply(comps, length, integer(1))
    d2 <- d[rep(seq_len(nrow(d)), n), , drop = FALSE]
    d2$analyte <- unlist(comps, use.names = FALSE)
    d2$concept <- rep(tp$prefix, n)
  } else {
    d2 <- d
    d2$analyte <- split_components(tp$value, "sort")
    d2$concept <- tp$prefix
  }

  ## The concept prefix becomes part of the set name, so "quinones" and
  ## "fatty acids" cannot collide on a shared component name.
  d2$analyte <- ifelse(is.na(d2$concept) | !nzchar(d2$concept), d2$analyte,
                       paste(d2$concept, d2$analyte, sep = "_"))
  d2 <- d2[!is.na(d2$analyte) & nzchar(d2$analyte), , drop = FALSE]
  if (!nrow(d2)) return(NULL)

  ## A semantic allowlist: map raw analytes onto curated classes and discard
  ## anything not listed. The compound-production vocabulary is 1,887 analytes of
  ## which 80% are single-species, largely named one-off natural products
  ## ("antibiotic A-195", "vernamycin A") and fragmented patent citations, so a
  ## frequency threshold would keep obscure compounds that happened to be studied
  ## three times. Only named classes survive.
  if (isTRUE(tolower(sp$restrict_to_alias %||% "no") %in% c("yes", "true", "1"))) {
    allow <- load_allowlist(sp$alias_table)
    if (is.null(allow)) stop("restrict_to_alias needs an alias_table: ", sp$family,
                             call. = FALSE)
    cls <- apply_allowlist(d2$analyte, allow)
    d2 <- d2[!is.na(cls), , drop = FALSE]
    cls <- cls[!is.na(cls)]
    if (!nrow(d2)) return(NULL)
    d2$analyte <- cls
  }

  d2 <- merge(d2, meta[, c("bacdive_id", "species")], by = "bacdive_id", all.x = TRUE)
  d2 <- d2[!is.na(d2$species), , drop = FALSE]
  if (!nrow(d2)) return(NULL)

  d2 <- d2[!duplicated(d2[, c("species", "bacdive_id", "analyte")]), , drop = FALSE]
  n_strains <- tapply(d2$bacdive_id, d2$species, function(z) length(unique(z)))
  agg <- aggregate(list(n_support = d2$bacdive_id),
                   by = list(species = d2$species, analyte = d2$analyte),
                   FUN = function(z) length(unique(z)))
  data.frame(
    species = agg$species, family = sp$family, set_infix = sp$set_infix,
    analyte = agg$analyte, context = "", call = "pos",
    n_support = agg$n_support, n_pos = agg$n_support, n_neg = NA_integer_,
    n_tested = as.integer(n_strains[agg$species]),
    predicted = FALSE, curated = FALSE, context_raw = NA_character_,
    stringsAsFactors = FALSE)
}
