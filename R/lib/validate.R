## lib/validate.R -- checks that would have caught the v7 defects.

## Gold-standard assertions from data/curation/gold_standard.csv.
##
## Four outcomes, not two, because "we got it wrong" and "BacDive is silent" are
## different facts and reporting them as one number misstates quality:
##
##   pass         expectation met
##   fail         BacDive HAS a call in this trait family for this species, and it
##                contradicts the expectation. This is the only real defect.
##   no_data      the species was harvested, but BacDive records nothing in this
##                trait family for it. A coverage ceiling, not an error. Three of
##                the most abundant gut commensals (F. prausnitzii,
##                R. intestinalis, R. bromii) have no oxygen tolerance in BacDive
##                at all, which is precisely why v7 had to borrow those calls from
##                the retired GPT table.
##   not_covered  the species is absent from the harvest entirely.
##
## The family behind a set name is recovered from the traits table rather than
## parsed out of the name, so this stays correct as naming changes.
check_gold_standard <- function(traits, sets, gold) {
  species_seen <- unique(traits$species)
  fam_of_set <- stats::setNames(traits$family, traits$set_name)
  fam_of_set <- fam_of_set[!duplicated(names(fam_of_set))]

  out <- lapply(seq_len(nrow(gold)), function(i) {
    g <- gold[i, ]
    members <- sets[[g$trait_set]] %||% character(0)
    tx <- unique(traits$ncbi_taxid[traits$species == g$species])
    tx <- tx[!is.na(tx)]
    present <- length(tx) > 0L && any(tx %in% members)

    fam <- unname(fam_of_set[g$trait_set])
    has_family_data <- !is.na(fam) &&
      any(traits$species == g$species & traits$family == fam)

    status <- if (!(g$species %in% species_seen)) {
      "not_covered"
    } else if (identical(g$expectation, "present")) {
      if (present) "pass" else if (has_family_data) "fail" else "no_data"
    } else {
      ## An "absent" expectation is satisfied by silence as much as by a
      ## contradicting call, so no_data does not apply.
      if (present) "fail" else "pass"
    }
    data.frame(species = g$species, trait_set = g$trait_set,
               expectation = g$expectation,
               observed = if (present) "present" else "absent",
               status = status, source = g$source, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

## Pairwise Jaccard overlap between sets. Near-identical pairs are the signature
## of an unmerged synonym: this is the check that finds the next
## D-glucose/glucose without anyone having to notice it by eye.
## Implemented as a sparse crossprod rather than a double loop. At 10,348 sets a
## nested loop is ~54M pairwise intersect() calls and does not finish in
## reasonable time; the sparse incidence matrix gives every intersection count in
## one multiplication and only the non-zero cells are ever examined.
set_redundancy <- function(sets, threshold = 0.9, min_size = 3L) {
  sizes <- vapply(sets, length, integer(1))
  keep <- sets[sizes >= min_size]
  nm <- names(keep)
  if (length(nm) < 2L) return(NULL)

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    warning("Matrix not installed; skipping the redundancy report", call. = FALSE)
    return(NULL)
  }

  taxa <- unique(unlist(keep, use.names = FALSE))
  i <- unlist(lapply(keep, function(s) match(s, taxa)), use.names = FALSE)
  j <- rep(seq_along(keep), vapply(keep, length, integer(1)))
  inc <- Matrix::sparseMatrix(i = i, j = j, x = 1,
                              dims = c(length(taxa), length(keep)))

  ## Intersection counts for every set pair that shares at least one taxon.
  inter <- Matrix::crossprod(inc)
  tri <- Matrix::which(Matrix::triu(inter, k = 1L) > 0, arr.ind = TRUE)
  if (!nrow(tri)) return(NULL)

  n <- vapply(keep, length, integer(1))
  a <- tri[, 1L]; b <- tri[, 2L]
  shared <- inter[cbind(a, b)]
  jac <- shared / (n[a] + n[b] - shared)
  sel <- jac >= threshold
  if (!any(sel)) return(NULL)

  out <- data.frame(
    set_a = nm[a[sel]], set_b = nm[b[sel]],
    n_a = n[a[sel]], n_b = n[b[sel]],
    n_shared = as.integer(shared[sel]), jaccard = round(jac[sel], 4),
    stringsAsFactors = FALSE
  )
  out[order(-out$jaccard, -out$n_shared), , drop = FALSE]
}

## Annotation-depth report. The point is to make the bias legible: v7 gave
## Escherichia coli 45 sets and Faecalibacterium prausnitzii 1, which means an
## enrichment result was partly a readout of how well studied a taxon is.
annotation_depth <- function(traits) {
  t2 <- traits[traits$call == "pos", , drop = FALSE]
  n_sets <- tapply(t2$set_name, t2$species, function(z) length(unique(z)))
  n_str  <- tapply(traits$n_tested, traits$species, function(z) max(z, na.rm = TRUE))
  sp <- names(n_sets)
  out <- data.frame(species = sp, n_sets = as.integer(n_sets[sp]),
                    max_strains = as.integer(n_str[sp]), stringsAsFactors = FALSE)
  out[order(-out$n_sets), , drop = FALSE]
}

coverage_report <- function(traits, sets, prov) {
  sizes <- vapply(sets, length, integer(1))
  data.frame(
    metric = c("n_sets", "n_unique_taxa", "n_species_with_calls",
               "n_measured_sets", "n_predicted_sets",
               "n_negative_polarity_sets",
               "median_set_size", "max_set_size", "n_singleton_sets",
               "n_families"),
    value = c(length(sets), length(unique(unlist(sets))),
              length(unique(traits$species)),
              sum(prov$evidence == "measured"), sum(prov$evidence == "predicted"),
              sum(grepl("_negative$", names(sets))),
              stats::median(sizes), max(sizes), sum(sizes == 1L),
              length(unique(traits$family))),
    stringsAsFactors = FALSE
  )
}
