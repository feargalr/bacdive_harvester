## lib/normalise.R -- analyte and value normalisation, driven by the curation
## tables in data/curation/. Every judgement call lives in a CSV, not in code,
## so the curation is reviewable and diffable rather than buried in a sequence of
## `df[df == "x"] <- "y"` statements.

NO_CONTEXT <- "<none>"

## ---- curation loading -------------------------------------------------------

## Reject a ragged CSV loudly. An unquoted comma in a rationale column silently
## invented a phantom alias row while this pipeline was being written, so this
## check is not hypothetical.
assert_rectangular <- function(path) {
  n <- utils::count.fields(path, sep = ",", quote = "\"", comment.char = "#")
  n <- n[!is.na(n)]
  off <- which(n != n[[1]])
  if (length(off)) {
    stop(sprintf("%s: %d line(s) have the wrong number of fields (line %s). ",
                 basename(path), length(off), paste(off, collapse = ", ")),
         "Quote any field containing a comma.", call. = FALSE)
  }
  invisible(TRUE)
}

load_curation <- function(name, required_cols) {
  path <- file.path(PIPE$curation_dir, name)
  if (!file.exists(path)) stop("missing curation file: ", path, call. = FALSE)
  assert_rectangular(path)
  x <- utils::read.csv(path, stringsAsFactors = FALSE, colClasses = "character",
                       na.strings = c("NA"), check.names = FALSE, comment.char = "#")
  miss <- setdiff(required_cols, names(x))
  if (length(miss)) {
    stop(sprintf("%s is missing column(s): %s", name, paste(miss, collapse = ", ")),
         call. = FALSE)
  }
  x
}

load_trait_spec <- function() {
  spec <- load_curation("trait_spec.csv",
    c("family", "kind", "section", "subsection", "analyte_field", "result_field",
      "context_field", "context_mode", "positive_pattern", "negative_pattern",
      "set_infix", "emit_negative", "alias_table", "exclude_analyte_pattern",
      "enabled"))
  ## Optional column: when yes, an analyte is kept ONLY if the alias table lists
  ## it. Turns an alias table into a semantic allowlist.
  if (is.null(spec$restrict_to_alias)) spec$restrict_to_alias <- "no"
  ## Optional: restrict a family to, or exclude from it, particular raw assay
  ## contexts. Lets one BacDive subsection feed two families -- `util` takes
  ## everything except hydrolysis, `hydrolysis` takes only that.
  if (is.null(spec$context_include)) spec$context_include <- ""
  if (is.null(spec$context_exclude)) spec$context_exclude <- ""
  ## A row whose fields have shifted still counts the right number of commas if
  ## the shift is compensated elsewhere, so check the flag's value, not just the
  ## shape: an unquoted comma once turned `enabled` into "yes DROPPED: ...".
  bad_flag <- !tolower(spec$enabled) %in% c("yes", "no", "true", "false", "1", "0")
  if (any(bad_flag)) stop("trait_spec.csv: `enabled` must be yes/no; check quoting in: ",
                          paste(spec$family[bad_flag], collapse = ", "), call. = FALSE)
  spec <- spec[tolower(spec$enabled) %in% c("yes", "true", "1"), , drop = FALSE]
  ok_kinds <- c("categorical", "assay", "panel", "label", "typed_label", "numeric")
  badk <- setdiff(unique(spec$kind), ok_kinds)
  if (length(badk)) stop("trait_spec.csv: unknown kind(s): ", paste(badk, collapse = ", "),
                         call. = FALSE)
  if (anyDuplicated(spec$family)) {
    stop("trait_spec.csv: duplicated family id(s): ",
         paste(unique(spec$family[duplicated(spec$family)]), collapse = ", "), call. = FALSE)
  }
  ## `assay` needs both an analyte and a result field; `categorical`/`label` need
  ## a result field; `panel` needs neither (the field name *is* the analyte).
  need_analyte <- spec$kind == "assay" & (is.na(spec$analyte_field) | !nzchar(spec$analyte_field))
  if (any(need_analyte)) stop("trait_spec.csv: assay rows need analyte_field: ",
                              paste(spec$family[need_analyte], collapse = ", "), call. = FALSE)
  need_result <- spec$kind %in% c("assay", "categorical", "label", "typed_label") &
    (is.na(spec$result_field) | !nzchar(spec$result_field))
  if (any(need_result)) stop("trait_spec.csv: rows need result_field: ",
                             paste(spec$family[need_result], collapse = ", "), call. = FALSE)
  rownames(spec) <- NULL
  spec
}

## An alias table maps a raw BacDive spelling to a canonical one. A blank
## canonical means "drop the token" (used for the default utilisation context).
load_alias_map <- function(name) {
  if (is.na(name) || !nzchar(name)) return(NULL)
  x <- load_curation(name, c("raw", "canonical"))
  if (anyDuplicated(tolower(x$raw))) {
    dup <- unique(tolower(x$raw)[duplicated(tolower(x$raw))])
    stop(sprintf("%s: duplicated raw value(s): %s", name, paste(dup, collapse = ", ")),
         call. = FALSE)
  }
  stats::setNames(ifelse(is.na(x$canonical), "", x$canonical), tolower(x$raw))
}

## A value table additionally carries a priority used to break ties when strains
## of one species disagree. Lower number wins.
## A value may additionally be marked `dominant`: if any strain reports it, the
## species takes it regardless of how many strains report something else. Used for
## capability traits, where observing one mode does not exclude the other.
load_value_map <- function(name) {
  if (is.na(name) || !nzchar(name)) return(NULL)
  x <- load_curation(name, c("raw", "canonical", "priority"))
  if (is.null(x$dominant)) x$dominant <- "no"
  x$priority <- suppressWarnings(as.integer(x$priority))
  if (anyNA(x$priority)) stop(name, ": non-integer priority", call. = FALSE)
  x$raw_key <- tolower(x$raw)
  if (anyDuplicated(x$raw_key)) {
    dup <- unique(x$raw_key[duplicated(x$raw_key)])
    stop(sprintf("%s: duplicated raw value(s): %s", name, paste(dup, collapse = ", ")),
         call. = FALSE)
  }
  x$dominant <- tolower(trimws(x$dominant)) %in% c("yes", "true", "1")
  x[, c("raw_key", "canonical", "priority", "dominant")]
}

## A semantic allowlist: an analyte is kept only if listed, and is renamed to the
## curated class. Rows may match exactly (default) or by regular expression, which
## lets a whole family of variants collapse to one biological concept -- every
## MK-n menaquinone to "menaquinone", every Q-n to "ubiquinone".
##
## The guiding rule (from the project owner): a TaxSEA set should describe a
## comparable biological property shared across taxa, not simply something
## somebody once reported about a strain. A frequency threshold cannot make that
## distinction; it would keep an obscure compound that happened to be studied
## three times and discard a real trait measured twice.
load_allowlist <- function(name) {
  if (is.na(name) || !nzchar(name)) return(NULL)
  x <- load_curation(name, c("raw", "canonical"))
  if (is.null(x$match)) x$match <- "exact"
  x$match <- tolower(trimws(ifelse(is.na(x$match), "exact", x$match)))
  bad <- setdiff(x$match, c("exact", "regex"))
  if (length(bad)) stop(name, ": unknown match mode(s): ", paste(bad, collapse = ", "),
                        call. = FALSE)
  x$canonical <- ifelse(is.na(x$canonical), "", x$canonical)
  x[, c("raw", "canonical", "match")]
}

## Returns the canonical class for each input, or NA to discard it.
apply_allowlist <- function(x, allow) {
  out <- rep(NA_character_, length(x))
  if (is.null(allow) || !nrow(allow)) return(out)

  ex <- allow[allow$match == "exact", , drop = FALSE]
  if (nrow(ex)) {
    hit <- match(tolower(x), tolower(ex$raw))
    out[!is.na(hit)] <- ex$canonical[hit[!is.na(hit)]]
  }
  ## Regex rows apply only where no exact row matched, in file order, so a
  ## specific exact rule always beats a broad pattern.
  rx <- allow[allow$match == "regex", , drop = FALSE]
  for (k in seq_len(nrow(rx))) {
    todo <- is.na(out) & !is.na(x)
    if (!any(todo)) break
    m <- todo & grepl(rx$raw[k], x, perl = TRUE, ignore.case = TRUE)
    out[m] <- rx$canonical[k]
  }
  out[!is.na(out) & !nzchar(out)] <- NA_character_
  out
}

## ---- analyte normalisation --------------------------------------------------

## Applied before the alias lookup so that spelling noise does not defeat it.
## Deliberately conservative: it does NOT strip D-/L- prefixes, because
## D- and L-arabinose (and xylose, arabitol, lactate, malate, tartrate) are
## genuinely different assays. Stereochemistry that should be collapsed is
## listed explicitly in metabolite_aliases.csv instead.
canon_whitespace <- function(x) {
  x <- gsub(" ", " ", x, fixed = TRUE)   # non-breaking space
  x <- gsub("\\s*-\\s*", "-", x)              # "alpha- Galactosidase" -> "alpha-Galactosidase"
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

## Returns the normalised spelling plus a flag saying whether it came from a
## curated alias row. The flag matters for display-label choice: a curated
## canonical form should win over an arbitrary panel spelling.
normalise_analyte <- function(x, alias_map = NULL) {
  ## clean_vocab strips HTML, folds Unicode and normalises punctuation before any
  ## alias lookup, so "alpha- Galactosidase" and "MK-8(H<sub>4</sub>)" are already
  ## in canonical form when the alias table is consulted.
  x <- canon_whitespace(clean_vocab(clean_str(x)))
  curated <- rep(FALSE, length(x))
  if (!is.null(alias_map)) {
    hit <- match(tolower(x), names(alias_map))
    repl <- alias_map[hit]
    use <- !is.na(hit) & nzchar(repl)
    x[use] <- unname(repl[use])
    curated <- use
  }
  data.frame(value = x, curated = curated, stringsAsFactors = FALSE)
}

## Canonicalise a categorical value and return its priority alongside, so the
## aggregator can break a tie deterministically.
normalise_value <- function(x, value_map = NULL) {
  x <- canon_whitespace(clean_vocab(clean_str(x)))
  if (is.null(value_map)) {
    return(data.frame(canonical = tolower(x), priority = rep(99L, length(x)),
                      dominant = rep(FALSE, length(x)), stringsAsFactors = FALSE))
  }
  hit <- match(tolower(x), value_map$raw_key)
  data.frame(
    canonical = ifelse(is.na(hit), tolower(x), value_map$canonical[hit]),
    priority  = ifelse(is.na(hit), 99L, value_map$priority[hit]),
    dominant  = ifelse(is.na(hit), FALSE, value_map$dominant[hit]),
    stringsAsFactors = FALSE
  )
}

normalise_context <- function(x, alias_map = NULL) {
  x <- canon_whitespace(clean_str(x))
  x[is.na(x)] <- NO_CONTEXT
  if (is.null(alias_map)) return(x)
  hit <- match(tolower(x), names(alias_map))
  out <- ifelse(is.na(hit), x, unname(alias_map[hit]))
  out[is.na(out)] <- ""
  out
}

## Collapse spellings that differ only by case, dataset-wide.
##
## normalise_analyte() leaves an unmatched spelling untouched, which preserves
## readability (DNase, H2S, N-acetyl-beta-glucosaminidase) but would otherwise
## let "Esterase Lipase" and "esterase lipase" become two sets -- exactly the
## split v7 shipped. So the merge key is the casefolded string, and the display
## label is chosen per group: a curated alias canonical if one is present,
## otherwise the most frequent spelling, ties broken alphabetically so the result
## is deterministic across runs.
##
## This must be applied ACROSS families that share a set_infix, not within one
## family: `enzymes` and the `API zym` panel both produce BacDive_Enzyme_* sets,
## and collapsing them separately leaves the split in place.
collapse_case <- function(analyte, curated = NULL) {
  if (is.null(curated)) curated <- rep(FALSE, length(analyte))
  key <- tolower(analyte)
  ok <- !is.na(key)
  uk <- unique(key[ok])
  labels <- stats::setNames(rep(NA_character_, length(uk)), uk)
  for (k in uk) {
    in_grp <- ok & key == k
    cand <- analyte[in_grp]
    cur  <- curated[in_grp]
    pool <- if (any(cur)) cand[cur] else cand
    tab <- sort(table(pool), decreasing = TRUE)
    top <- names(tab)[tab == tab[[1]]]
    labels[[k]] <- sort(top)[[1]]
  }
  out <- analyte
  out[ok] <- unname(labels[key[ok]])
  out
}

## ---- EC / ChEBI preference --------------------------------------------------
## Where BacDive supplies an identifier we prefer it as the merge key, because it
## is authoritative and immune to spelling. We keep a human-readable label for
## the set name and only use the identifier to detect that two spellings are the
## same thing.
build_id_bridge <- function(analyte, id) {
  keep <- !is.na(id) & nzchar(id) & !is.na(analyte)
  if (!any(keep)) return(NULL)
  tab <- table(id[keep], analyte[keep])
  ## For each identifier pick its most frequent spelling as the canonical label.
  best <- apply(tab, 1L, function(r) colnames(tab)[which.max(r)])
  stats::setNames(as.character(best), rownames(tab))
}

apply_id_bridge <- function(analyte, id, bridge) {
  if (is.null(bridge)) return(analyte)
  hit <- match(id, names(bridge))
  ifelse(!is.na(hit), unname(bridge[hit]), analyte)
}

## ---- set naming -------------------------------------------------------------
## Names must keep the BacDive prefix: TaxSEA finds and replaces this family with
## grepl("BacDive", names(TaxSEA_db)). Predicted traits get a Pred infix so a
## user can exclude them.
##
## Whitespace becomes "_", matching the other TaxSEA_db sources (e.g.
## GutMGene_producers_of_...). The analyte label itself keeps its spaces in the
## trait table and review exports; only the set identifier is normalised.
make_set_name <- function(infix, label, predicted = FALSE, polarity = NULL) {
  parts <- c(PIPE$set_prefix,
             if (predicted) PIPE$predicted_infix else NULL,
             if (!is.na(infix) && nzchar(infix)) infix else NULL,
             label,
             if (!is.null(polarity) && nzchar(polarity)) polarity else NULL)
  x <- paste(parts[nzchar(parts) & !is.na(parts)], collapse = "_")
  ## BacDive writes chemical locants with a space ("1, 2-propandiol"); close
  ## them up, and make any other comma-space a separator ("yes, in single cases").
  x <- gsub("([0-9]),[[:space:]]+([0-9])", "\\1,\\2", x)
  x <- gsub(",[[:space:]]+", "_", x)
  gsub("[[:space:]]+", "_", x)
}
