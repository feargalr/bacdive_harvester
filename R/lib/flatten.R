## lib/flatten.R -- turn a BacDive record (nested list) into a tidy long table.
##
## Why generic rather than named accessors: BacDive returns only the fields that
## have content, and any subsection may arrive either as a single named list
## (one entry) or as an unnamed list of named lists (several entries). The v7
## helpers handled that with a per-trait if/else, which is why adding a trait
## meant writing new extraction code and why an unexpected shape silently became
## "not reported". Here we walk the whole record once, losslessly, and let a
## declarative spec (data/curation/trait_spec.csv) decide what to do with it.

## Unit separator (0x1F): cannot occur in a BacDive key, so it is safe as a
## path delimiter.
SEP <- rawToChar(as.raw(31))

## Recursively reduce a nested list to a character matrix of (path, value).
## Leaves are atomic scalars. Atomic vectors of length > 1 are expanded with an
## index or name so nothing is dropped. Empty lists and NULLs vanish, which is
## correct: BacDive omits empty fields anyway.
.walk <- function(x, path = character(0)) {
  if (is.null(x) || length(x) == 0L) return(NULL)

  if (is.atomic(x)) {
    if (length(x) == 1L) {
      return(cbind(paste(path, collapse = SEP), as.character(x)))
    }
    nm <- names(x)
    keys <- if (is.null(nm) || !all(nzchar(nm))) as.character(seq_along(x)) else nm
    return(cbind(paste(c(path, keys), collapse = SEP), as.character(x)))
  }

  if (!is.list(x)) {
    return(cbind(paste(path, collapse = SEP), as.character(x)))
  }

  nm <- names(x)
  unnamed <- is.null(nm) || !any(nzchar(nm))
  out <- vector("list", length(x))
  for (i in seq_along(x)) {
    key <- if (unnamed || !nzchar(nm[[i]])) sprintf("[[%d]]", i) else nm[[i]]
    out[[i]] <- .walk(x[[i]], c(path, key))
  }
  do.call(rbind, out)
}

is_index_token <- function(x) grepl("^\\[\\[[0-9]+\\]\\]$", x)

.empty_long <- function() {
  data.frame(bacdive_id = character(0), section = character(0),
             subsection = character(0), entry = integer(0),
             field = character(0), value = character(0), ref = character(0),
             stringsAsFactors = FALSE)
}

## Flatten one record into a data frame:
##   section     top-level BacDive section, e.g. "Physiology and metabolism"
##   subsection  e.g. "metabolite utilization"; the empty string when the field
##               sits directly on the section (e.g. General$"BacDive-ID").
##               Deliberately "" and not NA: `long$subsection == "x"` on an NA
##               returns NA, which silently injects all-NA rows when used as a
##               subsetting index -- a trap for anyone querying this table.
##   entry       1-based index of the entry within the subsection
##   field       remaining path inside the entry, dot-joined
##   value       the leaf value, as character
##   ref         the entry's @ref, carried alongside every field of that entry
flatten_record <- function(rec, bacdive_id = NA_character_) {
  m <- .walk(rec)
  if (is.null(m) || !nrow(m)) return(.empty_long())

  parts <- strsplit(m[, 1], SEP, fixed = TRUE)
  n <- length(parts)

  section    <- character(n)
  subsection <- character(n)
  entry      <- integer(n)
  field      <- character(n)

  for (i in seq_len(n)) {
    p <- parts[[i]]
    section[i] <- p[[1]]
    r <- p[-1]

    if (!length(r)) {
      ## A scalar sitting directly at the top level of the record.
      subsection[i] <- ""
      entry[i] <- 1L
      field[i] <- section[i]
      next
    }

    if (is_index_token(r[[1]])) {
      ## The section itself was an unnamed list of entries.
      subsection[i] <- ""
      entry[i] <- as.integer(gsub("\\D", "", r[[1]]))
      tail <- r[-1]
      field[i] <- if (length(tail)) paste(tail, collapse = ".") else section[i]
      next
    }

    tail <- r[-1]
    if (!length(tail)) {
      ## A scalar sitting directly on the section, e.g. General$"BacDive-ID".
      subsection[i] <- ""
      entry[i] <- 1L
      field[i] <- r[[1]]
      next
    }

    subsection[i] <- r[[1]]
    if (is_index_token(tail[[1]])) {
      entry[i] <- as.integer(gsub("\\D", "", tail[[1]]))
      tail <- tail[-1]
    } else {
      entry[i] <- 1L
    }
    field[i] <- if (length(tail)) paste(tail, collapse = ".") else subsection[i]
  }

  df <- data.frame(
    bacdive_id = as.character(bacdive_id),
    section    = section,
    subsection = subsection,
    entry      = entry,
    field      = field,
    value      = m[, 2],
    stringsAsFactors = FALSE
  )

  ## Attach each entry's @ref to every field of that entry, so a downstream
  ## trait call can be traced back to a literature reference.
  df$ref <- NA_character_
  is_ref <- df$field == "@ref"
  if (any(is_ref)) {
    key <- paste(df$section, df$subsection, df$entry, sep = SEP)
    ref_map <- df$value[is_ref]
    names(ref_map) <- key[is_ref]
    ref_map <- ref_map[!duplicated(names(ref_map))]
    df$ref <- unname(ref_map[key])
    df <- df[!is_ref, , drop = FALSE]
  }

  df$value <- clean_str(df$value)
  df <- df[!is.na(df$value), , drop = FALSE]
  rownames(df) <- NULL
  df[, c("bacdive_id", "section", "subsection", "entry", "field", "value", "ref")]
}

## ---- per-strain metadata ----------------------------------------------------
## Pulled with named accessors because these few fields anchor everything else
## and we want a loud failure, not a silent NA, if BacDive renames them.

TAX_SECTION <- "Name and taxonomic classification"

## Species name taken from the record itself -- never from the query string.
## This is the guard the v7 pipeline lacked: request() truncates a taxon query to
## its first three word components and nothing checked what came back.
record_species <- function(rec) {
  tax <- rec[[TAX_SECTION]]
  if (is.null(tax)) return(NA_character_)
  sp <- clean_str(tax$species %||% NA_character_)[[1]]
  if (!is.na(sp)) return(sp)
  ## Fall back to genus + species epithet if `species` is absent.
  g <- clean_str(tax$genus %||% NA_character_)[[1]]
  e <- clean_str(tax[["species epithet"]] %||% NA_character_)[[1]]
  if (!is.na(g) && !is.na(e)) return(paste(g, e))
  NA_character_
}

record_subspecies <- function(rec) {
  tax <- rec[[TAX_SECTION]]
  if (is.null(tax)) return(NA_character_)
  clean_str(tax[["subspecies epithet"]] %||% NA_character_)[[1]]
}

record_is_type_strain <- function(rec) {
  tax <- rec[[TAX_SECTION]]
  if (is.null(tax)) return(NA)
  val <- tax[["type strain"]]
  if (is.null(val) || !length(val)) return(NA)
  ## v7 used `!is.null(val) && tolower(val) == "yes"`, which is a hard error on
  ## R >= 4.3 when the field has length > 1 -- and it sat outside the tryCatch,
  ## so the harvest loop would abort mid-run. any() makes it length-safe.
  any(tolower(as.character(unlist(val))) %in% c("yes", "true", "1"), na.rm = TRUE)
}

## Species-level NCBI taxid straight out of the record. Removes the dependency on
## name-based lookup, which is what broke on reclassified genera
## (Lactobacillus -> Lacticaseibacillus, Eubacterium rectale -> Agathobacter).
record_ncbi_taxid <- function(rec) {
  blk <- rec$General[["NCBI tax id"]]
  if (is.null(blk) || !length(blk)) return(NA_character_)
  entries <- if (!is.null(names(blk)) && any(nzchar(names(blk)))) list(blk) else blk
  lvl <- vapply(entries, function(e) {
    v <- clean_str(e[["Matching level"]])[[1]]
    if (is.na(v)) "" else tolower(v)
  }, character(1), USE.NAMES = FALSE)
  ids <- vapply(entries, function(e) {
    clean_str(as.character(e[["NCBI tax id"]]))[[1]]
  }, character(1), USE.NAMES = FALSE)
  hit <- which(lvl == "species" & !is.na(ids))
  if (length(hit)) return(ids[[hit[[1]]]])
  NA_character_
}

record_bacdive_id <- function(rec) {
  clean_str(as.character(rec$General[["BacDive-ID"]] %||% NA_character_))[[1]]
}

strain_meta <- function(rec) {
  data.frame(
    bacdive_id  = record_bacdive_id(rec),
    species     = record_species(rec),
    subspecies  = record_subspecies(rec),
    ncbi_taxid  = record_ncbi_taxid(rec),
    type_strain = record_is_type_strain(rec),
    stringsAsFactors = FALSE
  )
}

## ---- cache -> long table ----------------------------------------------------

flatten_cache <- function(cache_dir = PIPE$cache_dir) {
  files <- list.files(cache_dir, pattern = "[.]rds$", full.names = TRUE)
  files <- files[basename(files) != basename(PIPE$manifest_file)]
  if (!length(files)) stop("no cached genus files in ", cache_dir,
                           " -- run 01_harvest.R first", call. = FALSE)

  long_parts <- vector("list", length(files))
  meta_parts <- vector("list", length(files))

  for (i in seq_along(files)) {
    blob <- readRDS(files[[i]])
    recs <- blob$records
    if (!length(recs)) next
    ids <- names(recs)
    if (is.null(ids)) ids <- vapply(recs, function(r) record_bacdive_id(r) %||% NA_character_,
                                    character(1))
    long_parts[[i]] <- do.call(rbind, lapply(seq_along(recs), function(j) {
      flatten_record(recs[[j]], ids[[j]])
    }))
    meta_parts[[i]] <- do.call(rbind, lapply(recs, strain_meta))
    log_msg(sprintf("flattened %-28s %4d records", basename(files[[i]]), length(recs)))
  }

  long <- do.call(rbind, long_parts)
  meta <- do.call(rbind, meta_parts)
  rownames(long) <- NULL; rownames(meta) <- NULL

  ## Drop strains we cannot attribute to a species: without a species name the
  ## record cannot be aggregated or mapped to a taxid.
  bad <- is.na(meta$species)
  if (any(bad)) {
    log_msg(sprintf("dropping %d strain(s) with no species name in the record", sum(bad)))
    meta <- meta[!bad, , drop = FALSE]
    long <- long[long$bacdive_id %in% meta$bacdive_id, , drop = FALSE]
  }
  meta <- meta[!duplicated(meta$bacdive_id), , drop = FALSE]

  list(long = long, meta = meta)
}
