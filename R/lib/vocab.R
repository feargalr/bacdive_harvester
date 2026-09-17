## lib/vocab.R -- normalise BacDive's raw vocabulary before sets are generated.
##
## Several BacDive fields are semi-structured free text rather than controlled
## vocabulary, and treating them as categorical manufactures sets that are not
## distinct biological traits. Observed in the 75k-strain harvest:
##
##   * HTML leaks through:  "quinones: MK-8(H<sub>4</sub>)",
##     "MK-8(H<sub>4</sub>, &omega;-cycl)" -- 1,889 observation values
##   * component order is arbitrary, so the same profile becomes two sets:
##     "quinones: MK-9(H<sub>6</sub>), MK-9(H<sub>8</sub>)"   (56 strains)
##     "quinones: MK-9(H<sub>8</sub>), MK-9(H<sub>6</sub>)"   (30 strains)
##   * a typed prefix carries the real concept: "quinones: MK-7",
##     "Ability for Yeast lysis", "aggregates in chains"
##   * numeric ranges appear in a dozen spellings: 633 distinct halophily
##     concentrations including "0-2 %", "0-2.0 %", "0-2 %(w/v)", "0.0-2.0 %",
##     "02-03 %"
##
## The pipeline order is: decode HTML -> standardise Unicode -> standardise
## punctuation and case -> split structured lists -> sort unordered components ->
## map synonyms -> canonical name.

## ---- HTML and Unicode -------------------------------------------------------

HTML_ENTITIES <- c(
  "&amp;" = "&", "&lt;" = "<", "&gt;" = ">", "&quot;" = "\"", "&apos;" = "'",
  "&nbsp;" = " ", "&ndash;" = "-", "&mdash;" = "-", "&minus;" = "-",
  "&alpha;" = "alpha", "&beta;" = "beta", "&gamma;" = "gamma",
  "&delta;" = "delta", "&epsilon;" = "epsilon", "&omega;" = "omega",
  "&psi;" = "psi", "&iota;" = "iota", "&kappa;" = "kappa", "&lambda;" = "lambda",
  "&mu;" = "mu", "&sup2;" = "2", "&sup3;" = "3", "&deg;" = "degrees",
  "&times;" = "x", "&plusmn;" = "+/-", "&prime;" = "'"
)

## Greek and other non-ASCII characters BacDive uses inline. Mapped to the ASCII
## spellings already used by the enzyme/metabolite alias tables so that e.g.
## "beta-galactosidase" and a Unicode beta collapse to one analyte.
UNICODE_MAP <- c(
  "\u03b1" = "alpha", "\u03b2" = "beta", "\u03b3" = "gamma", "\u03b4" = "delta",
  "\u03b5" = "epsilon", "\u03c9" = "omega", "\u03bc" = "mu", "\u03c8" = "psi",
  "\u00df" = "beta",                                  # sharp s, used for beta
  "\u00b5" = "mu", "\u00b0" = "degrees", "\u00d7" = "x", "\u00b1" = "+/-",
  "\u2013" = "-", "\u2014" = "-", "\u2212" = "-",      # dashes
  "\u2018" = "'", "\u2019" = "'", "\u201c" = "\"", "\u201d" = "\"",
  "\u00a0" = " ", "\u2009" = " ", "\u200b" = "",
  "\u00b2" = "2", "\u00b3" = "3"
)

## Subscript/superscript markup carries meaning in quinone names -- MK-9(H2) is
## not MK-9 -- so the tag is stripped but the digit kept.
decode_html <- function(x) {
  x <- as.character(x)
  x <- gsub("<\\s*/?\\s*(sub|sup|i|b|em|strong|br|p)\\s*/?\\s*>", "", x,
            ignore.case = TRUE, perl = TRUE)
  for (e in names(HTML_ENTITIES)) {
    x <- gsub(e, HTML_ENTITIES[[e]], x, fixed = TRUE)
  }
  ## Numeric entities, then anything else tag-shaped.
  x <- gsub("&#x?[0-9a-fA-F]+;", "", x, perl = TRUE)
  x <- gsub("<[^>]*>", "", x, perl = TRUE)
  x
}

standardise_unicode <- function(x) {
  x <- as.character(x)
  for (u in names(UNICODE_MAP)) {
    x <- gsub(u, UNICODE_MAP[[u]], x, fixed = TRUE)
  }
  x
}

## Punctuation and spacing only -- case is handled dataset-wide by
## collapse_case() so that readable spellings survive where they are unambiguous.
standardise_punct <- function(x) {
  x <- gsub("\\s+", " ", x, perl = TRUE)
  x <- gsub("\\s*([,;])\\s*", "\\1 ", x, perl = TRUE)   # one space after a comma
  x <- gsub("\\s*-\\s*", "-", x, perl = TRUE)           # "alpha - GAL" -> "alpha-GAL"
  x <- gsub("\\(\\s+", "(", x, perl = TRUE)
  x <- gsub("\\s+\\)", ")", x, perl = TRUE)
  x <- gsub("[.,;]+$", "", x, perl = TRUE)              # trailing punctuation
  trimws(x)
}

## The single entry point applied to every raw value before anything else looks
## at it.
clean_vocab <- function(x) {
  standardise_punct(standardise_unicode(decode_html(x)))
}

## ---- structured values ------------------------------------------------------

## Many observation values are "<concept>: <value>" -- "quinones: MK-7". Split so
## the concept becomes the set family and the value becomes the analyte, instead
## of the whole string becoming one opaque set name.
split_typed_prefix <- function(x) {
  x <- clean_vocab(x)
  has <- grepl("^[A-Za-z][A-Za-z0-9 /_-]{2,30}:\\s*.+$", x, perl = TRUE)
  prefix <- ifelse(has, sub("^([^:]+):.*$", "\\1", x), NA_character_)
  value  <- ifelse(has, sub("^[^:]+:\\s*", "", x), x)
  data.frame(prefix = trimws(prefix), value = trimws(value),
             stringsAsFactors = FALSE)
}

## Split a comma/semicolon separated list into components. `mode`:
##   "explode" -- return every component separately (one set per quinone, which
##                is what makes MK-7 comparable across strains)
##   "sort"    -- keep the combination but order it canonically, so
##                "MK-9(H6), MK-9(H8)" and "MK-9(H8), MK-9(H6)" agree
split_components <- function(x, mode = c("explode", "sort")) {
  mode <- match.arg(mode)
  x <- clean_vocab(x)
  parts <- strsplit(x, "\\s*[,;]\\s*", perl = TRUE)
  if (mode == "sort") {
    return(vapply(parts, function(p) {
      p <- trimws(p); p <- p[nzchar(p)]
      if (!length(p)) return(NA_character_)
      paste(sort(unique(p)), collapse = ", ")
    }, character(1)))
  }
  lapply(parts, function(p) {
    p <- trimws(p); unique(p[nzchar(p)])
  })
}

## ---- numeric ranges (halophily, pH, temperature) ----------------------------

## Parse a BacDive concentration/range string to numeric low and high.
## Handles "10 %", "0-4 %", "0-2 %(w/v)", "0.0-2.0 %", "02-03 %", "6.5 %",
## ">10 %", "up to 5 %".
## Molar and g/L values are converted to % (w/v) so every concentration lands on
## one scale: 1 % (w/v) = 10 g/L, and 1 M NaCl = 58.44 g/L = 5.844 %.
NACL_MOLAR_MASS <- 58.44

parse_range <- function(x) {
  s <- clean_vocab(x)
  s <- gsub("\\(w/v\\)|\\(v/v\\)|\\(wt/vol\\)", "", s, ignore.case = TRUE, perl = TRUE)

  ## Detect the unit before stripping it, so a scale factor can be applied.
  unit_molar <- grepl("(^|[^A-Za-z])M($|[^A-Za-z])", s, perl = TRUE)
  unit_gl    <- grepl("g\\s*/\\s*L", s, ignore.case = TRUE, perl = TRUE)
  scale <- ifelse(unit_molar, 100 * NACL_MOLAR_MASS / 1000,
                  ifelse(unit_gl, 1 / 10, 1))

  s <- gsub("g\\s*/\\s*L", "", s, ignore.case = TRUE, perl = TRUE)
  s <- gsub("(^|[^A-Za-z])M($|[^A-Za-z])", "\\1\\2", s, perl = TRUE)
  s <- gsub("%", "", s, fixed = TRUE)
  s <- gsub("up to", "0-", s, ignore.case = TRUE, perl = TRUE)
  s <- gsub("(>=|>|at least)", "", s, perl = TRUE)
  s <- gsub("(<=|<)", "0-", s, perl = TRUE)
  ## "2 .0 %" -- a stray space inside the number.
  s <- gsub("([0-9])\\s+[.]\\s*([0-9])", "\\1.\\2", s, perl = TRUE)
  s <- trimws(s)

  num <- "[0-9]+(?:[.][0-9]+)?"
  lo <- rep(NA_real_, length(s)); hi <- rep(NA_real_, length(s))

  ## a-b  (the hyphen is a range separator, never a minus sign here)
  m2 <- regmatches(s, regexec(sprintf("^(%s)\\s*-\\s*(%s)$", num, num), s, perl = TRUE))
  ## single value
  m1 <- regmatches(s, regexec(sprintf("^(%s)$", num), s, perl = TRUE))

  for (i in seq_along(s)) {
    if (length(m2[[i]]) == 3L) {
      lo[i] <- as.numeric(m2[[i]][2]); hi[i] <- as.numeric(m2[[i]][3])
    } else if (length(m1[[i]]) == 2L) {
      lo[i] <- hi[i] <- as.numeric(m1[[i]][2])
    }
  }
  lo <- lo * scale
  hi <- hi * scale

  ## "02-03" style zero padding parses fine as numbers; guard inverted ranges.
  swap <- !is.na(lo) & !is.na(hi) & lo > hi
  if (any(swap)) { tmp <- lo[swap]; lo[swap] <- hi[swap]; hi[swap] <- tmp }
  data.frame(low = lo, high = hi, stringsAsFactors = FALSE)
}
