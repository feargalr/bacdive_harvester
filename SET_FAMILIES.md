# What the BacDive sets are

Reference build: **3,199 sets · 9,045 taxa · 24 enabled trait families**, from a
705-genus / 76,452-strain harvest (about three quarters of BacDive's strains).
**1,776** of those sets have at least 3 members and are what stage 06 ships to
TaxSEA.

Your own numbers will differ with the genus list and the BacDive release. The
counts quoted throughout this document come from that reference build and are
there to give a sense of scale and to justify the curation decisions, not as
guarantees.

## The curation principle

> **A TaxSEA set should describe a comparable biological property shared across
> taxa, not simply something somebody once reported about a strain.**

Every keep/drop decision below follows from this, and it is deliberately *not* a
numeric threshold. A frequency cut-off would retain an obscure compound that
happened to be studied three times and discard a real trait measured twice.

Two cases follow from it, treated differently:

| | Example | Action |
|---|---|---|
| **Comparable but sparsely measured** | `Uses_turanose` — a real assay few species were tested for | **Keep.** TaxSEA intersects sets with the observed taxa *before* size-filtering, so sparse sets never reach testing and cost nothing. They may become useful as BacDive grows. |
| **Not comparable at all** | `antibiotic A-195`; `Ability for Yeast lysis` (1,207 strain records but only 6 species); a polar-lipid prose fragment | **Drop, or salvage a curated subset via a semantic allowlist.** No amount of data makes these comparable across taxa. |

**29 families are disabled** on this principle, each with its rationale recorded in
`data/curation/trait_spec.csv`:

- laboratory-risk classification — `biosafety` (BSL_1 alone had 5,476 members)
- culture-dependent observables — colony colour / shape / size, incubation period,
  cell length / width, growth medium
- where a strain was isolated — `isolation_cat1`, `isolation_cat2`
- better derived from the metagenome itself — `resistance`, `sensitivity`
- panel codes with no readable meaning — 16 of the 17 `API *` panels
- idiosyncratic strain catalogues — raw `compound_prod` (1,887 analytes, 80%
  single-species)

## Review files

| File | One row per | Use it to |
|---|---|---|
| `review_sets_summary.tsv` | set (3,199) | triage: fill in `keep`, `merge_into`, `notes` |
| `review_analytes.tsv` | analyte (2,233) | spot unmerged synonyms: fill in `canonical_name` |
| `review_sets_long.tsv.gz` | (set, taxon), 449,769 rows | check actual membership |
| `05_redundant_sets.csv` | set pair with Jaccard ≥ 0.9 (41) | candidate unmerged synonyms |

## Reading a set name

Names are `BacDive_<Family>_<value>[_<context>][_negative]` and contain no spaces:
whitespace becomes `_`, and the comma-space BacDive writes in chemical locants is
closed up (`1, 2-propandiol` → `1,2-propandiol`).

```
BacDive_Oxygen_facultative_anaerobe  a categorical call
BacDive_Uses_glucose                 grows on / utilises glucose
BacDive_Uses_glucose_negative        tested and did NOT utilise it
BacDive_Uses_glucose_fermentation    context suffix: ferments glucose
BacDive_Uses_nitrate_reduction       context suffix: reduces nitrate
BacDive_Enzyme_catalase              positive enzyme activity
BacDive_Enzyme_esculin_hydrolysis    hydrolysis assays live with enzymes
BacDive_Temp_thermophile             derived from a numeric field
BacDive_Shape_coccoid                a derived union of coccus and ovoid
BacDive_Pred_*                       genome-based prediction (none in this build)
```

**1,372 sets are `_negative`.** The original scraper could not express these: it
recorded only `+` results, so "tested and negative" was indistinguishable from
"never tested". The `n_species_tested` column gives you the denominator.

## The families

### Headline categoricals — small, clean, well populated

| Family | Sets | Largest | Notes |
|---|---|---|---|
| `oxygen` | 4 | 3,213 | aerobe 3,213 · anaerobe 1,317 · facultative anaerobe 962 · microaerophile 628. Resolved as a **capability**, not a majority vote: `facultative anaerobe` is flagged `dominant`, because a facultative organism grown only aerobically is recorded "aerobe". Without that, *Citrobacter freundii*, *Enterobacter cloacae* and *Serratia marcescens* all came out aerobic. Dominance needs ≥10% of strains and ≥2 strains, or one stray annotation flips a species (*P. aeruginosa* became facultative on 1 of 142). |
| `gram` | 3 | 2,936 | positive / negative / variable |
| `motility` | 2 | 3,204 | yes / no. The original scraper shipped only `yes`. |
| `spore` | 2 | 1,755 | yes / no |
| `shape` | 12 | 4,168 | rod 4,168 · coccus 658 · ovoid 287 · filament 29 · curved 24 · spiral 22 … |
| `shape_union` | 1 | 944 | **`Shape_coccoid` = coccus ∪ ovoid.** The constituents stay separate: `ovoid` carries real information (it is largely *Acinetobacter*, which are coccobacilli, not cocci), so collapsing them would misclassify more than it fixes. The union gives the broader category that is usually the biological question, without destroying the source distinction. |
| `pathogen_human` | 2 | 408 | documented human pathogenicity: `yes` 183 · `yes_in_single_cases` 408 |

### Derived numeric traits — semi-structured text turned into biology

BacDive stores these as one entry per tested value, with a `type` saying whether it
is a growth point, optimum, minimum or maximum. Treated categorically they
manufacture hundreds of meaningless sets. Each is now parsed to a number, reduced
to one statistic per species (highest/lowest supporting growth, or median for a
measured quantity), and emitted as conventional thresholds and classes. Configured
by `numeric_traits.csv` + `numeric_classes.csv`. `optimum` is never read as a
maximum.

| Family | Sets | Raw distinct values | Result |
|---|---|---|---|
| `culture_temp` | 8 | 767, over 189k rows — the largest field in BacDive | mesophile 7,829 · thermophile 1,017 · psychrophile 75 · hyperthermophile 18 · thresholds ≥45/55/65/80 °C |
| `culture_temp_min` | 1 | — | psychrotolerant 1,011 |
| `halophily` | 10 | **633** (`0-2 %`, `0-2.0 %`, `0-2 %(w/v)`, `02-03 %`, `0-1.1 M`, `5-150 g/L` …) | Kushner classes (non-halophile 531 · halotolerant 1,223 · slight 1,229 · moderate 580 · extreme 97) + thresholds ≥2/5/10/15/20% NaCl. All 633 parse; M converted at 58.44 g/mol, g/L at ÷10. NaCl only — marine salts, MgCl₂ and KCl are different physiology and far too sparse. |
| `culture_ph_max` / `_min` | 4 | 456 | alkaliphile 1,693 · acidophile 955 · thresholds ≥pH 9/10 |
| `gc_content` | 3 | 1,268 | GC_low <40% 1,245 · GC_mid 2,124 · GC_high >60% 3,087 |

### Metabolism — the bulk

| Family | Prefix | Sets | Shipped (≥3) | Largest | Notes |
|---|---|---|---|---|---|
| `util` | `Uses_` | 2,185 | 1,348 | 2,938 | substrate use, by assay context (below). Median 5 taxa: a long sparse tail that is *comparable but rarely measured*, kept and simply never reaching testing. |
| `enzyme`, `api_zym`, `hydrolysis` | `Enzyme_` | 437 | 240 | 3,828 | enzyme activities. `API zym` is the one API panel using real analyte names and merges into this family; hydrolysis assays (esculin 2,086, hippurate 612 …) are enzyme activities and live here, not under `Uses_`. |
| `production` | `Prod_` | 392 | 53 | 3,980 | metabolites **produced**, not consumed — indole, H₂S, acetoin, nitrite |
| `metabolite_test` | `Test_` | 8 | 8 | 2,274 | Voges-Proskauer, methyl red, indole, citrate |

**Substrate context is collapsed, not discarded.** BacDive records the same
substrate under many assay descriptions ("builds acid from", "energy source",
"carbon source", …). These are collapsed to a small set of biologically distinct
contexts, and the raw description is preserved underneath in `context_raw`:

| Suffix | Meaning | Largest sets |
|---|---|---|
| *(none)* | growth on / utilisation of the substrate | glucose 2,938 · mannose 1,706 · maltose 1,680 |
| `_fermentation` | ferments it — kept separate from acid production | glucose 889 · sucrose 573 · mannitol 383 |
| `_acid` | builds acid from it | glucose 1,568 · maltose 1,323 · fructose 1,281 |
| `_gas` | builds gas from it | thiosulfate 54 · nitrite 20 · nitrate 11 |
| `_reduction` | reduces it | **nitrate 2,198** · nitrite 268 |
| `_respiration` | used in respiration | nitrate 517 |
| `_Nsource` | used as a nitrogen source — kept separate | L-asparagine 56 · L-arginine 54 |

Nitrate is the clearest case of why context matters: reducing nitrate (2,198
species), respiring it (517) and building gas from it (11) are different
physiology, and the original scraper merged them into one 219-member
`Utilizes_nitrate`.

Synonym work, all of which the original scraper shipped as split sets:
`D-glucose` + `glucose` + `alpha-D-glucose` → one set;
`esterase lipase (C 8)` + `esterase Lipase (C 8)` + `Esterase Lipase` → one
(`Enzyme_esterase_lipase`), merged across the `enzymes` subsection *and* the
`API zym` panel; `oxidase` + `cytochrome oxidase` + `cytochrome-c oxidase` → one.
D/L enantiomers are **kept distinct** where both are genuinely assayed (arabinose,
xylose, arabitol, lactate, malate, tartrate). A `family` of `api_zym|enzyme` in the
summary sheet is the cross-family merge working.

### Allowlisted families — salvaged subsets

| Family | Prefix | Sets | What survived |
|---|---|---|---|
| `observation` | `Observation_` | 5 | From 1,502 analytes (85% single-species). Kept: **menaquinone 1,698** and **ubiquinone 382** (101 MK-n and 8 Q-n variants collapsed; very few species carry both, so it is a real dichotomy in respiratory chemistry), and cell aggregation — chains 467 · clumps 255 · none 24. Discarded: polar-lipid profiles (172 analytes over 143 species, 64 of them prose fragments), yeast/*E. coli* lysis (1,207 strain records but 6 species), one-off observations. |
| `compound_class` | `Produces_` | 14 | From 1,887 raw analytes. steroid transformation 15 · siderophore 12 · antibiotic (unspecified) 11 · ethanol 11 · polyhydroxyalkanoate 11 · vitamin B12 10 · carotenoid 8 · pigment 8 · bacteriocin 6 · exopolysaccharide 4 · biosurfactant 3 · butanediol 3 · propanediol 3 · methanol 1. **Acetoin deliberately excluded** here — it is already covered by `Prod_acetoin` and `Test_voges-proskauer-test`. Enzyme-like entries needed no rehoming: β-lactamase, amylase, catalase, urease, protease, nitrate reductase, cellulase, chitinase and lipase were already in `Enzyme_`. |
| `tolerance` | `Tolerance_` | 3 | lysozyme 30 · tellurite 26 · bile 5. From 27 analytes, 70% single-species heavy metals and selective agents. |

### Kept as-is

`murein` (89 sets, largest 128) — peptidoglycan type is a comparable
chemotaxonomic property genuinely shared across taxa, not a one-off observation.
`nutrition_type` (14 sets) — chemoorganotroph, chemoheterotroph, methylotroph: a
small controlled vocabulary of comparable metabolic strategies.

## Limitations to keep in mind

**BacDive is silent for some key taxa.** 3,232 of 9,642 species with any call
(34%) have **no oxygen tolerance record at all** — including *Faecalibacterium
prausnitzii* and *Roseburia intestinalis*, across every harvested strain —
and *Ruminococcus bromii* is absent from the harvest altogether. This is a hard
ceiling, and it is why the original `BacDive_anaerobe` membership for those
species could only have come from the retired LLM-derived table.

**Reclassified taxa need the second pass.** *Agathobacter rectalis*,
*Mediterraneibacter gnavus*, *Enterocloster bolteae* and *Phocaeicola vulgatus* are
reached only because `01b_resolve_synonyms.R` finds their current genus and that
genus is then queried. The gold standard asserts two of them so a regression is
caught.

**Everything here is lab-measured.** Zero `confidence` fields and no
`Genome-based predictions` section appear in the 11.6M-row harvest, because
`BacDive::fetch()` does not pass the API's `predictions` parameter. An untapped
coverage lever.

**Use the evidence columns.** `median_n_tested` / `max_n_tested` say how many
strains stand behind a call; a median of 1 rests on single strains. Annotation
depth is very uneven — see `output/05_annotation_depth.csv`.

**Correlated sets.** Substrate sets are correlated with one another (a species
that ferments glucose often ferments mannose), so enrichment hits within `Uses_`
are not independent. `output/05_redundant_sets.csv` lists the 41 set pairs above
Jaccard 0.9 if you want to prune further.

**A documented divergence, not a failure.** BacDive records *Veillonella parvula*
as `oval-shaped` → canonicalised to `ovoid`, while the literature describes it as a
Gram-negative anaerobic coccus. The gold standard follows the source and asserts
`Shape_ovoid` plus the `Shape_coccoid` union, rather than overriding the data —
maintaining species-level overrides would let the validation set define the
database. Gold standard: **38 pass, 0 fail**, 3 no-data, 2 not-covered.
