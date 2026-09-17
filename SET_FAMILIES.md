# What the BacDive sets are

Reference build: **3,165 sets · 8,587 taxa · 25 trait families**, from a
592-genus / 75,062-strain harvest of ~73% of BacDive's 102,187 strains.

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
| **Comparable but sparsely measured** | `Util_turanose` — a real assay few species were tested for | **Keep.** TaxSEA intersects sets with the observed taxa *before* size-filtering, so sparse sets never reach testing and cost nothing. They may become useful as BacDive grows. |
| **Not comparable at all** | `antibiotic A-195`; `Ability for Yeast lysis` (1,207 strain records but only 6 species); a polar-lipid prose fragment | **Drop, or salvage a curated subset via a semantic allowlist.** No amount of data makes these comparable across taxa. |

**28 families are disabled** on this principle, each with its rationale recorded in
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
| `review_sets_summary.tsv` | set (5,268) | triage: fill in `keep`, `merge_into`, `notes` |
| `review_analytes.tsv` | analyte (3,570) | spot unmerged synonyms: fill in `canonical_name` |
| `review_sets_long.tsv.gz` | (set, taxon), 479k rows | check actual membership |
| `08_sets_reaching_testing.tsv` | the **389** sets that actually get tested on an IBD dataset | **curate this first** — the other ~4,900 cost nothing |

## Reading a set name

```
BacDive_Util_glucose                positive: species utilise glucose
BacDive_Util_glucose_negative       tested and did NOT utilise it
BacDive_Util_nitrate_gas            context suffix: "builds gas from" nitrate
BacDive_Enzyme_catalase             positive enzyme activity
BacDive_Temp_thermophile            derived from a numeric field
BacDive_Shape_coccoid               a derived union of coccus and ovoid
BacDive_Pred_*                      genome-based prediction (none in this build)
```

**2,444 sets are `_negative`.** v7 could not express these: it recorded only `+`
results, so "tested and negative" was indistinguishable from "never tested". The
`n_species_tested` column gives you the denominator.

## The families

### Headline categoricals — small, clean, well populated

| Family | Sets | Largest | Notes |
|---|---|---|---|
| `oxygen` | 4 | 3,095 | aerobe 3,173 · anaerobe 1,260 · facultative anaerobe 973 · microaerophile 639. Resolved as a **capability**, not a majority vote: `facultative anaerobe` is flagged `dominant`, because a facultative organism grown only aerobically is recorded "aerobe". Without that, *Citrobacter freundii*, *Enterobacter cloacae* and *Serratia marcescens* all came out aerobic. Dominance needs ≥10% of strains and ≥2 strains, or one stray annotation flips a species (*P. aeruginosa* became facultative on 1 of 142). |
| `gram` | 3 | 2,797 | positive / negative / variable |
| `motility` | 2 | 3,040 | yes / no. v7 shipped only `yes`, discarding 399 `no` calls. |
| `spore` | 2 | 1,649 | yes / no. Absent from v7 entirely. |
| `shape` | 12 | 3,925 | rod 3,925 · coccus 641 · ovoid 274 · filament · spiral · curved … |
| `shape_union` | 1 | 914 | **`Shape_coccoid` = coccus ∪ ovoid.** The constituents stay separate: `ovoid` carries real information (it is largely *Acinetobacter*, which are coccobacilli, not cocci), so collapsing them would misclassify more than it fixes. The union gives the broader category that is usually the biological question, without destroying the source distinction. |
| `pathogen_human` | 2 | 395 | documented human pathogenicity |

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
| `culture_temp` | 8 | 767, over 189k rows — the largest field in BacDive | mesophile 7,835 · thermophile 1,017 · psychrophile 76 · hyperthermophile 18 · thresholds ≥45/55/65/80 °C |
| `culture_temp_min` | 1 | — | psychrotolerant 975 |
| `halophily` | 10 | **633** (`0-2 %`, `0-2.0 %`, `0-2 %(w/v)`, `02-03 %`, `0-1.1 M`, `5-150 g/L` …) | Kushner classes + thresholds ≥2/5/10/15/20% NaCl. All 633 parse; M converted at 58.44 g/mol, g/L at ÷10. NaCl only — marine salts, MgCl₂ and KCl are different physiology and far too sparse. |
| `culture_ph_max` / `_min` | 4 | 456 | alkaliphile 1,477 · acidophile 875 · thresholds ≥pH 9/10 |
| `gc_content` | 3 | 1,268 | GC_low <40% 1,202 · GC_mid 2,036 · GC_high >60% 3,055 |

### Metabolism — the bulk

| Family | Sets | Largest | Notes |
|---|---|---|---|
| `util` | 4,457 | 2,601 | metabolite utilisation. Median 3 taxa: a long sparse tail that is *comparable but rarely measured*, so it is kept and simply never reaches testing. |
| `enzyme` | 191 | 3,697 | enzyme activities |
| `api_zym` | 38 | 3,877 | the one API panel using real analyte names; merged into the `Enzyme` family |
| `production` | 398 | 3,790 | metabolites **produced**, not consumed — butyrate, propionate, acetate, lactate, indole, H₂S. The family that speaks most directly to metabolomics. |
| `metabolite_test` | 8 | 2,172 | Voges-Proskauer, methyl red, indole, citrate |

Synonym work in `util`/`enzyme`, all of which v7 shipped as split sets:
`D-glucose` + `glucose` + `alpha-D-glucose` → one set (v7: 311/124/1);
`esterase lipase (C 8)` + `esterase Lipase (C 8)` + `Esterase Lipase` → one, merged
across the `enzymes` subsection *and* the `API zym` panel; `oxidase` +
`cytochrome oxidase` + `cytochrome-c oxidase` → one. D/L enantiomers are **kept
distinct** where both are genuinely assayed (arabinose, xylose, arabitol, lactate,
malate, tartrate). A `family` of `api_zym|enzyme` in the summary sheet is the
cross-family merge working.

Assay *context* is preserved, so `Util_nitrate_gas` ("builds gas from nitrate") is
no longer folded into nitrate utilisation as it was in v7's 219-member
`Utilizes_nitrate`.

### Allowlisted families — salvaged subsets

| Family | Sets | What survived |
|---|---|---|
| `observation` | 5 | From 1,502 analytes (85% single-species). Kept: **menaquinone 1,670** and **ubiquinone 367** (101 MK-n and 8 Q-n variants collapsed; only 8 species carry both, so it is a real dichotomy in respiratory chemistry), and cell aggregation — chains 436 · clumps 245 · none 11. Discarded: polar-lipid profiles (172 analytes over 143 species, 64 of them prose fragments), yeast/*E. coli* lysis (1,207 strain records but 6 species), one-off observations. |
| `compound_class` | 14 | From 1,887 raw analytes. steroid transformation 15 · siderophore 12 · antibiotic (unspecified) 11 · polyhydroxyalkanoate 11 · ethanol 9 · pigment 8 · vitamin B12 8 · carotenoid 7 · bacteriocin 6 · biosurfactant 3 · butanediol 3 · exopolysaccharide 3 · propanediol 2 · methanol 1. **Acetoin deliberately excluded** — it duplicates `Test_voges-proskauer-test`. Enzyme-like entries needed no rehoming: β-lactamase, amylase, catalase, urease, protease, nitrate reductase, cellulase, chitinase and lipase were already in `enzyme`. |
| `tolerance` | 3 | lysozyme 30 · tellurite 8 · bile 4. From 27 analytes, 70% single-species heavy metals and selective agents. |

### Kept as-is

`murein` (89 sets, largest 108) — peptidoglycan type is a comparable
chemotaxonomic property genuinely shared across taxa, not a one-off observation.
`nutrition_type` (13 sets) — chemoorganotroph, chemoheterotroph, methylotroph: a
small controlled vocabulary of comparable metabolic strategies.

## Limitations to keep in mind

**BacDive is silent for some key taxa.** 2,890 of 8,935 species (32%) have **no
oxygen tolerance record at all** — including *Faecalibacterium prausnitzii*,
*Roseburia intestinalis* and *Ruminococcus bromii*, across every harvested strain.
This is a hard ceiling, and it is why v7's `BacDive_anaerobe` membership for those
species could only have come from the retired GPT table.

**Everything here is lab-measured.** Zero `confidence` fields and no
`Genome-based predictions` section appear in the 11.3M-row harvest, because
`BacDive::fetch()` does not pass the API's `predictions` parameter. A clean claim
for the paper, and an untapped coverage lever.

**Use the evidence columns.** `median_n_tested` / `max_n_tested` say how many
strains stand behind a call; a median of 1 rests on single strains. Annotation
depth is very uneven — see `output/05_annotation_depth.csv`.

**Testing burden.** On HMP_2019_ibdmdb (1,201 IBD vs 426 controls, 288 observed
taxa): **389 of 5,268 sets reach testing** at `min_set_size = 3`, against 121 of
267 for v7. 65 are significant at FDR < 0.05, a 16.5% hit rate — down from 39%
before the API code panels were dropped, but still elevated because `util` sets are
correlated with one another. `output/05_redundant_sets.csv` lists 126 set pairs
above Jaccard 0.9 if you want to prune further.

**A documented divergence, not a failure.** BacDive records *Veillonella parvula*
as `oval-shaped` → canonicalised to `ovoid`, while the literature describes it as a
Gram-negative anaerobic coccus. The gold standard follows the source and asserts
`Shape_ovoid` plus the `Shape_coccoid` union, rather than overriding the data —
maintaining species-level overrides would let the validation set define the
database. Gold standard: **36 pass, 0 fail**, 3 no-data, 2 not-covered.
