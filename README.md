# bacdive_harvester

Build reproducible taxon sets from [BacDive](https://bacdive.dsmz.de), the
Bacterial Diversity Metadatabase, for use with
[TaxSEA](https://github.com/feargalr/TaxSEA) or any other set-enrichment tool.

Harvest once, then rebuild the sets offline as often as you like. Only the
harvest touches the network; everything downstream reads a local cache, so the
output is deterministic and the curation is reviewable.

## Install

Only one package is required, and only for the harvest step:

```r
## the DSMZ API client; it is on R-Forge, not CRAN
install.packages("BacDive", repos = "https://R-Forge.R-project.org")
```

`Matrix` is also used, but it ships with R.

Optional, each unlocking one extra step:

| Package | Needed for | Without it |
|---|---|---|
| `curatedMetagenomicData` | building the target list from public metagenomes | supply your own species list |
| `TaxSEA` | seeding the target list from `TaxSEA_db` | that source skipped |

Everything else — flattening, aggregation, set building, validation, export — is
base R.

## Quick start

Register for BacDive API access at <https://api.bacdive.dsmz.de>, then:

```bash
cp .Renviron.example .Renviron     # add your credentials, then restart R
Rscript R/01a_target_list.R        # which genera to harvest
Rscript R/01_harvest.R             # the only network step; resumable
Rscript run_all.R                  # stages 02-07
```

Credentials are read from the environment and never appear in the code.
`.Renviron` is gitignored.

### No Bioconductor?

Write a one-column CSV to `data/target_species.csv` with header `species` and
values like `Bacteroides uniformis`, then:

```bash
TAXSEA_USE_USER_LIST=1 Rscript R/01a_target_list.R
```

The genus list is derived from it and nothing else is needed.

## Pipeline

```
01a_target_list.R     species list        -> data/target_genera.csv
01b_resolve_synonyms.R reclassified names -> data/target_genera_extra.csv   [network: NCBI]
01_harvest.R          BacDive API         -> cache/<Genus>.rds              [network: BacDive]
02_flatten.R          records             -> output/02_strain_long.rds
03_aggregate.R        strains             -> output/03_species_traits.rds
04_build_sets.R       trait calls         -> output/04_taxon_sets.rds
05_validate.R         assertions, reports -> output/05_*.csv
06_merge_taxsea_db.R  sets                -> output/TaxSEA_db.rda, NCBI_ids.rda  [needs TaxSEA data]
07_export_review.R    sets                -> output/review_*.tsv
```

`01b_resolve_synonyms.R` is a second pass, not part of the main sequence. Run it
after a first harvest to find target species missed because they were
reclassified into a genus you did not query (in the reference build it found
*Agathobacter*, *Enterocloster*, *Mediterraneibacter*, *Phocaeicola* and the
*Lactobacillus* split). Merge its output into `data/target_genera.csv` and re-run
the harvest, which skips genera already cached. Repeat until it reports no new
genera; the reference build needed two rounds.

See `SET_FAMILIES.md` for what the sets are and why each family was kept, collapsed
or dropped.

Each stage reads the previous stage's cached artefact, so any subset can be
re-run. `05_validate.R` exits non-zero if a gold-standard assertion fails, which
makes it usable as a release gate.

## Building TaxSEA's data files

Stage 06 turns the sets into drop-in replacements for TaxSEA's `TaxSEA_db.rda` and
`NCBI_ids.rda`. It never touches a package in place: it writes to `output/`, and
copying the two files into TaxSEA's `data/` is a deliberate manual step.

```bash
TAXSEA_DATA_DIR=~/TaxSEA/data Rscript R/06_merge_taxsea_db.R
```

| Variable | Default | Meaning |
|---|---|---|
| `TAXSEA_DATA_DIR` | installed TaxSEA | directory holding the current `TaxSEA_db.rda` and `NCBI_ids.rda` |
| `TAXSEA_SHIP_MIN_MEMBERS` | `3` | smallest set shipped. TaxSEA intersects sets with the observed taxa before size-filtering, so a smaller set can never be tested |
| `TAXSEA_HARVESTER_COMMIT` | `git rev-parse` | recorded in the build info |

What it does:

- replaces every `BacDive_*` set and keeps every other source unchanged, refusing
  to run if a new name would collide with an existing one;
- extends `NCBI_ids` with the species names behind the shipped sets, plus former
  names resolved by `01b`, in both `Genus species` and `Genus_species` forms.
  Binomials only. **An existing taxid is never overwritten**: a name that maps
  today maps the same way afterwards, and disagreements are written to
  `NCBI_ids_conflicts.tsv` instead. A name already present with no taxid is a dead
  end in TaxSEA, so that is filled. In the reference build this took the share of
  shipped taxa reachable by name from 15% to 96%;
- writes `BacDive_set_provenance.tsv` (per set: family, evidence, strain and
  species counts) and `BacDive_build_info.tsv` (harvester commit, BacDive client
  version, harvest dates, counts) for the package to carry.

## Design decisions

**Harvest by genus, not species.** BacDive's taxon endpoint accepts a genus and
`retrieve()` paginates the whole result, so one call returns every strain of
every species in that genus. Fewer calls, more strains, and species you did not
think to ask for arrive for free.

**Use every strain.** A species-level call is a vote across all of its strains,
carrying `n_pos` / `n_neg` / `n_tested`. Reading one "representative" strain
makes the answer depend on API result ordering and cannot express a tested
negative.

**Take identity from the record.** Species name comes from
`Name and taxonomic classification$species` and the NCBI taxid from
`General$"NCBI tax id"` where `Matching level == "species"` — never from the
query string. BacDive's taxon endpoint truncates a query to its first three word
components and does not promise the result matches what you asked for.

**Curation lives in CSVs.** Every judgement call — synonym maps, oxygen priority,
numeric thresholds, semantic allowlists — is a row in `data/curation/`, not a
line of code. They are reviewable and diffable.

**One governing principle**, recorded at the top of `data/curation/trait_spec.csv`:

> A TaxSEA set should describe a comparable biological property shared across
> taxa, not simply something somebody once reported about a strain.

## Adding a trait

Add a row to `data/curation/trait_spec.csv`. No code change.

| kind | meaning | example |
|---|---|---|
| `categorical` | one value per strain, resolved by majority then curated priority | oxygen tolerance, Gram stain |
| `assay` | analyte + `+`/`-` result in the same entry | metabolite utilisation, enzymes |
| `panel` | flat block where each field name *is* the analyte | `API zym` |
| `label` | presence-only descriptor | nutrition type |
| `typed_label` | `<concept>: <value>` lists, split and normalised | quinones, murein |
| `numeric` | one entry per tested value, reduced to thresholds and classes | temperature, pH, NaCl, GC content |

## Data cleaning worth knowing about

BacDive is a curated database, but several fields are semi-structured free text.
The pipeline handles these; the details are documented in
`data/curation/*.csv`:

- HTML markup in values (`MK-8(H<sub>4</sub>)`, `&omega;-cycl`)
- component order varying within a value, so the same profile becomes two sets
- numeric ranges in many spellings and units (`0-2 %`, `0-2.0 %(w/v)`, `0-1.1 M`,
  `5-150 g/L`) — all parsed to one scale
- the same assay recorded under several substrate spellings
- assay context distinguishing genuinely different measurements: growth on a
  substrate, acid production from it, gas production from it and hydrolysis of it
  are not the same trait

## Tests

```bash
Rscript tests/run_tests.R
```

79 assertions, base R only, no network. Runs stages 02–06 against a synthetic
cache built from real BacDive record shapes
(`https://api.bacdive.dsmz.de/v2/example/fetch/24493`). Each assertion
corresponds to a specific defect, so a regression fails loudly.

## Output

| File | Contents |
|---|---|
| `output/04_taxon_sets.rds` | the sets, plus a `tested` universe per analyte |
| `output/04_set_provenance.csv` | per set: source family, evidence, strain and species counts |
| `output/review_sets_summary.tsv` | one row per set, with blank `keep` / `merge_into` / `notes` columns for manual triage |
| `output/review_analytes.tsv` | one row per analyte, for spotting unmerged synonyms |
| `output/review_sets_long.tsv` | full membership, one row per (set, taxon) |
| `output/05_validation.csv` | gold-standard results: pass / fail / no_data / not_covered |
| `output/05_redundant_sets.csv` | set pairs above Jaccard 0.9 — candidate unmerged synonyms |
| `output/TaxSEA_db.rda`, `output/NCBI_ids.rda` | stage 06: TaxSEA's data files with the BacDive sets replaced |
| `output/BacDive_set_provenance.tsv`, `output/BacDive_build_info.tsv` | stage 06: what shipped and how it was built |
| `output/NCBI_ids_additions.tsv`, `output/NCBI_ids_conflicts.tsv` | stage 06: every name added or filled, and every disagreement left alone |

## Legacy

The original single-pass scraper is preserved under `legacy/`, with notes on why
it was replaced.

## Licence

MIT. See `LICENSE`.
