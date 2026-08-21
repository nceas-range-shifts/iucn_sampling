# `R/` — pipeline code

Analysis code for the species-populations / sampling-priority pipeline, organized
by stage. Stages are numbered in execution order and orchestrated by
[`run_pipeline.R`](run_pipeline.R); see [`../docs/pipeline.md`](../docs/pipeline.md)
for the data flow between them.

Many stages come as a **`*_fn.R` / `*_run.R` pair**: `_fn` holds the functions
(the actual work), `_run` is a thin driver that loops those functions over
`taxa` and `years` (defaults from `config/`). Run a whole pipeline with
`run_pipeline()`, or call any `run_*()` / `build_*()` function on its own.

Everything joins on two keys: the **GBIF canonical name** (`name_gbif`, from the
Stage 0 crosswalk) and the **H3 cell id** (`h3_index`).

## Orchestration

| File | Role |
|------|------|
| [`run_pipeline.R`](run_pipeline.R) | Master driver. Sources `config/` + every stage and runs them in order. `run_pipeline()` for a full run on your data; `run_pipeline(sample = TRUE)` for the offline synthetic sample. |
| [`00_utils.R`](00_utils.R) | Small shared helpers (directory creation, empty-file removal, name standardization) used across stages, defined once instead of re-implemented per script. |

## Stage 0 — taxonomy (`01_taxonomy/`)

Runs first. Produces `data/processed/species_crosswalk.csv`, the one
authoritative table linking each species across its input name, GBIF key, eBird
code, and range-shapefile name. **API-gated** (GBIF/eBird).

| File | Role |
|------|------|
| [`resolve_species_names.R`](01_taxonomy/resolve_species_names.R) | The taxonomy engine. Turns a raw species list into a non-destructive crosswalk: `gna_parse` flags malformed names, then the GBIF backbone resolves synonyms/authorities as the authoritative source. |
| [`build_species_crosswalk.R`](01_taxonomy/build_species_crosswalk.R) | Stage 0 driver. Wraps `resolve_species_names()` and adds `range_shapefile_name` / `range_shapefile_path`, reconciling each range shapefile's vendor name to the canonical GBIF name. |

## Stage 1 — acquire (`02_acquire/`)

Downloads raw source data. **All API-gated** — skipped by
`run_pipeline(sample = TRUE)`.

| File | Role |
|------|------|
| [`download_gbif.R`](02_acquire/download_gbif.R) | Submits one filtered GBIF occurrence download per taxon, using the usageKeys already resolved in Stage 0. Needs `GBIF_USER/PWD/EMAIL` in `.Renviron`. |
| [`get_ebird_status.R`](02_acquire/get_ebird_status.R) | Bird range source. Downloads eBird Status & Trends seasonal occurrence products, crops them to the study region, and writes a multi-layer binary presence raster (one layer per species). |
| [`get_redlist_status.R`](02_acquire/get_redlist_status.R) | *Optional.* Looks up each species' IUCN Red List category, keyed off the canonical name. Needs `IUCN_REDLIST_KEY`. |

## Stage 2 — process (`03_process/`)

Turns raw data into the shared spatial layers and per-species tables the
indicators consume.

| File | Role |
|------|------|
| [`build_hex_grid.R`](03_process/build_hex_grid.R) | Builds the two shared spatial layers: every H3 cell covering the study region, and the subset inside conservation areas (with agency / GAP attributes). Region and resolution come from config. |
| [`parse_gbif.R`](03_process/parse_gbif.R) | Subsets each taxon's raw GBIF dataset into one parquet per species-year, filtering on the crosswalk's `name_gbif`. |
| [`clean_gbif.R`](03_process/clean_gbif.R) | Runs CoordinateCleaner + config-driven quality filters and H3 thinning (one record per hex cell) on each per-species-year parquet. |
| [`range_from_shapefiles_fn.R`](03_process/range_from_shapefiles_fn.R) | **Canonical range step.** Converts per-species range shapefiles into a wide table of H3 range cells (one 0/1 column per species), plus per-species totals overall and within conservation areas. |
| [`range_from_raster_fn.R`](03_process/range_from_raster_fn.R) | Bird alternative. Converts the eBird presence raster into the same per-species hex range-cell table the shapefile path produces, so downstream treats all taxa identically. |
| [`range_from_shapefiles_run.R`](03_process/range_from_shapefiles_run.R) | Range-building driver. Dispatches by each taxon's `range_source`: `"shapefile"` → shapefile fn, `"raster"` → raster fn. Both write the same `range_grid_cells.parquet` shape. |
| [`calc_richness_in_grid.R`](04_indicators/calc_richness_in_grid.R) | *(listed under indicators below — richness per hex cell.)* |

## Stage 3 — indicators (`04_indicators/`)

Measures sampling coverage and proposes new sampling locations.

| File | Role |
|------|------|
| [`calc_richness_in_grid.R`](04_indicators/calc_richness_in_grid.R) | Summarizes cleaned occurrences to per-hex-cell species richness and record counts, joined onto the region grid so the result carries geometry. |
| [`find_sampled_cells_fn.R`](04_indicators/find_sampled_cells_fn.R) | Per species: how many range cells are sampled, how many fall inside conservation areas, and in/out-of-range observation counts. |
| [`find_sampled_cells_run.R`](04_indicators/find_sampled_cells_run.R) | Driver: loops the sampled-cells workflow over taxa and years. |
| [`find_optimal_sampling_locs_fn.R`](04_indicators/find_optimal_sampling_locs_fn.R) | Identifies species sampled below `target_sampling_frac` and lists candidate hex cells where extra sampling would close the gap — the inputs to the set-cover. |
| [`find_optimal_sampling_locs_run.R`](04_indicators/find_optimal_sampling_locs_run.R) | Driver over taxa × years. |
| [`greedy_setcover_fn.R`](04_indicators/greedy_setcover_fn.R) | Deterministic greedy weighted set-cover: picks a minimal set of hex cells that brings every undersampled species up to its target, with a conservation-area balance adjustment. Produces the final sampling plan. |
| [`greedy_setcover_run.R`](04_indicators/greedy_setcover_run.R) | Driver: runs the set-cover per taxon (and optionally `"combined"`) for a year. |

## Running a single stage

```r
source("R/run_pipeline.R")     # sources config + all stages
build_species_crosswalk()      # Stage 0
build_hex_grid()               # once per region
run_gbif_parse(taxa = "reptile", years = 2024)
run_clean_gbif(taxa = "reptile", years = 2024)
run_build_ranges(taxa = "reptile")
run_richness_in_grid(taxa = "reptile", years = 2024)
run_all_find_sampled_cells(taxa = "reptile", years = 2024)
run_find_optimal_sampling_locs(taxa = "reptile", years = 2024)
run_greedy_setcover(taxa = "reptile", year = 2024)   # or combined = TRUE
```
