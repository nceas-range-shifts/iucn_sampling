# CHANGELOG

## 2026-11-17 21:10 UTC — Consolidate OBIS download implementations

**Archived:** `R/02_acquire/02_download_obis_v2.R` and `02_download_obis_original_api_only.R` → `R/02_acquire/_archive/`

**Reason:** Consolidated to single canonical `02_download_obis.R` which supports both AWS S3 GeoParquet (primary) and API fallback. The archived versions were intermediate prototypes; the current version subsumes their functionality.


## 2026-11-18 PST — Switch range_from_polygons to h3o

**Archived:** `R/03_process/04a_range_from_polygons_fn.R` (h3jsr version) → `R/03_process/_archive/04a_range_from_polygons_h3jsr_fn.R`
**Archived:** `R/03_process/04c_range_from_polygons_h3o_fn.R` (experimental prototype) → `R/03_process/_archive/`
**Edited:** `R/03_process/04a_range_from_polygons_fn.R` now uses `h3o::sfc_to_cells()`; function names unchanged, so `04_range_from_spatial_run.R` needed no edit.
**Edited:** `R/03_process/05_benchmark_h3o_vs_h3jsr.R` (h3o call fixed to `sfc_to_cells`), `R/03_process/README_h3o_benchmark.md` (rewritten with results).
**Note:** range manifests hash the range code, so all polygon taxa will rebuild on the next run.

## 2026-11-18 PST — Switch build_hex_grid to h3o; add hex-grid benchmark

**Edited:** `R/03_process/03_build_hex_grid.R` — h3jsr calls (`get_children`/`get_res0`, `polygon_to_cells`, `cell_to_polygon`) replaced with h3o (`get_children`, `sfc_to_cells`, `sf::st_as_sfc()` on H3 vectors).
**Added:** `R/03_process/06_benchmark_hex_grid_h3o_vs_h3jsr.R` → `logs/benchmark_hex_grid_h3o_vs_h3jsr.csv`.
**Note:** the hex-grid manifest hashes inputs, not code, so existing grid layers stay cached; use `build_hex_grid(force = TRUE)` to regenerate with h3o (not yet done).

## 2026-11-18 09:37 PST — Add caching and progress option to richness-in-grid

**Edited:** `R/04_indicators/01_calc_richness_in_grid.R` — `calc_richness_in_grid()` / `run_richness_in_grid()` gain `force` and `progress` ("terse"/"verbose"/"quiet") args. Per taxon-year manifest `<taxon output dir>/logs/richness_manifest_<taxon>_<year>.json` (hash of clean-parquet md5s, hex-grid md5, function code). Both functions now return a status tibble instead of the sf / `TRUE`.
**Edited:** `R/run_pipeline.R` — call is now `run_richness_in_grid(force = FALSE)`; stale TODO note removed.
**Edited:** `R/README.md` — richness row notes caching/progress.
**Note:** no existing manifests, so the first run after this change recomputes every taxon-year.

## 2026-11-18 10:15 PST — Add caching and progress option to sampled-cells, optimal-locs and greedy set-cover

**Edited:** `R/00_utils.R` — new shared helpers `.stage_cached()` and `.write_stage_manifest()` (manifest = settings hash; reusable only if hash matches and outputs exist; manifest written last).
**Edited:** `R/04_indicators/02a_find_sampled_cells_fn.R`, `02_find_sampled_cells_run.R`, `03a_find_optimal_sampling_locs_fn.R`, `03_find_optimal_sampling_locs_run.R`, `04a_greedy_setcover_fn.R`, `04_greedy_setcover_run.R` — `force` and `progress` ("terse"/"verbose"/"quiet") args. Manifests: `<taxon output dir>/logs/{sampled_cells,optimal_locs}_manifest_<taxon>_<year>.json`; `<greedy output dir>/greedy_manifest.json`. Hashes cover upstream input md5s, relevant `sampling_params`, and function code.
**Edited:** `R/run_pipeline.R` — stage 3 calls pass `force = FALSE`.
**Note:** `find_optimal_sampling_locs()` and `greedy_setcover()` now return a `status` element ("written"/"cached"); cached runs return only `status` plus paths, not the in-memory tables. No existing manifests, so the first run recomputes everything.

## 2026-11-18 10:40 PST — `force` defaults to `config$force` in all stage functions

**Edited:** every `force = FALSE` default in `R/01_taxonomy`, `R/02_acquire`, `R/03_process`, `R/04_indicators` (not the benchmark helpers or internal `.stage_cached()`) is now `force = config$force`. Explicit `force = ...` arguments in `R/run_pipeline.R` still override it.

## 2026-11-18 — Sample mode: point conservation-area config at sample gpkg

**Edited:** `R/run_pipeline.R` — `apply_sample_config()` now clears `conservation_areas$layer` (the WDPA gdb layer name `WDPA_poly_Nov2026` was being passed to the single-layer sample gpkg, causing "Cannot open layer") and overrides the field names and values (`gap_field`, `gap_values`, `realm_*`, `unit_field`, `agency_*`) to match `data/sample/spatial_layers/conservation_areas.gpkg`.
**Note:** Field mapping inferred from the sample file's attributes; not yet run end-to-end.

## 2026-11-18 — Sample config split into data/sample/config.R

**Added:** `data/sample/config.R` — standalone `config` (plus taxon path helpers) for the offline sample: sample conservation-area schema, sample paths under `data/sample/`, `taxa = c("reptile","mammal")`, OBIS off.
**Edited:** `R/run_pipeline.R` — `apply_sample_config()` now `sys.source()`s `data/sample/config.R` and overwrites the global `config` and taxon helpers, instead of patching individual keys of the real config. Per-taxon `taxa_config` path overrides retained.
**Note:** Stage-0 crosswalk copy and the `taxa_config` path overrides are unchanged. Not yet run end-to-end.

## 2026-11-18 — Sample configs moved to config/sample/

**Moved:** `data/sample/config.R` removed; replaced by `config/sample/config.R`.
**Added:** `config/sample/taxa_config.R` (polygon taxa `reptile`, `mammal`; `polygon_name_field = "SName"` from the range shapefiles; taxa taken from `data/sample/species_crosswalk.csv`) and `config/sample/sampling_config.R` (same tuning knobs and `make_sampling_config` as the full run).
**Edited:** `R/run_pipeline.R` — `apply_sample_config()` sources all three `config/sample/` files into a scratch env and assigns `config`, `taxa_config`, `sampling_params`, `make_sampling_config` and the `taxon_*` helpers into the global env.
**Note:** Not yet run end-to-end.

## 2026-11-18 — Sample config uses GBIF download API

**Edited:** `config/sample/config.R` — added `gbif_source = 'download'` so sample runs exercise the GBIF download API path (requires `GBIF_USER`, `GBIF_PWD`, `GBIF_EMAIL`).
**Note:** Stage 1 is skipped when `run_pipeline(sample = TRUE)` is called, so this setting only takes effect if acquire is run explicitly for sample data.

## 2026-11-18 — Guard empty undersampled-species set in find_optimal_sampling_locs

**Edited:** `R/04_indicators/03a_find_optimal_sampling_locs_fn.R` — when no species are undersampled, return an empty candidate table (`h3_index`, `species`, `pres_abs`) instead of calling `tidyr::pivot_longer(-h3_index)` on an empty column selection, which errored with "`cols` must select at least one column".
**Note:** Confirmed on sample data: 0 undersampled species for reptile and mammal at focal year 2024. Not yet rerun end-to-end.

## 2026-11-18 — Pin undersamp_species column types in greedy set-cover

**Edited:** `R/04_indicators/04a_greedy_setcover_fn.R` — `read_csv()` of `undersamp_spp` now sets explicit column types. An empty CSV (no undersampled species) was read back with every column as character, so `cells_needed - final_deficit` errored with "non-numeric argument to binary operator".
**Note:** Cause confirmed on sample data: the on-disk CSV has 0 rows and all columns read as chr. The in-memory table was numeric, so the earlier type check missed it.

## 2026-11-18 — Empty-selection output in greedy set-cover

**Edited:** `R/04_indicators/04a_greedy_setcover_fn.R` — when no batches are selected, `selected_cells` is now an empty table with the expected columns (`order`, `h3_index`, `cons_area`, `n_spp_covered`, `batch_id`, `spp_remaining`) instead of a column-less table. Previously the later `left_join(by = "h3_index")` failed with "Join columns in `y` must be present in the data".
**Note:** Sample run selects 0 cells (no undersampled species), so this branch is exercised but the greedy algorithm itself is not tested by the sample.

## 2026-11-18 15:23 PST — Ignore sample-run outputs

**Edited:** `.gitignore` — added `data/sample/output/` and `data/sample/processed/` under a new section after the sample force-include block. Both are regenerated by `run_pipeline(sample = TRUE)`; neither was tracked.

## 2026-11-18 15:24 PST — Deleted untracked GBIF pending-request records

**Deleted:** `data/sample/gbif/mammal/gbif_pending.json` (download_key 0000626-261118162104854) and `data/sample/gbif/reptile/gbif_pending.json` (download_key 0000625-261118162104854). Untracked; created 2026-11-18 23:02 UTC by a sample-mode GBIF submission.

## 2026-11-20 11:01 PST — Crosswalk stage: progress messages and manifest cache

**Edited:** `R/01_taxonomy/01_build_species_crosswalk.R` — added `progress = "terse"|"verbose"|"quiet"` (via `.progress_msg()`) with per-taxon GBIF, polygon and WoRMS match counts; `attach_worms_aphia()` now requires a `taxon` argument; `build_taxon_crosswalk()` takes an optional `raw` list; manifest cache at `species_crosswalk_manifest.json` beside the crosswalk CSV.

## 2026-11-20 11:01 PST — OBIS download: combined S3 query, AphiaID-or-name match

**Edited:** `R/02_acquire/02_download_obis.R` — `run_obis_downloads()` rewritten (plan/cache, one combined S3 query split by taxon, per-taxon API fallback); new helpers `get_obis_aphia_ids()`, `prepare_obis_taxon()`, `obis_summary_row()`, `obis_error_row()`, `write_obis_taxon()`; `fetch_obis_aws(aphia_ids, species_names, ...)` signature changed; `get_obis_names()` drops species with no `aphia_id`; `obis_request_hash()` takes `aphia_ids`; `download_one_taxon_obis()` is now a wrapper over `run_obis_downloads()`; new `batch_size` argument. Result tibble gains `n_no_aphia`.
**Backup:** pre-refactor version at `/tmp/pa-31947/7a73bc5a-0795-4d14-ab5c-f68180cfbf7c/02_download_obis.R.bak` (scratch, not in repo).
**Note:** test outputs went to scratch only; `data/raw/obis/` not touched.

## 2026-11-20 18:08 PST — OBIS download: match_by argument

**Edited:** `R/02_acquire/02_download_obis.R` — new `match_by = c("both", "name", "id")` on `run_obis_downloads()` (default "both"), passed through `download_one_taxon_obis()`, `prepare_obis_taxon()`, `obis_request_hash()` and `fetch_obis_aws()`; header usage examples updated. API fallback warns when `match_by = "id"` (it queries by name only).

## 2026-11-21 PST — Range table tail rewritten (matrix/tabulate)

**Edited:** `R/03_process/04a_range_from_polygons_fn.R` (wide table, parquet frame and totals via `match` + matrix + `tabulate`); `R/03_process/04b_range_from_raster_fn.R` (totals via `colSums`). No files renamed or archived.
**Note:** abalones outputs were regenerated during testing (values unchanged); backup of the previous outputs at `/tmp/pa-31947/b86c5e80-4aef-419d-94e0-a999dda9daa2/abalones_bak/` (scratch, not in repo). Range manifests are now stale for all polygon and raster taxa because the hashed range code changed: the next `run_build_ranges()` rebuilds them.
