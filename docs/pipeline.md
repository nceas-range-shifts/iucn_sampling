# Pipeline stages and data flow

Each stage reads the outputs of the previous one. The **species crosswalk from
Stage 0 is the join key** that ties every data source together, so it runs
first. API-gated steps are marked ⚡ (skipped in `run_pipeline(sample = TRUE)`).

```
                 config/  (region, years, CRS, taxa, thresholds, paths)
                    │
Stage 0  TAXONOMY   ▼
  build_species_crosswalk()  ⚡GBIF/eBird
    → data/processed/species_crosswalk.csv
      one row per species: name_input · name_accepted · name_gbif ·
      gbif_usage_key · ebird_code · range_shapefile_name · needs_review
                    │  (canonical name + keys used by every stage below)
Stage 1  ACQUIRE    ▼
  run_gbif_downloads()        ⚡GBIF   → data/raw/gbif/<taxon>/occurrence.parquet
  get_ebird_status (birds)    ⚡eBird  → bird presence raster
  run_redlist_status()        ⚡IUCN   → optional status table
                    │
Stage 2  PROCESS    ▼
  build_hex_grid()            → hex_grid.gpkg, cons_areas_hex.gpkg
  run_occurrence_parse()            → per-species-year *_occs.parquet
  run_clean_occurrence()            → *_clean.parquet (CoordinateCleaner + H3 thinning)
  run_build_ranges()          → range_grid_cells.parquet  (shapefiles OR raster)
                                 + spp_cells_in_region / _in_cons_area.csv
  run_richness_in_grid()      → hex_species_richness_<year>.gpkg
                    │
Stage 3  INDICATORS ▼
  run_all_find_sampled_cells()        → sampled_cells_<taxon>_<year>.csv
  run_find_optimal_sampling_locs()    → undersamp_cells / undersamp_spp
  run_greedy_setcover()               → greedy_selected_cells.gpkg  ← final plan
```

## The join key

Everything keys on the **GBIF canonical name** (`name_gbif`):

- GBIF occurrences carry `species = name_gbif` (parse filters on it).
- Range shapefiles carry their own `SName`; Stage 0 resolves `SName` → GBIF key
  → `name_gbif` and stores it as `range_shapefile_name`, so range cells and
  occurrence data line up on one canonical name.
- eBird layers are renamed from `ebird_name` → `name_gbif` in `range_from_raster`.

If a species is flagged `needs_review` in the crosswalk (fuzzy GBIF match, no
eBird match, or no range shapefile matched), check it there before trusting
downstream numbers.

## Cell identity

One column, `h3_index`, is the hex-cell id everywhere: the grid, range cells,
cleaned occurrences, and richness all use it, so joins are unambiguous.

## Running stages individually

Every stage exposes a function taking `taxa` / `years` (defaults from config):

```r
source("R/run_pipeline.R")     # sources config + all stages
build_species_crosswalk()      # Stage 0
build_hex_grid()               # once per region
run_occurrence_parse(taxa = "reptile", years = 2024)
run_clean_occurrence(taxa = "reptile", years = 2024)
run_build_ranges(taxa = "reptile")
run_richness_in_grid(taxa = "reptile", years = 2024)
run_all_find_sampled_cells(taxa = "reptile", years = 2024)
run_find_optimal_sampling_locs(taxa = "reptile", years = 2024)
run_greedy_setcover(taxa = "reptile", year = 2024)              # or combined = TRUE
```
