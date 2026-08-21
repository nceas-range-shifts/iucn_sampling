# Adapting the pipeline to your own data

You should only ever need to **edit config and swap input files** — not the
analysis code in `R/`. This page lists every config knob and the exact shape of
each input file.

## 1. Config knobs

### `config/config.R` (global)
| Field | What to set |
|---|---|
| `region$country_code`, `region$admin_name`, `region$admin_level` | Your study region. `admin_name` scopes GBIF/eBird and selects the boundary. |
| `region$boundary_source` | `"rnaturalearth"` to pull the admin polygon, or a path to a boundary shp/gpkg. |
| `years`, `focal_year`, `compare_years` | Years to process; the single snapshot year; the two years compared. |
| `crs$working / equal_area / display` | Working CRS is lon/lat for H3; equal_area for stats; display for maps. |
| `hex_resolution` | H3 resolution (8 ≈ 0.7 km² cells). |
| `gbif_filters` | Occurrence-download predicates (basis of record, etc.). |
| `cleaning` | Coordinate-quality thresholds. |
| `conservation_areas` | Path to your protected-areas layer **and the field names** in it (GAP status, realm, unit, agency). |
| `taxa` | Which taxa to run (must match `taxa_config`). |

### `config/taxa_config.R` (per taxon)
Add or edit one entry per taxon:
- `species_source` — `csv_path` (or `sheet_url`/`sheet_name`) and `sci_name_col`.
  `taxon_filter` filters a shared multi-taxon list by a `taxon` column.
- `range_source` — `"shapefile"` (per-species range shapefiles) or `"raster"`
  (eBird Status & Trends abundance, birds).
- `shapefile_dir` + `shapefile_name_field` — folder of range shapefiles and the
  attribute holding the species name (default `"SName"`).
- `gbif_dataset_path` — where the raw GBIF download for this taxon lands.
- `gbif_columns_keep` — shared; keep as-is unless your GBIF export differs.

### `config/sampling_config.R`
`sampling_params` holds the thresholds: `well_sampled_threshold` (records/cell),
`target_sampling_frac` (coverage target), and the greedy knobs. Usually you only
touch these to tune the analysis.

## 2. Input file shapes

### Species list (CSV)
At minimum a scientific-name column (name it in `sci_name_col`). A shared list
may include a `taxon` column so one file serves all taxa.
```
scientific_name,taxon
Sceloporus occidentalis,reptile
Sylvilagus bachmani,mammal
```

### Range shapefiles (one per species)
A polygon layer per species in `shapefile_dir`, each with a name attribute
(`shapefile_name_field`, e.g. `SName`) holding that species' name. Names need
not match GBIF exactly — Stage 0 resolves them to the GBIF key and records the
link in the crosswalk.

### Conservation-areas layer (shp/gpkg)
A protected-areas polygon layer with the fields named in
`config$conservation_areas` (status, realm, unit, agency). Only strictly
protected terrestrial units (by your `gap_values`/`realm_value`) are kept.

### GBIF occurrences
Normally downloaded by `run_gbif_downloads()` (needs GBIF keys). To supply your
own, place a parquet at each taxon's `gbif_dataset_path` with at least the
columns in `gbif_columns_keep` (notably `species`, `year`,
`decimalLatitude/Longitude`, `basisOfRecord`, `individualCount`).

### eBird (birds, `range_source = "raster"`)
`get_ebird_status.R` downloads eBird Status & Trends products (needs
`EBIRDST_KEY`) and writes a presence raster; `range_from_raster.R` converts it
to range cells. Or point `taxa_config$bird$ebird_richness_raster` at your own
multi-layer presence raster.

## 3. Checklist to move to a new study
1. Replace the species list, range shapefiles, and conservation-areas layer.
2. Set `region`, `years`, and `crs` in `config/config.R`.
3. Set each taxon's `species_source`, `range_source`, and `shapefile_dir`.
4. Add API keys to `.Renviron`.
5. `run_pipeline()` — or run stage functions individually.
