# Species Populations Pipeline (template)

[![License: MIT](https://img.shields.io/github/license/king0708/biodiversity-sampling-indicators)](LICENSE)

A reusable R pipeline for building **species-population and sampling-priority
indicators** from biodiversity data. Given a species list and their range
shapefiles (plus, optionally, GBIF / eBird / IUCN access), it:

1. resolves every species name into one authoritative crosswalk,
2. downloads and cleans GBIF occurrence records,
3. turns species ranges into a hexagonal (H3) grid,
4. measures how well each species' range is currently sampled, and
5. runs a greedy set-cover to propose a **minimal set of new sampling
   locations** that brings undersampled species up to a coverage target.

It is a template: adapt it to a new region, taxa, or set of species **by editing
config files and swapping input files — not by editing analysis code**. A small
committed sample dataset lets you run the whole thing offline first.

> Adapted from a California biodiversity case study. The defaults (California,
> H3 resolution 8, California-Albers CRS) are just defaults — all live in
> `config/`.

---

## Repository layout

```
config/     config.R (region, years, CRS, paths), taxa_config.R, sampling_config.R
R/
  00_utils.R
  01_taxonomy/    resolve_species_names.R, build_species_crosswalk.R   <- runs FIRST
  02_acquire/     download_gbif.R, get_ebird_status.R, get_redlist_status.R
  03_process/     parse_gbif.R, clean_gbif.R, build_hex_grid.R,
                  range_from_shapefiles_fn/_run.R, range_from_raster_fn.R
  04_indicators/  calc_richness_in_grid.R, find_sampled_cells_*,
                  find_optimal_sampling_locs_*, greedy_setcover_*
  run_pipeline.R  <- orchestrates all stages
data/     sample/ (committed demo inputs), raw/ + processed/ (git-ignored)
output/   indicators + sampling plans (git-ignored)
docs/     adapting-to-your-data.md, pipeline.md
```

Name resolution is **Stage 0** on purpose: the crosswalk it writes
(`data/processed/species_crosswalk.csv`) is the single join key linking each
species across its input name, GBIF key, eBird code, and range-shapefile name.
Every later stage joins on it. See `docs/pipeline.md`.

---

## Prerequisites

- **R ≥ 4.4**. Open `case-study-template.Rproj` so the working directory is the
  project root (paths use `here::here()`).
- **System libraries** for the geospatial/columnar stack: GDAL, GEOS, PROJ
  (for `sf`/`terra`) and a working `V8` (for `h3jsr`).
- **Packages**: `renv::restore()` (preferred, pinned) or `Rscript install.R`.
- **API accounts** (only for a full run on your own data — the sample needs none):
  - GBIF — https://www.gbif.org/user/profile
  - eBird Status & Trends — https://science.ebird.org/en/status-and-trends/download-data
  - IUCN Red List (optional) — https://apiv3.iucnredlist.org/api/v3/token

Copy the environment template and fill in your keys:

```bash
cp .Renviron.example .Renviron   # then edit .Renviron and restart R
```

`.Renviron` is git-ignored — never commit real keys.

---

## Quick start — run the sample (offline, no keys)

```r
# from the project root (open the .Rproj first)
source("data/sample/make_sample_data.R")     # generate the tiny demo inputs
source("R/run_pipeline.R")
run_pipeline(sample = TRUE)
```

This runs reptiles + mammals end-to-end with no network access and writes
results under `output/indicators/species_populations/`, ending in a greedy
sampling plan (`.../sampling_optimized/<year>/greedy_selected_cells.gpkg`).

## Run on your own data

1. Put your inputs in place (see `docs/adapting-to-your-data.md` for exact file
   shapes): a species-list CSV, per-species range shapefiles, and a
   conservation-areas layer.
2. Edit `config/config.R` (region, years, CRS) and `config/taxa_config.R`
   (per-taxon species source, `range_source`, shapefile dir).
3. Add your API keys to `.Renviron`.
4. Run:

```r
source("R/run_pipeline.R")
run_pipeline()            # full run: resolves names, downloads GBIF, ...
```

Stages can also be run individually — every `*_run.R` / `run_*()` function takes
`taxa` and `years` arguments. See `docs/pipeline.md`.

---

## What comes out

Per taxon and year, under `output/indicators/species_populations/<taxon>/`:

- `hex_species_richness_<year>.gpkg` — richness + record counts per hex cell
- `spp_cells_in_region.csv`, `spp_cells_in_cons_area.csv` — range size summaries
- `sampling_optimized/sampled_cells_<taxon>_<year>.csv` — per-species sampling status
- `sampling_optimized/<year>/greedy_selected_cells.gpkg` — **the proposed new
  sampling locations**, plus `greedy_species_summary.csv`

Run logs (per-species success/failure) are written under `.../logs/`.

## Report

`reports/species_populations_report.qmd` is a small Quarto report that reads the
pipeline outputs and shows (1) how many species meet their sampling target, by
taxon, and (2) a sampling-priority map (a kernel-density surface of the selected
cells over the study region). Render it after a pipeline run:

```bash
quarto render reports/species_populations_report.qmd
```

It uses `config$focal_year` and whichever taxa have outputs, so it works on the
sample run as-is. The rendered HTML is git-ignored (regenerate as needed).

## Data sources & citation

The committed sample dataset under [`data/sample/`](data/sample/) is **synthetic**
(see [`data/sample/README.md`](data/sample/README.md)) and must **not** be cited as
real data.

When you run the pipeline on **real** inputs, cite the sources you actually use.
The pipeline is designed to consume:

- **Species range maps — CDFW CWHR.** The original case study used the California
  Wildlife Habitat Relationships (CWHR) System Species Distribution Maps. Per the
  dataset's own metadata, cite as:

  > California Department of Fish and Wildlife, California Interagency Wildlife
  > Task Group. 2014. *California Wildlife Habitat Relationships (CWHR) System
  > Species Distribution Maps*, CWHR version 9.0. Biogeographic Data Branch,
  > Sacramento, CA. Accessed \<date\>.
  > https://www.wildlife.ca.gov/Data/CWHR

  (Each shapefile carries the species' scientific name in the `SName` field, which
  Stage 0 reconciles to the GBIF canonical name.)

- **Occurrence records — GBIF.** Every GBIF download issues its own DOI; cite that
  specific download:

  > GBIF.org (\<DD Month YYYY\>) *GBIF Occurrence Download*
  > https://doi.org/10.15468/dl.\<XXXXXX\>

- **Bird ranges — eBird Status & Trends.** Use the versioned citation shown on the
  [eBird Status & Trends download page](https://science.ebird.org/en/status-and-trends)
  (Fink et al., Cornell Lab of Ornithology; `https://doi.org/10.2173/ebirdst.<version>`),
  matching the data version you downloaded.

- **Conservation status — IUCN Red List** (optional):

  > IUCN (\<YYYY\>). *The IUCN Red List of Threatened Species.* Version \<YYYY-N\>.
  > https://www.iucnredlist.org

Fill in the bracketed access dates, DOIs, and versions for your own run.

## License

MIT — see [LICENSE](LICENSE).
