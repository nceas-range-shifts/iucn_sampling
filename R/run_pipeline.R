### =============================================================================
### R/run_pipeline.R  —  master driver
### =============================================================================
#
### Sources the config + every stage and runs them in order:
###   Stage 0  taxonomy   -> data/processed/species_crosswalk.csv
###   Stage 1  acquire    -> raw GBIF / eBird / IUCN            (API-gated)
###   Stage 2  process    -> hex grid, parsed+cleaned occs, range cells, richness
###   Stage 3  indicators -> sampled cells, undersampled spp, greedy sampling plan
#
### Two entry points:
###   run_pipeline()               full run on YOUR data (needs API keys + inputs)
###   run_pipeline(sample = TRUE)  offline run on the committed sample dataset
###                                (reptiles + mammals; no API calls). Run
###                                data/sample/make_sample_data.R first.
#
### Open case-study-template.Rproj (sets the working directory) before sourcing.

library(here)

### --- Load config + code -------------------------------------------------------
source(here("config", "config.R"))
source(here("config", "taxa_config.R"))
source(here("config", "sampling_config.R"))
source(here("R", "00_utils.R"))

source(here("R", "01_taxonomy", "build_species_crosswalk.R"))
source(here("R", "02_acquire", "download_gbif.R"))
source(here("R", "02_acquire", "get_ebird_status.R"))
source(here("R", "02_acquire", "get_redlist_status.R"))
source(here("R", "03_process", "parse_gbif.R"))
source(here("R", "03_process", "clean_gbif.R"))
source(here("R", "03_process", "build_hex_grid.R"))
source(here("R", "03_process", "range_from_shapefiles_run.R"))
source(here("R", "04_indicators", "calc_richness_in_grid.R"))
source(here("R", "04_indicators", "find_sampled_cells_run.R"))
source(here("R", "04_indicators", "find_optimal_sampling_locs_run.R"))
source(here("R", "04_indicators", "greedy_setcover_run.R"))

### --- Sample-mode overrides ----------------------------------------------------
### Repoint config + taxa_config at data/sample and narrow to the offline taxa.
apply_sample_config <- function() {
  s <- function(...) here("data", "sample", ...)
  config$taxa                    <<- c("reptile", "mammal")
  config$years                   <<- c(2014, 2024)
  config$focal_year              <<- 2024
  config$region$boundary_source  <<- s("spatial_layers", "region.gpkg")
  config$conservation_areas$path <<- s("spatial_layers", "conservation_areas.gpkg")
  config$species_list_default    <<- s("species_list.csv")
  for (t in c("reptile", "mammal")) {
    taxa_config[[t]]$species_source$csv_path <<- s("species_list.csv")
    taxa_config[[t]]$shapefile_dir           <<- s("ranges", t)
    taxa_config[[t]]$gbif_dataset_path       <<- s("gbif", t, "occurrence.parquet")
  }
  message("Sample mode: taxa = reptile, mammal; inputs from data/sample/.")
}

### --- Orchestration ------------------------------------------------------------
run_pipeline <- function(sample = FALSE, do_acquire = !sample) {
  if (isTRUE(sample)) apply_sample_config()

  message("== Stage 0: taxonomy ==")
  if (isTRUE(sample)) {
    ### Name resolution calls the GBIF/eBird APIs, so the sample ships a
    ### pre-built crosswalk; copy it into place instead of rebuilding.
    ensure_parent_dir(config$paths$crosswalk)
    file.copy(here("data", "sample", "species_crosswalk.csv"),
              config$paths$crosswalk, overwrite = TRUE)
    message("Sample mode: using committed data/sample/species_crosswalk.csv")
  } else {
    build_species_crosswalk()
  }

  if (isTRUE(do_acquire)) {
    message("== Stage 1: acquire (API-gated) ==")
    run_gbif_downloads()
    ### get_ebird_status / run_redlist_status are opt-in; call directly if needed.
  } else {
    message("== Stage 1: acquire skipped (using provided inputs) ==")
  }

  message("== Stage 2: process ==")
  build_hex_grid()
  run_gbif_parse()
  run_clean_gbif()
  run_build_ranges()
  run_richness_in_grid()

  message("== Stage 3: indicators ==")
  run_all_find_sampled_cells()
  run_find_optimal_sampling_locs()
  run_greedy_setcover()

  message("== Pipeline complete. Outputs in output/ ==")
  invisible(TRUE)
}
