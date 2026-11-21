### =============================================================================
### R/run_pipeline.R  —  master driver
### =============================================================================
###
### Sources the config + every stage and runs them in order:
###   Stage 0  taxonomy   -> data/processed/species_crosswalk.csv
###   Stage 1  acquire    -> raw GBIF / eBird / IUCN            (API-gated)
###   Stage 2  process    -> hex grid, parsed+cleaned occs, range cells, richness
###   Stage 3  indicators -> sampled cells, undersampled spp, greedy sampling plan
###
### Two entry points:
###   run_pipeline()               full run on YOUR data (needs API keys + inputs)
###   run_pipeline(sample = TRUE)  offline run on the committed sample dataset
###                                (reptiles + mammals; no API calls). Run
###                                data/sample/make_sample_data.R first.
###
### Open case-study-template.Rproj (sets the working directory) before sourcing.

library(here)

### --- Load config + code -------------------------------------------------------
source(here("config", "config.R"))
source(here("config", "taxa_config.R"))
source(here("config", "sampling_config.R"))
source(here("R", "00_utils.R"))

source(here("R", "01_taxonomy",   "01_build_species_crosswalk.R"))
source(here("R", "02_acquire",    "01_download_gbif.R"))
source(here("R", "02_acquire",    "02_download_obis.R"))
source(here("R", "02_acquire",    "03_get_ebird_status.R"))
source(here("R", "02_acquire",    "04_get_redlist_status.R"))
source(here("R", "03_process",    "01_parse_occurrence.R"))
source(here("R", "03_process",    "02_clean_occurrence.R"))
source(here("R", "03_process",    "03_build_hex_grid.R"))
source(here("R", "03_process",    "04_range_from_spatial_run.R"))
source(here("R", "04_indicators", "01_calc_richness_in_grid.R"))
source(here("R", "04_indicators", "02_find_sampled_cells_run.R"))
source(here("R", "04_indicators", "03_find_optimal_sampling_locs_run.R"))
source(here("R", "04_indicators", "04_greedy_setcover_run.R"))

## --- Sample-mode overrides ----------------------------------------------------
## Replace config, taxa_config and sampling config with the sample versions in
## config/sample/. Each file is sourced into a scratch env, then copied over the
## real objects, so the real config stays untouched on disk.
apply_sample_config <- function() {
  sample_env <- new.env()
  for (f in c("config.R", "taxa_config.R", "sampling_config.R")) {
    sys.source(here("config", "sample", f), envir = sample_env)
  }
  for (obj in c("config", "taxa_config", "sampling_params", "make_sampling_config",
                "taxon_output_dir", "taxon_processed_dir", "taxon_year_occ_dir")) {
    assign(obj, get(obj, envir = sample_env), envir = globalenv())
  }
  message("Sample mode: taxa = ", paste(config$taxa, collapse = ", "),
          "; configs from config/sample/.")
}

### --- Orchestration ------------------------------------------------------------
run_pipeline <- function(sample = FALSE) {
  if (isTRUE(sample)) apply_sample_config()

  message("== Stage 0: taxonomy ==")
      
  build_species_crosswalk(force = config$force, progress = 'terse')

  message("== Stage 1: acquire (API-gated) ==")
  system.time({run_gbif_downloads(force = config$force)})
    ### slow, but hopefully should run just as fast across all taxa
    ### as it does with just two taxa (single-shot query)
    ### 1.5 hrs abs and mamms 11/19/26
  system.time({
    ### Only taxa with an obis_dataset_path in taxa_config are fetched from OBIS.
    ### use_aws = TRUE attempts AWS S3 GeoParquet first (fast); falls back to API.
    run_obis_downloads(force = TRUE, match_by = 'both')
  })
  
  message("== Stage 2: process ==")
  system.time({build_hex_grid()}) 
    ### 5 min level 4 - h3o doesn't speed up the conservation polygon intersect step, which is the slow one
    ### * is it line 209: hits <- s2::s2_intersects_matrix(hex_s2, cons_s2)
    ### * maybe can filter out tiny polygons to speed things up?
    ### * but in the grand scheme, this should only need to be run once, so low priority

  run_occurrence_parse() ### 45 s for mammals/abalone

  run_clean_occurrence() ### 5 min for mammals/abalne

  system.time({run_build_ranges()})
    ### using h3o instead of h3jsr: abalone quick; a few seconds
    ### Mammals (more spp, larger ranges): 5 minutes, way better than before

  system.time({run_richness_in_grid()})
    ### 1 minute for abalones and mammals together when not cached

  message("== Stage 3: indicators ==")
  system.time({run_all_find_sampled_cells()})
    ### 70 sec abalones and mammals

  system.time({run_find_optimal_sampling_locs()})
    ### 21 sec abalones and mammals

  system.time({run_greedy_setcover()})
    ### reasonable for abalones; but 30 mintes overall for ab+mamm

  message("== Pipeline complete. Outputs in output/ ==")
  return(invisible(TRUE))
}
