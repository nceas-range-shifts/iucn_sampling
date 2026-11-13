### =============================================================================
### config/sampling_config.R  —  sampling-optimization configuration
### =============================================================================
###
### Replaces the original project's four near-identical files (sampling_config.R,
### _2014.R, _2023.R, _cell_threshold.R) that differed only by a hardcoded year
### and an output-filename suffix. Here it is ONE builder parameterized by year.
###
### Requires config/config.R and config/taxa_config.R sourced first.
###
### Usage:
###   cfg    <- make_sampling_config(config$focal_year)
###   cfg    <- make_sampling_config(2014, suffix = "cell_thresh")  # variant run
###   params <- sampling_params

### Tuning knobs — the single source of truth (the originals scattered these and
### even hardcoded `50` as a literal inside functions).
sampling_params <- list(
  well_sampled_threshold = 50,    ### min records per cell to count as "well-sampled"
  target_sampling_frac   = 0.10,  ### coverage target: 10% of each species' range
  cons_balance_weight    = 0.20,  ### greedy: weight toward conservation areas
  batch_fraction         = 0.35,  ### greedy: cells added per batch
  deficit_percentile     = 0.15,  ### greedy: deficit cut-off per batch
  hex_resolution         = config$hex_resolution
)

### Build the per-taxon path set for one year. `suffix` (optional) is inserted
### into the undersampled-output filenames so an alternate run (e.g. a
### cell-threshold variant) does not overwrite the main outputs.
make_sampling_config <- function(year, suffix = "", cfg = config,
                                 taxa_cfg = taxa_config) {
  year <- as.character(year)
  tag  <- if (nzchar(suffix)) paste0("_", suffix) else ""

  build_one <- function(taxon) {
    tcfg     <- taxa_cfg[[taxon]]
    proc_dir <- taxon_processed_dir(taxon, cfg)
    out_dir  <- taxon_output_dir(taxon, cfg)

    list(
      taxon        = taxon,
      range_source = tcfg$range_source,
      ### parsed/cleaned occurrences for this taxon-year
      occ_dir      = taxon_year_occ_dir(taxon, year, cfg),
      ### per-species range -> hex cells (from shapefiles or the eBird raster).
      ### Ranges are static across years, so this path is not year-keyed.
      range_cells_path = file.path(proc_dir, "range_grid_cells.parquet"),
      ### richness + observation counts per hex cell
      hex_richness_path = file.path(out_dir,
                                    paste0("hex_species_richness_", year, ".gpkg")),
      ### per-species total range cells (static; written by the range step)
      cell_list_path    = file.path(out_dir, "spp_cells_in_region.csv"),
      ### sampling-optimization outputs
      output_dir        = file.path(out_dir, "sampling_optimized"),
      undersamp_cells_out = file.path(proc_dir,
                                      paste0("undersamp_cells_", year, tag, ".parquet")),
      undersamp_spp_out   = file.path(proc_dir,
                                      paste0("undersamp_spp_", year, tag, ".csv"))
    )
  }

  stats::setNames(lapply(cfg$taxa, build_one), cfg$taxa)
}
