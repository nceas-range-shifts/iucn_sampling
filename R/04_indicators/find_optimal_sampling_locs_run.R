# =============================================================================
# R/04_indicators/find_optimal_sampling_locs_run.R  —  driver (taxa x years)
# =============================================================================
# Replaces run_find_optimal_samp_locs.R, which reassigned `taxon` and re-
# source()'d the whole script six times. Here it is a plain loop.
#
# Requires config + sampling_config + R/00_utils.R.

source(here::here("R", "04_indicators", "find_optimal_sampling_locs_fn.R"))

run_find_optimal_sampling_locs <- function(taxa = config$taxa, years = config$focal_year,
                                           cfg = config, params = sampling_params) {
  grid <- tidyr::expand_grid(taxon = taxa, year = years)
  purrr::walk2(grid$taxon, grid$year,
               ~find_optimal_sampling_locs(.x, .y, cfg, params, write = TRUE))
  invisible(TRUE)
}
