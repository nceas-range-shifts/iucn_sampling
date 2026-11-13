### =============================================================================
### R/04_indicators/find_sampled_cells_run.R  —  driver (taxa x years)
### =============================================================================
### Thin runner: loops the sampled-cells workflow over taxa and years instead of
### the original's hand-edited spp_group/target_year literals.
###
### Requires config + sampling_config + R/00_utils.R.

source(here::here("R", "04_indicators", "find_sampled_cells_fn.R"))

run_all_find_sampled_cells <- function(taxa = config$taxa, years = config$years,
                                       cfg = config, params = sampling_params) {
  grid <- tidyr::expand_grid(taxon = taxa, year = years)
  purrr::map_dfr(seq_len(nrow(grid)), function(i) {
    run_find_sampled_cells(grid$taxon[i], grid$year[i], cfg, params, write_output = TRUE) |>
      dplyr::mutate(taxon = grid$taxon[i], .before = 1)
  })
}
