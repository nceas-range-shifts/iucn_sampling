### =============================================================================
### R/04_indicators/calc_richness_in_grid.R  —  species richness per hex cell
### =============================================================================
###
### Summarise cleaned occurrences to per-hex-cell species richness and record
### counts, joined onto the region hex grid so the result carries geometry.
### Written as a GPKG per taxon-year (the sampling step reads it back).
###
### Reworked from calc_gbif_spp_rich_in_grid.R: uses the canonical `h3_index`
### cell column, joins to the grid geometry (the original wrote a geometry-less
### tibble via st_write), and derives paths from sampling_config.
###
### Requires config + sampling_config + R/00_utils.R.

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(purrr)
  library(sf)
  library(here)
})

### One taxon-year: read every *_clean.parquet, summarise, join to grid, write.
calc_richness_in_grid <- function(taxon, year, cfg = config, samp_cfg = NULL) {
  if (is.null(samp_cfg)) samp_cfg <- make_sampling_config(year, cfg = cfg)[[taxon]]
  occ_dir <- taxon_year_occ_dir(taxon, year, cfg)

  files <- list.files(occ_dir, pattern = "_clean\\.parquet$", full.names = TRUE)
  if (length(files) == 0) {
    warning(glue::glue("[{taxon} {year}] no _clean.parquet files in {occ_dir}"))
    return(invisible(NULL))
  }

  occ <- purrr::map_dfr(files, ~arrow::read_parquet(
    .x, col_select = c("h3_index", "species", "tot_obs")))

  summary <- occ |>
    dplyr::filter(!is.na(species), !is.na(h3_index)) |>
    dplyr::group_by(h3_index) |>
    dplyr::summarise(n_species = dplyr::n_distinct(species),
                     n_records = sum(tot_obs, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(h3_index = as.character(h3_index))

  grid <- sf::st_read(cfg$paths$hex_grid, quiet = TRUE)
  richness_sf <- dplyr::left_join(grid, summary, by = "h3_index") |>
    dplyr::mutate(n_species = tidyr::replace_na(n_species, 0L),
                  n_records = tidyr::replace_na(n_records, 0L))

  ensure_parent_dir(samp_cfg$hex_richness_path)
  sf::st_write(richness_sf, samp_cfg$hex_richness_path, delete_dsn = TRUE, quiet = TRUE)
  invisible(richness_sf)
}

run_richness_in_grid <- function(taxa = config$taxa, years = config$years, cfg = config) {
  grid <- tidyr::expand_grid(taxon = taxa, year = years)
  purrr::walk2(grid$taxon, grid$year, ~calc_richness_in_grid(.x, .y, cfg))
  invisible(TRUE)
}
