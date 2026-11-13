### =============================================================================
### R/04_indicators/find_optimal_sampling_locs_fn.R  —  undersampled spp + cells
### =============================================================================
###
### Identify species whose regional range is sampled below target_sampling_frac,
### and list candidate hex cells where extra sampling would close the gap. These
### are the inputs the greedy set-cover consumes.
###
### Outputs (per taxon-year, via sampling_config):
###   undersamp_cells_out (parquet): h3_index, species, pres_abs (1 = candidate)
###   undersamp_spp_out   (csv)    : species, cells, n, frac_sampled, cells_needed
###
### Reworked from find_optimal_sampling_locs.R: a function of (taxon, year) instead
### of a 384-line monolith reading globals; the ~220 dead commented lines and the
### stop()-references-undefined-`richness_gpkg` bug are gone; the well-sampled and
### conservation flags are derived inline (no separate annotated-grid file).
###
### Requires config + sampling_config + R/00_utils.R.

suppressPackageStartupMessages({
  library(sf)
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(glue)
  library(here)
})

find_optimal_sampling_locs <- function(taxon, year, cfg = config,
                                       params = sampling_params, write = TRUE) {
  samp_cfg <- make_sampling_config(year, cfg = cfg)[[taxon]]

  ### Well-sampled cells overall (from the richness gpkg for this taxon-year).
  richness <- sf::st_read(samp_cfg$hex_richness_path, quiet = TRUE) |> sf::st_drop_geometry()
  well_sampled_cells <- richness |>
    dplyr::filter(n_records > params$well_sampled_threshold) |>
    dplyr::pull(h3_index)

  ### Range cells (wide) + per-species total range size.
  range_cells <- arrow::read_parquet(samp_cfg$range_cells_path)
  cell_list   <- readr::read_csv(samp_cfg$cell_list_path, show_col_types = FALSE)

  ### Observations per species (from cleaned occurrences).
  occ_dir <- taxon_year_occ_dir(taxon, year, cfg)
  species_files <- file.path(occ_dir,
    paste0(standardize_species_name(cell_list$species), "_clean.parquet"))
  occ_data <- purrr::map_dfr(species_files, function(f) {
    if (file.exists(f)) arrow::read_parquet(f, col_select = c("h3_index", "species", "tot_obs"))
  })

  ### Undersampled species: sampled fraction of range <= target.
  occupied <- dplyr::rename(cell_list, cells = tot_range_cells)
  undersamp_species <- occ_data |>
    dplyr::count(species) |>
    dplyr::right_join(occupied, by = "species") |>
    dplyr::filter(cells != 0) |>
    dplyr::mutate(n = ifelse(is.na(n), 0, n), frac_sampled = n / cells) |>
    dplyr::filter(frac_sampled <= params$target_sampling_frac) |>
    dplyr::mutate(cells_needed = (params$target_sampling_frac * cells) - n) |>
    dplyr::filter(!is.na(species))

  ### Candidate cells: in range, no occurrence yet, not already well-sampled.
  occ_filtered <- dplyr::filter(occ_data, species %in% undersamp_species$species)
  undersamp_cells <- range_cells |>
    dplyr::select(h3_index, dplyr::all_of(undersamp_species$species)) |>
    tidyr::pivot_longer(-h3_index, names_to = "species", values_to = "in_range") |>
    dplyr::left_join(
      occ_filtered |> dplyr::select(h3_index, species) |> dplyr::mutate(has_occ = 1L),
      by = c("h3_index", "species")
    ) |>
    dplyr::mutate(in_range = ifelse(!is.na(has_occ), 0, in_range)) |>
    dplyr::select(-has_occ) |>
    dplyr::filter(!(h3_index %in% well_sampled_cells), in_range == 1) |>
    dplyr::rename(pres_abs = in_range)

  if (isTRUE(write)) {
    ensure_parent_dir(samp_cfg$undersamp_cells_out)
    arrow::write_parquet(undersamp_cells, samp_cfg$undersamp_cells_out)
    readr::write_csv(undersamp_species, samp_cfg$undersamp_spp_out)
    message(glue("[{taxon} {year}] {nrow(undersamp_species)} undersampled species, ",
                 "{dplyr::n_distinct(undersamp_cells$h3_index)} candidate cells"))
  }

  list(undersamp_species = undersamp_species, undersamp_cells = undersamp_cells,
       samp_cfg = samp_cfg)
}
