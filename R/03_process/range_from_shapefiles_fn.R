### =============================================================================
### R/03_process/range_from_shapefiles_fn.R  —  CANONICAL range step (functions)
### =============================================================================
#
### Turn per-species range shapefiles into a wide table of H3 range cells:
### one row per hex cell, one 0/1 column per species (keyed by the GBIF canonical
### name), plus per-species totals overall and within conservation areas.
#
### Reworked from get_hex_grid_from_shape.R, which overwrote its shapefile list
### three times (so only mammals ever ran) and hardcoded the mammal output path.
### Here it is a function of `taxon`, looped by the _run driver.
#
### The species identity comes from the crosswalk, not the raw shapefile: each
### shapefile's name field (SName) was matched to name_gbif in Stage 0, so range
### cells and occurrence/richness data share ONE canonical key and join cleanly.
#
### Requires config + taxa_config + R/00_utils.R.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(glue)
  library(sf)
  library(h3jsr)
  library(arrow)
  library(here)
})

### Map each range shapefile -> its canonical species name via the crosswalk.
### Returns a tibble: shapefile_path, range_shapefile_name, species (name_gbif).
.shapefile_species_map <- function(taxon, tcfg, cfg = config) {
  cw <- readr::read_csv(cfg$paths$crosswalk, show_col_types = FALSE) |>
    dplyr::filter(taxon == !!taxon, !is.na(range_shapefile_name), !is.na(name_gbif)) |>
    dplyr::transmute(range_shapefile_name, species = name_gbif,
                     shapefile_path = range_shapefile_path)
  if (nrow(cw) == 0) {
    warning(glue("[{taxon}] no shapefile-linked species in crosswalk."))
  }
  cw
}

### Read one range shapefile, transform to the working CRS, return its H3 cells.
.cells_for_shapefile <- function(path, species, cfg = config) {
  range_sf <- sf::st_read(path, quiet = TRUE) |>
    sf::st_transform(cfg$crs$working)
  cells <- h3jsr::polygon_to_cells(range_sf, res = cfg$hex_resolution) |> unlist()
  tibble::tibble(species = species, h3_index = as.character(unique(cells)))
}

### Build the wide range-cell table for one taxon and write it + summaries.
build_range_from_shapefiles <- function(taxon, cfg = config, taxa_cfg = taxa_config,
                                        samp_cfg = NULL) {
  tcfg <- taxa_cfg[[taxon]]
  if (!identical(tcfg$range_source, "shapefile")) {
    message(glue("[{taxon}] range_source is '{tcfg$range_source}', not shapefile; skipped."))
    return(invisible(NULL))
  }
  if (is.null(samp_cfg)) samp_cfg <- make_sampling_config(cfg$focal_year, cfg = cfg)[[taxon]]

  ### Full CA hex grid and conservation-area cells (from build_hex_grid.R).
  hex_grid <- sf::st_read(cfg$paths$hex_grid, quiet = TRUE) |> sf::st_drop_geometry()
  cons_ids <- sf::st_read(cfg$paths$cons_areas_hex, quiet = TRUE)$h3_index

  spp_map <- .shapefile_species_map(taxon, tcfg, cfg)

  long <- purrr::pmap_dfr(
    list(spp_map$shapefile_path, spp_map$species),
    function(path, species) tryCatch(
      .cells_for_shapefile(path, species, cfg),
      error = function(e) { warning(glue("failed {basename(path)}: {conditionMessage(e)}")); NULL }
    )
  )

  wide <- long |>
    dplyr::mutate(presence = 1L) |>
    dplyr::distinct() |>
    tidyr::pivot_wider(names_from = species, values_from = presence, values_fill = 0L)

  ### Join onto the full grid so every cell is represented.
  range_cells <- hex_grid |>
    dplyr::left_join(wide, by = "h3_index") |>
    dplyr::mutate(dplyr::across(-h3_index, ~tidyr::replace_na(., 0L)))

  ensure_parent_dir(samp_cfg$range_cells_path)
  arrow::write_parquet(range_cells, samp_cfg$range_cells_path)

  ### Per-species totals: whole region and within conservation areas.
  tot <- range_cells |>
    dplyr::summarize(dplyr::across(-h3_index, ~sum(.x, na.rm = TRUE))) |>
    tidyr::pivot_longer(dplyr::everything(),
                        names_to = "species", values_to = "tot_range_cells")
  in_cons <- range_cells |>
    dplyr::filter(h3_index %in% cons_ids) |>
    dplyr::summarize(dplyr::across(-h3_index, ~sum(.x, na.rm = TRUE))) |>
    tidyr::pivot_longer(dplyr::everything(),
                        names_to = "species", values_to = "range_cells_in_cons_area")

  out_dir <- taxon_output_dir(taxon, cfg)
  ensure_dir(out_dir)
  readr::write_csv(tot, file.path(out_dir, "spp_cells_in_region.csv"))
  readr::write_csv(in_cons, file.path(out_dir, "spp_cells_in_cons_area.csv"))

  invisible(range_cells)
}
