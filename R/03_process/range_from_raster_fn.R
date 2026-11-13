### =============================================================================
### R/03_process/range_from_raster_fn.R  —  raster range source (birds)
### =============================================================================
###
### Convert the eBird bird presence raster (from get_ebird_status.R) into the same
### per-species hex range-cell table the shapefile path produces, so downstream
### indicators treat all taxa identically.
###
### Reworked from the bird zonal block embedded in get_hex_grid.R (the mangled
### double-assignment there is dropped). Species are keyed by the GBIF canonical
### name (name_gbif) via the crosswalk, so range cells join to occurrence data.
###
### Requires config + taxa_config + sampling_config + R/00_utils.R.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(sf)
  library(terra)
  library(arrow)
  library(here)
})

### Map eBird layer names -> GBIF canonical names using the crosswalk.
.ebird_to_gbif <- function(cfg = config) {
  readr::read_csv(cfg$paths$crosswalk, show_col_types = FALSE) |>
    dplyr::filter(taxon == "bird", !is.na(ebird_name), !is.na(name_gbif)) |>
    dplyr::select(ebird_name, name_gbif)
}

build_range_from_raster <- function(taxon = "bird", cfg = config,
                                    taxa_cfg = taxa_config, samp_cfg = NULL) {
  tcfg <- taxa_cfg[[taxon]]
  if (!identical(tcfg$range_source, "raster")) {
    message(glue::glue("[{taxon}] range_source is not raster; skipped.")); return(invisible(NULL))
  }
  if (is.null(samp_cfg)) samp_cfg <- make_sampling_config(cfg$focal_year, cfg = cfg)[[taxon]]

  presence <- terra::rast(tcfg$ebird_richness_raster)

  grid <- sf::st_read(cfg$paths$hex_grid, quiet = TRUE) |> sf::st_transform(terra::crs(presence))
  grid$zone_id <- seq_len(nrow(grid))
  zone_r <- terra::rasterize(terra::vect(grid), presence, field = "zone_id")

  ### Max presence per hex zone across all species layers.
  zm <- terra::zonal(presence, zone_r, fun = "max", na.rm = TRUE)
  zm <- dplyr::as_tibble(zm)
  ### First column is the zone id; rename to join back to grid.
  names(zm)[1] <- "zone_id"

  wide <- sf::st_drop_geometry(grid) |>
    dplyr::select(h3_index, zone_id) |>
    dplyr::left_join(zm, by = "zone_id") |>
    dplyr::select(-zone_id) |>
    dplyr::mutate(dplyr::across(-h3_index, ~tidyr::replace_na(as.integer(.x > 0), 0L)))

  ### Rename species columns from eBird names to GBIF canonical names.
  lut <- .ebird_to_gbif(cfg)
  ren <- setNames(lut$ebird_name, lut$name_gbif)
  present <- intersect(names(wide), lut$ebird_name)
  if (length(present) > 0) {
    inv <- setNames(lut$name_gbif, lut$ebird_name)
    names(wide)[match(present, names(wide))] <- inv[present]
  }

  ensure_parent_dir(samp_cfg$range_cells_path)
  arrow::write_parquet(wide, samp_cfg$range_cells_path)

  ### Per-species totals (region + conservation areas), same outputs as shapefiles.
  cons_ids <- sf::st_read(cfg$paths$cons_areas_hex, quiet = TRUE)$h3_index
  tot <- wide |> dplyr::summarize(dplyr::across(-h3_index, ~sum(.x, na.rm = TRUE))) |>
    tidyr::pivot_longer(dplyr::everything(), names_to = "species", values_to = "tot_range_cells")
  in_cons <- wide |> dplyr::filter(h3_index %in% cons_ids) |>
    dplyr::summarize(dplyr::across(-h3_index, ~sum(.x, na.rm = TRUE))) |>
    tidyr::pivot_longer(dplyr::everything(), names_to = "species",
                        values_to = "range_cells_in_cons_area")
  out_dir <- ensure_dir(taxon_output_dir(taxon, cfg))
  readr::write_csv(tot, file.path(out_dir, "spp_cells_in_region.csv"))
  readr::write_csv(in_cons, file.path(out_dir, "spp_cells_in_cons_area.csv"))

  invisible(wide)
}
