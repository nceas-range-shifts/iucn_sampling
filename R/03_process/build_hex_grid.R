# =============================================================================
# R/03_process/build_hex_grid.R  —  region H3 grid + conservation-area cells
# =============================================================================
#
# Produces the two shared spatial layers the rest of the pipeline joins on:
#   config$paths$hex_grid        — every H3 cell covering the study region
#   config$paths$cons_areas_hex  — the subset of cells inside conservation areas,
#                                  with managing-agency / GAP attributes.
#
# Reworked from get_hex_grid.R: region and resolution come from config; the
# conservation-layer field names are configurable; the bird eBird-raster zonal
# block (a range source, not a grid step) is moved to range_from_raster.R; and
# the mangled double-assignment bug at get_hex_grid.R:144-147 is dropped.
#
# Requires config + R/00_utils.R.  Usage:  build_hex_grid()

suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(h3jsr)
  library(rnaturalearth)
  library(here)
})

# Load the study-region boundary as an sf polygon in the working CRS.
get_region_boundary <- function(cfg = config) {
  src <- cfg$region$boundary_source
  if (identical(src, "rnaturalearth")) {
    b <- rnaturalearth::ne_states(country = "United States of America") |>
      dplyr::filter(name == cfg$region$admin_name) |>
      dplyr::select(name)
  } else {
    b <- sf::st_read(src, quiet = TRUE)
  }
  sf::st_transform(b, cfg$crs$working)
}

# Load and filter the conservation-areas layer to strictly-protected terrestrial
# units, dissolving to one geometry per management unit.
get_conservation_areas <- function(cfg = config) {
  ca <- cfg$conservation_areas
  sfp <- sf::st_read(ca$path, quiet = TRUE) |> sf::st_zm(drop = TRUE)

  keep <- sfp[[ca$gap_field]] %in% ca$gap_values
  if (!is.null(ca$realm_field)) keep <- keep & sfp[[ca$realm_field]] == ca$realm_value
  sfp <- sfp[keep, ]

  sfp |>
    dplyr::group_by(dplyr::across(dplyr::all_of(
      c(ca$unit_field, ca$gap_field, ca$agency_field, ca$agency_level_field)))) |>
    dplyr::summarize(.groups = "drop") |>
    sf::st_transform(cfg$crs$working)
}

build_hex_grid <- function(cfg = config) {
  region <- get_region_boundary(cfg)

  # Region H3 grid.
  ids <- h3jsr::polygon_to_cells(region, res = cfg$hex_resolution, simple = TRUE)
  hex_grid <- h3jsr::cell_to_polygon(ids, simple = FALSE) |>
    dplyr::rename(h3_index = h3_address)
  ensure_parent_dir(cfg$paths$hex_grid)
  sf::st_write(hex_grid, cfg$paths$hex_grid, delete_dsn = TRUE, quiet = TRUE)

  # Conservation-area cells (attributes summarised per cell).
  ca_cfg <- cfg$conservation_areas
  cons <- get_conservation_areas(cfg)

  old_s2 <- sf::sf_use_s2(); sf::sf_use_s2(FALSE); on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  cons_cells <- hex_grid |>
    sf::st_join(cons, join = sf::st_intersects) |>
    dplyr::filter(!is.na(.data[[ca_cfg$unit_field]])) |>
    dplyr::group_by(h3_index) |>
    dplyr::summarize(
      unit_name  = paste(unique(na.omit(.data[[ca_cfg$unit_field]])), collapse = "; "),
      mng_agency = paste(unique(na.omit(.data[[ca_cfg$agency_field]])), collapse = "; "),
      gap_status = paste(sort(unique(na.omit(.data[[ca_cfg$gap_field]]))), collapse = ", "),
      do_union = FALSE, .groups = "drop"
    )
  sf::st_write(cons_cells, cfg$paths$cons_areas_hex, delete_dsn = TRUE, quiet = TRUE)

  message(glue::glue("[hex_grid] {nrow(hex_grid)} region cells, ",
                     "{nrow(cons_cells)} conservation cells written."))
  invisible(list(hex_grid = hex_grid, cons_cells = cons_cells))
}
