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
  ### options for boundary_source: 
  ### * "global" for entire globe
  ### * bbox in degrees (4-element vector, optional named xmin, xmax, ymin, ymax)
  ### * "rnaturalearth" with region name to pull polygon
  ### * a path to a spatial file readable by st_read()
  src <- cfg$region$boundary_source
  
  ### set up boundary_source = "global" option
  if (identical(src, "global")) {
    ### Return world bounds (entire globe)
    boundary <- sf::st_sf(
      name = "World",
      geometry = sf::st_sfc(
        sf::st_polygon(list(rbind(
          c(-180, -90), c(180, -90), c(180, 90), c(-180, 90), c(-180, -90)
        ))),
        crs = 4326
      )
    )
    return(boundary)
  }

  ### set up for boundary_source = "bbox" in deg xmin/xmax/ymin/ymax
  if (identical(src, "bbox")) {
    ### Create polygon from bounding box
    bbox <- cfg$region$bbox

    ### if bbox vec is named, use names in whatever order given; 
    ### otherwise set names in assumed order
    bbox_names <- c('xmin', 'xmax', 'ymin', 'ymax')
    if(length(bbox) != 4) stop('bbox must include 4 values: xmin, xmax, ymin, ymax in degrees')
    if(any(!names(bbox) %in% bbox_names)) stop('bbox must include 4 values: xmin, xmax, ymin, ymax in degrees')
    if(is.null(names(bbox))) names(bbox) <- bbox_names
    polygon <- sf::st_polygon(list(rbind(
      c(bbox['xmin'], bbox['ymin']),
      c(bbox['xmin'], bbox['ymax']),
      c(bbox['xmax'], bbox['ymax']),
      c(bbox['xmax'], bbox['ymin']),
      c(bbox['xmin'], bbox['ymin'])  ### close ring
    )))
    boundary <- sf::st_sf(name = cfg$region$admin_name, geometry = sf::st_sfc(polygon, crs = 4326))
    return(boundary)
  }

  ### set up boundary_source = "rnaturalearth" and filter by state
  if (identical(src, "rnaturalearth")) {
    b <- rnaturalearth::ne_states(country = "United States of America") |>
      ### update: name %in% cfg$region$admin_name so multiple regions
      ### can be added?
      dplyr::filter(name %in% cfg$region$admin_name) |>
      dplyr::select(name)
  } else {
    b <- sf::st_read(src, quiet = TRUE)
  }
  b_out <- sf::st_transform(b, cfg$crs$working)
  return(b_out)
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
