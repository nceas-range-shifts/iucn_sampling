### =============================================================================
### R/02_acquire/get_ebird_status.R  —  bird range source (eBird Status & Trends)
### =============================================================================
###
### Birds use range_source = "raster": their ranges come from eBird Status &
### Trends seasonal occurrence products rather than range shapefiles. This script
### downloads those products for the crosswalk's bird species, crops them to the
### study region, and writes a multi-layer binary presence raster (one layer per
### species) at taxa_config$bird$ebird_richness_raster.
###
### range_from_raster.R then converts that raster to per-species hex range cells,
### producing the same range_grid_cells.parquet shape as the shapefile path.
###
### API-gated: needs EBIRDST_KEY and EBIRDST_DATA_DIR in .Renviron.
### Reworked from get_ebird_status.R: species come from the crosswalk (ebird_code),
### region/CRS from config; the duplicated for-loop / functional richness blocks
### and the ~/Downloads raster override are dropped.
###
### Requires config + taxa_config + R/00_utils.R + a built crosswalk.

suppressPackageStartupMessages({
  library(ebirdst)
  library(dplyr)
  library(readr)
  library(purrr)
  library(sf)
  library(terra)
  library(rnaturalearth)
  library(glue)
  library(here)
})

### eBird scientific names for this study's birds that have Status & Trends models.
ebird_species_to_process <- function(cfg = config) {
  cw <- readr::read_csv(cfg$paths$crosswalk, show_col_types = FALSE) |>
    dplyr::filter(taxon == "bird", !is.na(ebird_name))
  dplyr::filter(cw, ebird_name %in% ebirdst::ebirdst_runs$scientific_name)$ebird_name
}

### Download one seasonal occurrence product per species (skips failures).
download_ebird_status <- function(species, product = "occurrence_seasonal_mean_3km") {
  purrr::map_lgl(species, function(sp) {
    tryCatch({ ebirdst::ebirdst_download_status(sp, pattern = product); TRUE },
             error = function(e) { message(glue("skip {sp}: {conditionMessage(e)}")); FALSE })
  })
}

### Build a binary presence raster (one layer per species) cropped to the region.
build_bird_presence_raster <- function(cfg = config, taxa_cfg = taxa_config,
                                       product = "occurrence_seasonal_mean_3km",
                                       threshold = 0) {
  species <- ebird_species_to_process(cfg)
  download_ebird_status(species, product)

  region <- rnaturalearth::ne_states(country = "United States of America") |>
    dplyr::filter(name == cfg$region$admin_name) |>
    sf::st_transform(cfg$crs$working) |>
    terra::vect()

  layers <- purrr::map(species, function(sp) {
    r <- tryCatch(ebirdst::load_raster(sp, product = "occurrence", period = "seasonal"),
                  error = function(e) NULL)
    if (is.null(r)) return(NULL)
    r <- terra::app(r, max, na.rm = TRUE)          ### collapse seasons -> max occurrence
    r <- terra::crop(r, terra::project(region, r), mask = TRUE)
    r <- r > threshold                              ### binary presence
    names(r) <- sp
    r
  })
  layers <- purrr::compact(layers)
  presence <- terra::rast(layers)

  out <- taxa_cfg$bird$ebird_richness_raster
  ensure_parent_dir(out)
  terra::writeRaster(presence, out, overwrite = TRUE)
  message(glue("[ebird] {terra::nlyr(presence)} bird presence layers -> {out}"))
  invisible(presence)
}
