# =============================================================================
# R/03_process/clean_gbif.R  —  clean + spatially thin GBIF occurrences
# =============================================================================
#
# For each per-species-year parquet from parse_gbif.R: run CoordinateCleaner,
# apply quality filters (thresholds from config$cleaning), diagnostic checks,
# and H3 thinning (one record per hex cell at config$hex_resolution). Writes
# `<species>_clean.parquet` next to the input and a per-taxon-year run log.
#
# Reworked from the original clean_gbif_occs.R: thresholds and CRS/resolution
# come from config; the ~13 hand-toggled input/output pairs are replaced by a
# loop over taxa x years; every species' outcome is recorded in a run log
# instead of a bare cat()/NULL.
#
# Requires config + R/00_utils.R.  Usage:  res <- run_clean_gbif()

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(glue)
  library(arrow)
  library(sf)
  library(CoordinateCleaner)
  library(h3jsr)
  library(here)
})

# GBIF SIMPLE_PARQUET uses lower-case column names; normalise the ones we use to
# the canonical camelCase, renaming only those actually present (the original
# rename() errored if any were missing).
.normalise_gbif_cols <- function(df) {
  map <- c(
    gbifID = "gbifid", datasetKey = "datasetkey", occurrenceID = "occurrenceid",
    scientificName = "scientificname", occurrenceStatus = "occurrencestatus",
    individualCount = "individualcount", publishingOrgKey = "publishingorgkey",
    decimalLatitude = "decimallatitude", decimalLongitude = "decimallongitude",
    coordinateUncertaintyInMeters = "coordinateuncertaintyinmeters",
    coordinatePrecision = "coordinateprecision", eventDate = "eventdate",
    basisOfRecord = "basisofrecord"
  )
  for (canonical in names(map)) {
    lower <- map[[canonical]]
    if (!canonical %in% names(df) && lower %in% names(df)) {
      names(df)[names(df) == lower] <- canonical
    }
  }
  df
}

# Clean + thin one species file. Returns a one-row status tibble (never stops).
clean_one_file <- function(input_path, output_dir, cfg = config) {
  clean_name <- basename(input_path) |>
    str_remove("_occs\\.parquet$") |>
    janitor::make_clean_names()

  res <- tryCatch({
    occ <- arrow::read_parquet(input_path) |> .normalise_gbif_cols()
    cl  <- cfg$cleaning

    d <- occ |>
      dplyr::filter(!basisOfRecord %in% cl$exclude_basis_of_record) |>
      dplyr::filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) |>
      cc_dupl() |> cc_zero() |> cc_equ() |> cc_val() |>
      cc_cap(buffer = 2000) |> cc_cen(buffer = 2000) |>
      cc_gbif(buffer = 2000) |> cc_inst(buffer = 500) |>
      dplyr::filter(coordinateUncertaintyInMeters < cl$max_coord_uncertainty_m |
                      is.na(coordinateUncertaintyInMeters)) |>
      dplyr::mutate(coordinatePrecision = as.numeric(coordinatePrecision)) |>
      dplyr::filter(coordinatePrecision > cl$min_coord_precision |
                      is.na(coordinatePrecision)) |>
      dplyr::mutate(individualCount = as.numeric(individualCount)) |>
      dplyr::filter(individualCount > 0 | is.na(individualCount)) |>
      dplyr::filter(individualCount < cl$max_individual_count | is.na(individualCount))

    if (nrow(d) == 0) {
      return(tibble::tibble(species = clean_name, rows = 0L,
                            file_path = NA_character_, status = "empty_after_clean"))
    }

    # H3 thinning: one record per hex cell.
    d_sf <- st_as_sf(d, coords = c("decimalLongitude", "decimalLatitude"),
                     crs = cfg$crs$working)
    d_thin <- d_sf |>
      dplyr::mutate(h3_index = h3jsr::point_to_cell(d_sf, res = cfg$hex_resolution)) |>
      dplyr::group_by(h3_index) |>
      dplyr::mutate(tot_obs = dplyr::n_distinct(gbifID)) |>
      dplyr::slice_sample(n = 1) |>
      dplyr::ungroup() |>
      dplyr::mutate(decimalLongitude = st_coordinates(geometry)[, "X"],
                    decimalLatitude  = st_coordinates(geometry)[, "Y"]) |>
      st_drop_geometry()

    out_path <- file.path(ensure_dir(output_dir), paste0(clean_name, "_clean.parquet"))
    arrow::write_parquet(d_thin, out_path)
    tibble::tibble(species = clean_name, rows = nrow(d_thin),
                   file_path = out_path, status = "written")
  }, error = function(e) {
    tibble::tibble(species = clean_name, rows = NA_integer_,
                   file_path = NA_character_,
                   status = paste0("error: ", conditionMessage(e)))
  })
  res
}

# Clean one taxon-year: every parsed occ file in that folder.
clean_taxon_year <- function(taxon, year, cfg = config) {
  in_dir <- taxon_year_occ_dir(taxon, year, cfg)
  files <- list.files(in_dir, pattern = "_occs\\.parquet$", full.names = TRUE)
  if (length(files) == 0) return(tibble::tibble())
  purrr::map_dfr(files, clean_one_file, output_dir = in_dir, cfg = cfg) |>
    dplyr::mutate(taxon = taxon, year = year, .before = 1)
}

run_clean_gbif <- function(taxa = config$taxa, years = config$years, cfg = config) {
  grid <- tidyr::expand_grid(taxon = taxa, year = years)
  results <- purrr::pmap_dfr(list(grid$taxon, grid$year),
                             function(taxon, year) clean_taxon_year(taxon, year, cfg))
  log_path <- file.path(cfg$paths$output, "indicators", "species_populations",
                        "logs", "clean_gbif.csv")
  ensure_parent_dir(log_path)
  readr::write_csv(results, log_path)
  results
}
