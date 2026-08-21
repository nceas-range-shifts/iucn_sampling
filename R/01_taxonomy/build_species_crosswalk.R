# =============================================================================
# R/01_taxonomy/build_species_crosswalk.R  —  STAGE 0 driver
# =============================================================================
#
# Turns each taxon's raw input species list into ONE authoritative crosswalk
# that every downstream stage joins on. Runs FIRST in the pipeline.
#
# What it adds beyond resolve_species_names():
#   - range_shapefile_name / range_shapefile_path for shapefile-range taxa.
#     Range shapefiles carry the vendor's own species name (the shapefile_name_
#     field, e.g. "SName"), which may not match the GBIF canonical name. We
#     resolve those shapefile names through the GBIF backbone and join them to
#     the crosswalk on the shared gbif_usage_key -- the same key-bridge the
#     resolver already uses for eBird. This is the join key that lets range
#     cells (keyed by shapefile name) line up with occurrence/richness data
#     (keyed by GBIF canonical name).
#
# Output: data/processed/species_crosswalk.csv (one row per species, all taxa).
#
# Requires config/config.R + config/taxa_config.R sourced first.
#
# Usage:
#   source("config/config.R"); source("config/taxa_config.R")
#   source("R/00_utils.R")
#   source("R/01_taxonomy/build_species_crosswalk.R")
#   crosswalk <- build_species_crosswalk()

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(glue)
  library(janitor)
  library(sf)
  library(here)
})

source(here::here("R", "01_taxonomy", "resolve_species_names.R"))

# Load one taxon's raw species list from its configured CSV or Google Sheet.
# When reading a shared CSV that carries a `taxon` column, filter to this taxon.
load_species_list <- function(taxon, tcfg) {
  src <- tcfg$species_source
  has_sheet <- !is.null(src$sheet_url) && !is.na(src$sheet_url) && nzchar(src$sheet_url)

  if (has_sheet) {
    googlesheets4::gs4_deauth()
    df <- googlesheets4::read_sheet(ss = src$sheet_url, sheet = src$sheet_name)
  } else if (!is.null(src$csv_path) && file.exists(src$csv_path)) {
    df <- readr::read_csv(src$csv_path, show_col_types = FALSE)
  } else {
    stop(glue("No species source found for taxon '{taxon}' (csv_path missing)."))
  }

  # Optional filter when a shared list holds multiple taxa.
  if (!is.null(src$taxon_filter)) {
    df_clean <- janitor::clean_names(df)
    if ("taxon" %in% names(df_clean)) {
      df <- df[tolower(df_clean$taxon) == tolower(src$taxon_filter), , drop = FALSE]
    }
  }
  df
}

# Read the species-name attribute out of every range shapefile in a directory.
# Only the attribute table is needed, so geometry is dropped immediately.
harvest_shapefile_names <- function(shapefile_dir, name_field) {
  if (is.null(shapefile_dir) || !dir.exists(shapefile_dir)) {
    return(tibble::tibble(range_shapefile_name = character(),
                          range_shapefile_path = character()))
  }
  shp <- list.files(shapefile_dir, pattern = "\\.shp$", full.names = TRUE)
  purrr::map_dfr(shp, function(f) {
    out <- tryCatch({
      atts <- sf::st_read(f, quiet = TRUE) |> sf::st_drop_geometry()
      if (!name_field %in% names(atts)) {
        warning(glue("shapefile '{basename(f)}' has no field '{name_field}'; skipped"))
        return(NULL)
      }
      nm <- unique(as.character(atts[[name_field]]))
      nm <- nm[!is.na(nm) & nm != ""]
      if (length(nm) == 0) return(NULL)
      tibble::tibble(range_shapefile_name = nm, range_shapefile_path = f)
    }, error = function(e) {
      warning(glue("could not read '{basename(f)}': {conditionMessage(e)}"))
      NULL
    })
    out
  })
}

# Resolve shapefile names to GBIF keys and attach them to the crosswalk on
# gbif_usage_key. Rows in the crosswalk with no matching shapefile get NA.
attach_shapefile_link <- function(crosswalk, shp_names) {
  if (nrow(shp_names) == 0) {
    return(dplyr::mutate(crosswalk,
                         range_shapefile_name = NA_character_,
                         range_shapefile_path = NA_character_))
  }

  keys <- rgbif::name_backbone_checklist(shp_names$range_shapefile_name) |>
    dplyr::transmute(
      range_shapefile_name = verbatim_name,
      gbif_usage_key = as.numeric(speciesKey)
    )

  link <- shp_names |>
    dplyr::left_join(keys, by = "range_shapefile_name") |>
    dplyr::filter(!is.na(gbif_usage_key)) |>
    dplyr::distinct(gbif_usage_key, .keep_all = TRUE)

  crosswalk |>
    dplyr::left_join(link, by = "gbif_usage_key")
}

# Build one taxon's crosswalk (names resolved + shapefile linkage where used).
build_taxon_crosswalk <- function(taxon, tcfg, cfg = config) {
  raw <- load_species_list(taxon, tcfg)

  cw <- resolve_species_names(
    spp_df = raw,
    sci_name_col = tcfg$species_source$sci_name_col,
    taxon_group = taxon           # stamp group so reptiles aren't left NA
  )

  if (identical(tcfg$range_source, "shapefile")) {
    shp_names <- harvest_shapefile_names(tcfg$shapefile_dir, tcfg$shapefile_name_field)
    cw <- attach_shapefile_link(cw, shp_names)

    # Flag species whose range shapefile we could not match -> review.
    cw <- cw |>
      dplyr::mutate(
        needs_review = needs_review | is.na(range_shapefile_name),
        review_reason = dplyr::if_else(
          is.na(range_shapefile_name) & !override,
          stringr::str_trim(paste(review_reason, "no range shapefile matched", sep = "; ")),
          review_reason
        )
      )
  } else {
    cw <- dplyr::mutate(cw,
                        range_shapefile_name = NA_character_,
                        range_shapefile_path = NA_character_)
  }

  dplyr::mutate(cw, taxon = taxon, .before = 1)
}

# Build the master crosswalk across all configured taxa and write it once.
build_species_crosswalk <- function(taxa = config$taxa,
                                    taxa_cfg = taxa_config,
                                    cfg = config,
                                    write = TRUE) {
  crosswalk <- purrr::map_dfr(taxa, function(taxon) {
    tcfg <- taxa_cfg[[taxon]]
    if (is.null(tcfg)) stop(glue("Taxon '{taxon}' not found in taxa_config."))
    build_taxon_crosswalk(taxon, tcfg, cfg)
  })

  if (isTRUE(write)) {
    out <- cfg$paths$crosswalk
    ensure_parent_dir(out)
    readr::write_csv(crosswalk, out)
    message(glue(
      "[crosswalk] {nrow(crosswalk)} species across {length(taxa)} taxa -> {out} ",
      "({sum(crosswalk$needs_review)} flagged for review)"
    ))
  }

  crosswalk
}
