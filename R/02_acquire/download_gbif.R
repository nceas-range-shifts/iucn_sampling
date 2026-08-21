# =============================================================================
# R/02_acquire/download_gbif.R  —  submit GBIF occurrence downloads
# =============================================================================
#
# Reads the master crosswalk built in Stage 0, takes each taxon's already-
# resolved GBIF usageKeys, and submits one filtered occurrence download per
# taxon. Name resolution does NOT happen here (that is Stage 0's job) — this
# script only downloads.
#
# API-gated: needs GBIF credentials in .Renviron (GBIF_USER/PWD/EMAIL).
#
# Requires config + taxa_config + R/00_utils.R sourced first.
#
# Usage:
#   res <- run_gbif_downloads()                 # all configured taxa
#   res <- run_gbif_downloads(taxa = "reptile", year_range = c(2014, 2015))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(glue)
  library(rgbif)
  library(here)
})

# Load the master crosswalk (built by Stage 0), optionally filtered to a taxon.
load_crosswalk <- function(taxon = NULL, cfg = config) {
  if (!file.exists(cfg$paths$crosswalk)) {
    stop(glue("Crosswalk not found at {cfg$paths$crosswalk}. Run Stage 0 ",
              "(build_species_crosswalk()) first."))
  }
  cw <- readr::read_csv(cfg$paths$crosswalk, show_col_types = FALSE)
  if (!is.null(taxon)) cw <- dplyr::filter(cw, taxon == !!taxon)
  cw
}

# GBIF usageKeys to download. Flagged (needs_review) rows are included by
# default; set include_flagged = FALSE to drop them.
get_gbif_keys <- function(crosswalk, include_flagged = TRUE) {
  keys <- dplyr::filter(crosswalk, !is.na(gbif_usage_key))
  if (!isTRUE(include_flagged)) keys <- dplyr::filter(keys, !needs_review)
  keys |> dplyr::distinct(gbif_usage_key) |> dplyr::pull(gbif_usage_key)
}

build_year_predicates <- function(year_range) {
  if (length(year_range) != 2) stop("year_range must be c(start_year, end_year).")
  start_year <- as.integer(year_range[[1]]); end_year <- as.integer(year_range[[2]])
  if (start_year > end_year) stop("year_range start must be <= end.")
  list(rgbif::pred_gte("year", start_year), rgbif::pred_lte("year", end_year))
}

submit_gbif_download <- function(gbif_keys, cfg = config) {
  f <- cfg$gbif_filters
  preds <- c(
    list(
      rgbif::pred("hasGeospatialIssue", f$has_geo_issue),
      rgbif::pred("hasCoordinate", f$has_coordinate),
      rgbif::pred_in("basisOfRecord", f$basis_of_record),
      rgbif::pred_in("taxonKey", gbif_keys),
      rgbif::pred("country", cfg$region$country_code),
      rgbif::pred("stateProvince", cfg$region$admin_name)
    ),
    build_year_predicates(range(cfg$years)),
    list(format = f$download_format)
  )
  do.call(rgbif::occ_download, preds)
}

# Download one taxon; returns a one-row result tibble (a run-log entry).
download_one_taxon <- function(taxon, cfg = config, taxa_cfg = taxa_config,
                               wait = TRUE, include_flagged = TRUE) {
  tcfg <- taxa_cfg[[taxon]]
  if (is.null(tcfg)) stop(glue("Taxon '{taxon}' not in taxa_config."))

  crosswalk <- load_crosswalk(taxon, cfg)
  gbif_keys <- get_gbif_keys(crosswalk, include_flagged = include_flagged)
  if (length(gbif_keys) == 0) stop(glue("No GBIF keys resolved for '{taxon}'."))

  dl <- submit_gbif_download(gbif_keys, cfg)

  out_dir <- dirname(tcfg$gbif_dataset_path)
  if (isTRUE(wait)) {
    rgbif::occ_download_wait(dl)
    ensure_dir(out_dir)
    rgbif::occ_download_get(dl, path = out_dir)
  }

  tibble::tibble(
    taxon = taxon, status = "ok",
    n_species = nrow(crosswalk),
    n_flagged = sum(crosswalk$needs_review, na.rm = TRUE),
    n_gbif_keys = length(gbif_keys),
    download_key = dl$key, output_dir = out_dir,
    waited = wait, run_time_utc = as.character(Sys.time())
  )
}

# Driver over taxa, with per-taxon error capture so one failure doesn't abort.
run_gbif_downloads <- function(taxa = config$taxa, cfg = config,
                               taxa_cfg = taxa_config, wait = TRUE,
                               include_flagged = TRUE) {
  purrr::map_dfr(taxa, function(taxon) {
    tryCatch(
      download_one_taxon(taxon, cfg, taxa_cfg, wait, include_flagged),
      error = function(e) tibble::tibble(
        taxon = taxon, status = "error",
        n_species = NA_integer_, n_flagged = NA_integer_, n_gbif_keys = NA_integer_,
        download_key = NA_character_, output_dir = NA_character_,
        waited = wait, run_time_utc = as.character(Sys.time()),
        error_message = conditionMessage(e)
      )
    )
  })
}
