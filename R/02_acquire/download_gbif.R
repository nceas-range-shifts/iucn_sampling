### =============================================================================
### R/02_acquire/download_gbif.R  —  submit GBIF occurrence downloads
### =============================================================================
###
### Reads the master crosswalk built in Stage 0, takes each taxon's already-
### resolved GBIF usageKeys, and submits one filtered occurrence download per
### taxon. Name resolution does NOT happen here (that is Stage 0's job) — this
### script only downloads.
###
### API-gated: needs GBIF credentials in config$credentials (default: read from
### .Renviron GBIF_USER/GBIF_PWD/GBIF_EMAIL).
###
### Requires config + taxa_config + R/00_utils.R sourced first.
###
### Usage:
###   res <- run_gbif_downloads()                 # all configured taxa
###   res <- run_gbif_downloads(taxa = "reptile", year_range = c(2014, 2015))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(glue)
  library(rgbif)
  library(here)
})

### Load the master crosswalk (built by Stage 0), optionally filtered to a taxon.
load_crosswalk <- function(taxon = NULL, cfg = config) {
  if (!file.exists(cfg$paths$crosswalk)) {
    stop(glue("Crosswalk not found at {cfg$paths$crosswalk}. Run Stage 0 ",
              "(build_species_crosswalk()) first."))
  }
  cw <- readr::read_csv(cfg$paths$crosswalk, show_col_types = FALSE)
  if (!is.null(taxon)) cw <- dplyr::filter(cw, taxon == !!taxon)
  cw
}

### GBIF usageKeys to download. Flagged (needs_review) rows are included by
### default; set include_flagged = FALSE to drop them.
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

### Bounding box (xmin, xmax, ymin, ymax) for the configured region, or NULL for
### no spatial restriction. Unnamed config bboxes are read as xmin, xmax, ymin,
### ymax (same convention as get_region_boundary() in build_hex_grid.R).
region_bbox <- function(cfg = config) {
  src <- cfg$region$boundary_source
  if (identical(src, "global")) return(NULL)

  if (identical(src, "bbox")) {
    bb <- cfg$region$bbox
    nms <- c("xmin", "xmax", "ymin", "ymax")
    if (length(bb) != 4) stop("region$bbox must have 4 values: xmin, xmax, ymin, ymax.")
    if (is.null(names(bb))) names(bb) <- nms
    if (!all(nms %in% names(bb))) stop("region$bbox names must be xmin, xmax, ymin, ymax.")
    bb <- bb[nms]
  } else {
    ### Polygon sources: restrict to the boundary's extent.
    ext <- sf::st_bbox(sf::st_transform(get_region_boundary(cfg), 4326))
    bb <- c(xmin = ext[["xmin"]], xmax = ext[["xmax"]],
            ymin = ext[["ymin"]], ymax = ext[["ymax"]])
  }

  ### A bbox spanning the whole globe is no restriction.
  if (bb[["xmin"]] <= -180 && bb[["xmax"]] >= 180 &&
      bb[["ymin"]] <= -90  && bb[["ymax"]] >= 90) return(NULL)
  bb
}

### GBIF geometry predicate for the region, or NULL when none is needed.
### GBIF requires counter-clockwise polygon rings.
build_region_predicate <- function(cfg = config) {
  bb <- region_bbox(cfg)
  if (is.null(bb)) return(NULL)
  wkt <- glue::glue(
    "POLYGON(({bb[['xmin']]} {bb[['ymin']]}, {bb[['xmax']]} {bb[['ymin']]}, ",
    "{bb[['xmax']]} {bb[['ymax']]}, {bb[['xmin']]} {bb[['ymax']]}, ",
    "{bb[['xmin']]} {bb[['ymin']]}))"
  )
  rgbif::pred_within(as.character(wkt))
}

submit_gbif_download <- function(gbif_keys, cfg = config) {
  f <- cfg$gbif_filters
  preds <- c(
    list(
      rgbif::pred("hasGeospatialIssue", f$has_geo_issue),
      rgbif::pred("hasCoordinate", f$has_coordinate),
      rgbif::pred_in("basisOfRecord", f$basis_of_record),
      rgbif::pred_in("taxonKey", gbif_keys)
    ),
    ### Region filters are optional. admin_name is only a label, so it is not
    ### sent to GBIF. Country is added only if set; the geometry filter is
    ### omitted for a global region.
    if (!is.null(cfg$region$country_code))
      list(rgbif::pred("country", cfg$region$country_code)),
    if (!is.null(region_pred <- build_region_predicate(cfg)))
      list(region_pred),
    build_year_predicates(range(cfg$years)),
    list(format = f$download_format),
    ### Credentials from config; empty values are dropped so rgbif falls back to
    ### its own GBIF_USER/GBIF_PWD/GBIF_EMAIL lookup.
    Filter(function(x) !is.null(x) && nzchar(x),
           list(user = cfg$credentials$gbif_user,
                pwd = cfg$credentials$gbif_pwd,
                email = cfg$credentials$gbif_email))
  )
  do.call(rgbif::occ_download, preds)
}

### --- Download manifest -------------------------------------------------------
### One JSON file per taxon output dir records what was downloaded and from which
### request, so unchanged requests are not resubmitted.
manifest_path <- function(out_dir) file.path(out_dir, "gbif_download_manifest.json")

### Hash of everything that determines the download's content. Changes to the
### species keys, filters, years, or region produce a new hash.
gbif_request_hash <- function(gbif_keys, cfg = config) {
  rlang::hash(list(
    taxon_keys = sort(unique(as.numeric(gbif_keys))),
    ### `enabled` is excluded so toggling the source doesn't change the request.
    filters    = cfg$gbif_filters[setdiff(names(cfg$gbif_filters), "enabled")],
    years      = range(cfg$years),
    country    = cfg$region$country_code,
    region     = region_bbox(cfg)
  ))
}

read_gbif_manifest <- function(out_dir) {
  p <- manifest_path(out_dir)
  if (!file.exists(p)) return(NULL)
  tryCatch(jsonlite::read_json(p, simplifyVector = TRUE), error = function(e) NULL)
}

### Returns NULL if the existing download is reusable, else a short reason string.
gbif_download_stale_reason <- function(manifest, request_hash, out_dir,
                                       max_age_days = Inf) {
  if (is.null(manifest)) return("no manifest")
  if (!identical(manifest$request_hash, request_hash)) return("request changed")
  zip <- file.path(out_dir, manifest$file)
  if (is.null(manifest$file) || !file.exists(zip)) return("downloaded file missing")
  age <- as.numeric(difftime(Sys.time(), as.POSIXct(manifest$downloaded_utc, tz = "UTC"),
                             units = "days"))
  if (is.na(age) || age > max_age_days) {
    return(glue("older than max_age_days ({round(age, 1)} d)"))
  }
  NULL
}

### Download one taxon; returns a one-row result tibble (a run-log entry).
### Skips submission when the manifest matches the current request and is fresh
### (see `max_age_days`); `force = TRUE` always re-downloads.
download_one_taxon <- function(taxon, cfg = config, taxa_cfg = taxa_config,
                               wait = TRUE, include_flagged = TRUE,
                               max_age_days = Inf, force = FALSE) {
  tcfg <- taxa_cfg[[taxon]]
  if (is.null(tcfg)) stop(glue("Taxon '{taxon}' not in taxa_config."))

  crosswalk <- load_crosswalk(taxon, cfg)
  gbif_keys <- get_gbif_keys(crosswalk, include_flagged = include_flagged)
  if (length(gbif_keys) == 0) stop(glue("No GBIF keys resolved for '{taxon}'."))

  out_dir <- dirname(tcfg$gbif_dataset_path)
  req_hash <- gbif_request_hash(gbif_keys, cfg)

  summary_row <- function(status, key) {
    tibble::tibble(
      taxon = taxon, status = status,
      n_species = nrow(crosswalk),
      n_flagged = sum(crosswalk$needs_review, na.rm = TRUE),
      n_gbif_keys = length(gbif_keys),
      download_key = key, output_dir = out_dir,
      waited = wait, run_time_utc = as.character(Sys.time())
    )
  }

  if (!isTRUE(force)) {
    manifest <- read_gbif_manifest(out_dir)
    reason <- gbif_download_stale_reason(manifest, req_hash, out_dir, max_age_days)
    if (is.null(reason)) {
      message(glue("[gbif] {taxon}: using CACHED download {manifest$download_key} ",
                   "(downloaded {manifest$downloaded_utc} UTC; no new download submitted)."))
      return(summary_row("skipped", manifest$download_key))
    }
    message(glue("[gbif] {taxon}: downloading ({reason})."))
  }

  dl <- submit_gbif_download(gbif_keys, cfg)
  key <- as.character(dl)

  if (isTRUE(wait)) {
    rgbif::occ_download_wait(dl)
    ensure_dir(out_dir)
    got <- rgbif::occ_download_get(dl, path = out_dir, overwrite = TRUE)
    jsonlite::write_json(
      list(taxon = taxon, download_key = key, request_hash = req_hash,
           file = basename(as.character(got)),
           n_gbif_keys = length(gbif_keys),
           downloaded_utc = format(Sys.time(), tz = "UTC", "%Y-%m-%d %H:%M:%S")),
      manifest_path(out_dir), auto_unbox = TRUE, pretty = TRUE
    )
  }

  summary_row("ok", key)
}

### Driver over taxa, with per-taxon error capture so one failure doesn't abort.
### max_age_days: re-download if the stored download is older than this
### (default: config$gbif_max_age_days, else 180; Inf = never).
run_gbif_downloads <- function(taxa = config$taxa, cfg = config,
                               taxa_cfg = taxa_config, wait = TRUE,
                               include_flagged = TRUE,
                               max_age_days = cfg$gbif_max_age_days %||% 180,
                               force = FALSE) {
  if (!source_enabled("gbif", cfg)) {
    message("[gbif] disabled (config$gbif_filters$enabled = FALSE); skipping.")
    return(tibble::tibble())
  }
  purrr::map_dfr(taxa, function(taxon) {
    tryCatch(
      download_one_taxon(taxon, cfg, taxa_cfg, wait, include_flagged,
                         max_age_days, force),
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
