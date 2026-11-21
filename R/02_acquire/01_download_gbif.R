### =============================================================================
### R/02_acquire/01_download_gbif.R  —  submit GBIF occurrence downloads
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
### Two ways to get the data, chosen by config$gbif_source:
###  'download': downloads are submitted in a batch and polled together (up to 3 in
###    flight, the GBIF per-user limit); in-flight downloads are recorded in
###    gbif_pending.json so a rerun re-attaches instead of resubmitting.
###  'snapshot': one DuckDB pass over the monthly GBIF snapshot on AWS S3 (date
###    pinned in config$gbif_snapshot_date), split into per-taxon datasets. No
###    credentials or queue; matches on specieskey. See run_gbif_snapshot().
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
  library(DBI)
})

### Load the master crosswalk (built by Stage 0), optionally filtered to a taxon.
load_crosswalk <- function(taxon = NULL, cfg = config) {
  if (!file.exists(cfg$paths$crosswalk)) {
    stop(glue("Crosswalk not found at {cfg$paths$crosswalk}. Run Stage 0 ",
              "(build_species_crosswalk()) first."))
  }
  ### gbif_usage_key is read as character: keys are alphanumeric (e.g. "5384V"), and
  ### an all-digit subset would otherwise be guessed as numeric.
  cw <- readr::read_csv(cfg$paths$crosswalk, show_col_types = FALSE,
                        col_types = readr::cols(gbif_usage_key = readr::col_character(),
                                                .default = readr::col_guess()))
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
### ymax (same convention as get_region_boundary() in 03_build_hex_grid.R).
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

### TRUE if the extracted dataset exists (a non-empty partitioned directory, or a
### single parquet file).
gbif_dataset_present <- function(dataset_path) {
  if (dir.exists(dataset_path)) return(length(list.files(dataset_path)) > 0)
  file.exists(dataset_path)
}

### Unzip a GBIF download into out_dir (producing dataset_path), replacing any
### previous copy, then delete the zip. The zip is kept if extraction fails.
extract_gbif_zip <- function(zip_path, out_dir, dataset_path) {
  zip_path <- as.character(zip_path)
  ### Extract to a staging dir first so a failed unzip doesn't destroy the old data.
  stage <- file.path(out_dir, ".extract_tmp")
  unlink(stage, recursive = TRUE)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  utils::unzip(zip_path, exdir = stage)

  staged <- file.path(stage, basename(dataset_path))
  if (!gbif_dataset_present(staged)) {
    stop(glue("{basename(dataset_path)} not found inside {basename(zip_path)}; ",
              "zip kept for inspection."))
  }
  unlink(dataset_path, recursive = TRUE)
  ### Move everything (dataset plus any citation/metadata files) into out_dir.
  for (f in list.files(stage, full.names = TRUE, all.files = TRUE, no.. = TRUE)) {
    file.rename(f, file.path(out_dir, basename(f)))
  }
  file.remove(zip_path)
  invisible(dataset_path)
}

### Returns NULL if the existing download is reusable, else a short reason string.
### The zip is deleted after extraction, so reusability is judged by whether the
### extracted dataset (occurrence.parquet) is present and non-empty.
gbif_download_stale_reason <- function(manifest, request_hash, dataset_path,
                                       max_age_days = Inf) {
  if (is.null(manifest)) return("no manifest")
  if (!identical(manifest$request_hash, request_hash)) return("request changed")
  if (!gbif_dataset_present(dataset_path)) return("occurrence dataset missing")
  age <- as.numeric(difftime(Sys.time(), as.POSIXct(manifest$downloaded_utc, tz = "UTC"),
                             units = "days"))
  if (is.na(age) || age > max_age_days) {
    return(glue("older than max_age_days ({round(age, 1)} d)"))
  }
  NULL
}

### --- Batch submit + concurrent wait -------------------------------------------
### GBIF builds each download server-side, so waiting on taxa one at a time
### makes total time the SUM of every queue wait. Instead we keep up to
### `max_concurrent` downloads in flight and poll them together.
### GBIF allows 3 incomplete (preparing/running) downloads per user, dropping to 1
### when the platform is busy; extra submissions are refused with HTTP 420/429.
### (rgbif::occ_download_queue() is not used: it hard-codes 3, cannot resume,
### and never stops polling a FAILED download.)

gbif_log <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z", tz = "America/Los_Angeles"),
          " [gbif] ", ...)
}

### Written at submit time, removed once the download is extracted. Lets a rerun
### re-attach to an in-flight download instead of resubmitting after a crash.
pending_path <- function(out_dir) file.path(out_dir, "gbif_pending.json")

gbif_failed_states <- c("killed", "cancelled", "failed", "file_erased")

### Resolve one taxon into a job record. state is one of:
###   skipped (cached), active (re-attached to an in-flight key), queued (to submit).
plan_gbif_taxon <- function(taxon, cfg, taxa_cfg, include_flagged, max_age_days, force) {
  tcfg <- taxa_cfg[[taxon]]
  if (is.null(tcfg)) stop(glue("Taxon '{taxon}' not in taxa_config."))

  crosswalk <- load_crosswalk(taxon, cfg)
  gbif_keys <- get_gbif_keys(crosswalk, include_flagged = include_flagged)
  if (length(gbif_keys) == 0) stop(glue("No GBIF keys resolved for '{taxon}'."))

  dataset_path <- tcfg$gbif_dataset_path
  out_dir <- dirname(dataset_path)
  req_hash <- gbif_request_hash(gbif_keys, cfg)

  job <- list(
    taxon = taxon, state = "queued", key = NA_character_, error = NA_character_,
    gbif_keys = gbif_keys, dataset_path = dataset_path, out_dir = out_dir,
    req_hash = req_hash, n_species = nrow(crosswalk),
    n_flagged = sum(crosswalk$needs_review, na.rm = TRUE)
  )

  if (!isTRUE(force)) {
    manifest <- read_gbif_manifest(out_dir)
    reason <- gbif_download_stale_reason(manifest, req_hash, dataset_path, max_age_days)
    if (is.null(reason)) {
      gbif_log(glue("{taxon}: using CACHED download {manifest$download_key} ",
                    "(downloaded {manifest$downloaded_utc} UTC; nothing submitted)."))
      job$state <- "skipped"
      job$key <- manifest$download_key
      return(job)
    }
    gbif_log(glue("{taxon}: needs download ({reason})."))
  }

  ### Re-attach to an in-flight identical request, even when force = TRUE (it is
  ### the same request, already paid for in queue time).
  pend <- tryCatch(jsonlite::read_json(pending_path(out_dir), simplifyVector = TRUE),
                   error = function(e) NULL)
  if (!is.null(pend) && identical(pend$request_hash, req_hash)) {
    st <- tolower(rgbif::occ_download_meta(pend$download_key)$status)
    if (!st %in% gbif_failed_states) {
      gbif_log(glue("{taxon}: re-attaching to pending download {pend$download_key} ({st})."))
      job$state <- "active"
      job$key <- pend$download_key
      return(job)
    }
    gbif_log(glue("{taxon}: pending download {pend$download_key} is {st}; resubmitting."))
    unlink(pending_path(out_dir))
  }
  job
}

### Get + extract a SUCCEEDED download and write its manifest.
finalize_gbif_job <- function(job) {
  ensure_dir(job$out_dir)
  got <- rgbif::occ_download_get(job$key, path = job$out_dir, overwrite = TRUE)
  extract_gbif_zip(got, job$out_dir, job$dataset_path)
  jsonlite::write_json(
    list(taxon = job$taxon, download_key = job$key, request_hash = job$req_hash,
         dataset = basename(job$dataset_path),
         n_gbif_keys = length(job$gbif_keys),
         downloaded_utc = format(Sys.time(), tz = "UTC", "%Y-%m-%d %H:%M:%S")),
    manifest_path(job$out_dir), auto_unbox = TRUE, pretty = TRUE
  )
  unlink(pending_path(job$out_dir))
  invisible(job)
}

### Download one taxon; returns a one-row result tibble. Thin wrapper kept for
### callers that want a single taxon; all logic lives in run_gbif_downloads().
download_one_taxon <- function(taxon, cfg = config, taxa_cfg = taxa_config,
                               wait = TRUE, include_flagged = TRUE,
                               max_age_days = Inf, force = config$force) {
  run_gbif_downloads(taxa = taxon, cfg = cfg, taxa_cfg = taxa_cfg, wait = wait,
                     include_flagged = include_flagged, max_age_days = max_age_days,
                     force = force)
}

### Driver over taxa, with per-taxon error capture so one failure doesn't abort.
### max_age_days: re-download if the stored download is older than this
### (default: config$gbif_max_age_days, else 180; Inf = never).
### max_concurrent: downloads in flight at once (GBIF cap is 3; auto-lowered if
###   GBIF refuses a submission, restored when a download finishes).
### poll_secs: seconds between status checks. max_wait_mins: stop polling after
###   this long; unfinished downloads stay pending and are re-attached next run.
### wait = FALSE: submit what fits under the cap, write pending files, return.
run_gbif_downloads <- function(taxa = config$taxa, cfg = config,
                               taxa_cfg = taxa_config, wait = TRUE,
                               include_flagged = TRUE,
                               max_age_days = cfg$gbif_max_age_days %||% 180,
                               force = config$force,
                               max_concurrent = cfg$gbif_max_concurrent %||% 3,
                               poll_secs = cfg$gbif_poll_secs %||% 30,
                               max_wait_mins = cfg$gbif_max_wait_mins %||% 360,
                               snapshot_parts = NULL) {
  if (!source_enabled("gbif", cfg)) {
    message("[gbif] disabled (config$gbif_filters$enabled = FALSE); skipping.")
    return(tibble::tibble())
  }
  ### for new alphanumeric speciesKey taxonomic keys from Catalog of Life
  check_package_version('rgbif', '3.9')

  ### config$gbif_source = 'snapshot': one pass over the monthly AWS snapshot
  ### instead of the download API (see run_gbif_snapshot()). snapshot_parts limits
  ### the scan to the first N parts, for testing only.
  if (identical(cfg$gbif_source %||% "download", "snapshot")) {
    return(run_gbif_snapshot(taxa, cfg, taxa_cfg, include_flagged, max_age_days, force,
                             part_limit = snapshot_parts))
  }

  ### 1. Plan: resolve every taxon (cached / re-attach / to submit).
  jobs <- stats::setNames(lapply(taxa, function(taxon) {
    tryCatch(
      plan_gbif_taxon(taxon, cfg, taxa_cfg, include_flagged, max_age_days, force),
      error = function(e) list(
        taxon = taxon, state = "error", key = NA_character_,
        gbif_keys = NULL, out_dir = NA_character_,
        n_species = NA_integer_, n_flagged = NA_integer_,
        error = conditionMessage(e)
      )
    )
  }), taxa)

  ids <- function(s) names(Filter(function(j) j$state == s, jobs))
  rate_limit_re <- "\\b(420|429)\\b|too many|rate.?limit|concurrent|simultaneous"
  cap <- max_concurrent
  t_start <- Sys.time()

  ### 2./3. Submit up to `cap`, poll together, finalize each as it succeeds.
  repeat {
    queued <- ids("queued")
    active <- ids("active")

    while (length(queued) > 0 && length(active) < cap) {
      tx <- queued[[1]]
      res <- tryCatch(submit_gbif_download(jobs[[tx]]$gbif_keys, cfg),
                      error = function(e) e)
      if (inherits(res, "error")) {
        if (grepl(rate_limit_re, conditionMessage(res), ignore.case = TRUE)) {
          cap <- max(1, length(active))
          gbif_log(glue("submission refused (rate limit); cap now {cap}: ",
                        "{conditionMessage(res)}"))
          break
        }
        jobs[[tx]]$state <- "error"
        jobs[[tx]]$error <- conditionMessage(res)
        gbif_log(glue("{tx}: submit FAILED: {conditionMessage(res)}"))
      } else {
        jobs[[tx]]$key <- as.character(res)
        jobs[[tx]]$state <- "active"
        ensure_dir(jobs[[tx]]$out_dir)
        jsonlite::write_json(
          list(taxon = tx, download_key = jobs[[tx]]$key,
               request_hash = jobs[[tx]]$req_hash,
               submitted_utc = format(Sys.time(), tz = "UTC", "%Y-%m-%d %H:%M:%S")),
          pending_path(jobs[[tx]]$out_dir), auto_unbox = TRUE, pretty = TRUE
        )
        gbif_log(glue("{tx}: submitted {jobs[[tx]]$key}"))
      }
      queued <- ids("queued")
      active <- ids("active")
    }

    if (!isTRUE(wait)) break
    if (length(active) == 0 && length(queued) == 0) break
    if (as.numeric(difftime(Sys.time(), t_start, units = "mins")) > max_wait_mins) {
      gbif_log(glue("max_wait_mins ({max_wait_mins}) reached; leaving ",
                    "{length(active)} pending and {length(queued)} unsubmitted."))
      break
    }

    ### Poll once; a failed status call is treated as transient and retried.
    for (tx in active) {
      st <- tryCatch(tolower(rgbif::occ_download_meta(jobs[[tx]]$key)$status),
                     error = function(e) NA_character_)
      if (is.na(st)) next
      if (st == "succeeded") {
        fin <- tryCatch(finalize_gbif_job(jobs[[tx]]), error = function(e) e)
        if (inherits(fin, "error")) {
          jobs[[tx]]$state <- "error"
          jobs[[tx]]$error <- conditionMessage(fin)
          gbif_log(glue("{tx}: download ok but get/extract FAILED: {conditionMessage(fin)}"))
        } else {
          jobs[[tx]]$state <- "done"
          gbif_log(glue("{tx}: extracted {jobs[[tx]]$key}"))
        }
        cap <- max_concurrent
      } else if (st %in% gbif_failed_states) {
        jobs[[tx]]$state <- "error"
        jobs[[tx]]$error <- glue("GBIF download {jobs[[tx]]$key} ended as '{st}'")
        unlink(pending_path(jobs[[tx]]$out_dir))
        gbif_log(glue("{tx}: {jobs[[tx]]$error}"))
        cap <- max_concurrent
      }
    }

    gbif_log(glue("progress: {length(ids('done'))} done, {length(ids('active'))} in flight, ",
                  "{length(ids('queued'))} queued, {length(ids('skipped'))} cached, ",
                  "{length(ids('error'))} failed ",
                  "({round(difftime(Sys.time(), t_start, units = 'mins'), 1)} min elapsed)"))
    if (length(ids("active")) > 0 || length(ids("queued")) > 0) Sys.sleep(poll_secs)
  }

  ### 4. One summary row per taxon.
  purrr::map_dfr(taxa, function(taxon) {
    j <- jobs[[taxon]]
    tibble::tibble(
      taxon = taxon,
      status = switch(j$state, skipped = "skipped", done = "ok", error = "error",
                      active = "pending", queued = "not_submitted"),
      n_species = j$n_species, n_flagged = j$n_flagged,
      n_gbif_keys = length(j$gbif_keys) %||% NA_integer_,
      download_key = j$key, output_dir = j$out_dir,
      waited = wait, run_time_utc = as.character(Sys.time()),
      error_message = j$error
    )
  })
}

### --- Snapshot path: one pass over the monthly GBIF snapshot on AWS ------------
### Alternative to the download API (config$gbif_source = 'snapshot'). GBIF
### publishes a full occurrence snapshot each month as Parquet in a public S3
### bucket (anonymous read). We scan it ONCE with DuckDB for the union of every
### queued taxon's names, then split the hits into each taxon's gbif_dataset_path
### (a directory holding part-0.parquet, which arrow::open_dataset() reads like a
### download). No queue, no credentials; the cost is one long scan (~1 h from a
### home connection for the full snapshot), the same whether 1 or 10 taxa.
###
### Differences from the download API, all by design:
###  - Matching is on the snapshot `specieskey` (alphanumeric, Catalog of Life
###    based) against the crosswalk's gbif_usage_key. There is NO fallback to
###    species names: crosswalk rows without a key are logged and not extracted.
###    The parse stage still matches on name_gbif, so downstream is unchanged.
###  - Data are as of the snapshot date (config$gbif_snapshot_date, pinned for
###    reproducibility; 'latest' resolves to the newest and is recorded).
###  - Geospatial-issue exclusion is rebuilt from the `issue` array (see
###    gbif_snapshot_geo_issues; UNVERIFIED against GBIF's own hasGeospatialIssue).

gbif_snapshot_geo_issues <- c("ZERO_COORDINATE", "COORDINATE_OUT_OF_RANGE",
                              "COORDINATE_INVALID", "COUNTRY_COORDINATE_MISMATCH")

gbif_snapshot_bucket <- function(cfg) cfg$gbif_snapshot_bucket %||% "gbif-open-data-us-east-1"

### Snapshot dates available in the bucket (anonymous S3 listing), ascending.
gbif_snapshot_dates <- function(cfg) {
  resp <- httr2::request(glue("https://{gbif_snapshot_bucket(cfg)}.s3.amazonaws.com/")) |>
    httr2::req_url_query(`list-type` = 2, prefix = "occurrence/", delimiter = "/") |>
    httr2::req_retry(max_tries = 3) |>
    httr2::req_perform() |>
    httr2::resp_body_string()
  hits <- regmatches(resp, gregexpr("occurrence/[0-9]{4}-[0-9]{2}-[0-9]{2}/", resp))[[1]]
  sort(unique(gsub("occurrence/|/", "", hits)))
}

### Pinned date from config, or the newest when set to 'latest'.
gbif_snapshot_resolve <- function(cfg) {
  want <- cfg$gbif_snapshot_date %||% "latest"
  avail <- gbif_snapshot_dates(cfg)
  if (length(avail) == 0) stop("No GBIF snapshots found in the bucket listing.")
  if (identical(want, "latest")) {
    want <- max(avail)
    gbif_log(glue("gbif_snapshot_date = 'latest' resolved to {want}; pin it in config ",
                  "for reproducibility."))
  }
  if (!want %in% avail) {
    stop(glue("Snapshot {want} not found. Available: {paste(tail(avail, 6), collapse = ', ')}"))
  }
  want
}

### citation.txt shipped with the snapshot (NA if it cannot be fetched).
gbif_snapshot_citation <- function(date, cfg) {
  tryCatch(
    httr2::request(glue("https://{gbif_snapshot_bucket(cfg)}.s3.amazonaws.com/",
                        "occurrence/{date}/citation.txt")) |>
      httr2::req_perform() |>
      httr2::resp_body_string() |>
      trimws(),
    error = function(e) NA_character_
  )
}

### Hash of everything that determines the extracted content.
gbif_snapshot_hash <- function(sp_keys, date, cols, cfg) {
  rlang::hash(list(
    source  = "snapshot", date = date, match = "specieskey",
    keys    = sort(unique(sp_keys)),
    columns = sort(cols),
    ### `enabled` and `download_format` do not affect snapshot content.
    filters = cfg$gbif_filters[setdiff(names(cfg$gbif_filters), c("enabled", "download_format"))],
    years   = range(cfg$years),
    country = cfg$region$country_code,
    region  = region_bbox(cfg)
  ))
}

### Resolve one taxon into a job record (skipped = cached, queued = needs scan).
plan_snapshot_taxon <- function(taxon, cfg, taxa_cfg, include_flagged, max_age_days,
                                force, date, cols) {
  tcfg <- taxa_cfg[[taxon]]
  if (is.null(tcfg)) stop(glue("Taxon '{taxon}' not in taxa_config."))

  cw <- load_crosswalk(taxon, cfg)
  if (!isTRUE(include_flagged)) cw <- dplyr::filter(cw, !needs_review)
  ### Match on specieskey only; rows without a key are reported, never name-matched.
  no_key <- is.na(cw$gbif_usage_key) | cw$gbif_usage_key == ""
  if (any(no_key)) {
    gbif_log(glue("{taxon}: {sum(no_key)} crosswalk row(s) have no gbif_usage_key and will ",
                  "NOT be extracted: ",
                  "{paste(utils::head(dplyr::coalesce(cw$name_gbif, cw$name_input)[no_key], 10), collapse = '; ')}"))
  }
  sp_keys <- sort(unique(cw$gbif_usage_key[!no_key]))
  if (length(sp_keys) == 0) stop(glue("No gbif_usage_key values in crosswalk for '{taxon}'."))
  sp_names <- sort(unique(stats::na.omit(dplyr::coalesce(cw$name_gbif, cw$name_accepted)[!no_key])))

  dataset_path <- tcfg$gbif_dataset_path
  out_dir <- dirname(dataset_path)
  req_hash <- gbif_snapshot_hash(sp_keys, date, cols, cfg)

  job <- list(
    taxon = taxon, state = "queued", error = NA_character_, n_records = NA_integer_,
    sp_keys = sp_keys, sp_names = sp_names,
    key_names = stats::setNames(
      dplyr::coalesce(cw$name_gbif, cw$name_accepted, cw$name_input)[!no_key],
      cw$gbif_usage_key[!no_key]),
    dataset_path = dataset_path, out_dir = out_dir,
    req_hash = req_hash, n_species = nrow(cw),
    n_flagged = sum(cw$needs_review, na.rm = TRUE)
  )

  if (!isTRUE(force)) {
    manifest <- read_gbif_manifest(out_dir)
    reason <- gbif_download_stale_reason(manifest, req_hash, dataset_path, max_age_days)
    if (is.null(reason)) {
      gbif_log(glue("{taxon}: using CACHED snapshot extract ",
                    "({manifest$snapshot_date %||% 'unknown date'}; nothing scanned)."))
      job$state <- "skipped"
      job$n_records <- manifest$n_records %||% NA_integer_
      return(job)
    }
    gbif_log(glue("{taxon}: needs snapshot extract ({reason})."))
  }
  job
}

### Driver for the snapshot path; same per-taxon summary-row contract as
### run_gbif_downloads(). part_limit scans only the first N parts (for testing;
### point taxa_cfg/cfg$paths$data_raw at a scratch dir so real outputs are untouched).
### Chunk checkpoints live in <data_raw>/gbif/_snapshot_work/<date>_<hash>/ and are
### reused on rerun when `resume` (they are keyed by the exact request, so they
### cannot be stale); deleted after a fully successful run.
run_gbif_snapshot <- function(taxa, cfg, taxa_cfg, include_flagged, max_age_days, force,
                              part_limit = NULL,
                              n_chunks = cfg$gbif_snapshot_chunks %||% 20,
                              threads = cfg$gbif_snapshot_threads %||% 6,
                              memory_limit = cfg$gbif_snapshot_memory %||% "12GB",
                              resume = cfg$gbif_snapshot_resume %||% TRUE) {
  if (!requireNamespace("duckdb", quietly = TRUE)) {
    stop("config$gbif_source = 'snapshot' needs the duckdb package.")
  }
  date <- gbif_snapshot_resolve(cfg)
  bucket <- gbif_snapshot_bucket(cfg)

  ### Columns written: those the parse stage keeps (union over taxa) plus the
  ### columns it filters on. Snapshot column names are all lower case.
  keep <- unique(unlist(lapply(taxa, function(t) taxa_cfg[[t]]$gbif_columns_keep)))
  cols <- cfg$gbif_snapshot_columns %||%
    if (length(keep) > 0) unique(c(tolower(keep), "species", "specieskey", "year", "basisofrecord", "license"))
  cols_sql <- if (is.null(cols)) "*" else paste(cols, collapse = ", ")

  db <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(db, shutdown = TRUE), add = TRUE)
  q <- function(x) paste(DBI::dbQuoteString(db, x), collapse = ", ")

  ### 1. Plan.
  jobs <- stats::setNames(lapply(taxa, function(taxon) {
    tryCatch(
      plan_snapshot_taxon(taxon, cfg, taxa_cfg, include_flagged, max_age_days, force,
                          date, cols %||% "*"),
      error = function(e) list(
        taxon = taxon, state = "error", sp_names = NULL, out_dir = NA_character_,
        n_species = NA_integer_, n_flagged = NA_integer_, n_records = NA_integer_,
        error = conditionMessage(e)
      )
    )
  }), taxa)
  ids <- function(s) names(Filter(function(j) j$state == s, jobs))
  queued <- ids("queued")
  citation <- NA_character_

  if (length(queued) > 0) {
    citation <- gbif_snapshot_citation(date, cfg)
    f <- cfg$gbif_filters
    yrs <- range(cfg$years)
    bb <- region_bbox(cfg)
    all_keys <- unique(unlist(lapply(jobs[queued], `[[`, "sp_keys")))

    ### 2. One scan over the snapshot parts, in resumable chunks.
    scan_res <- tryCatch({
      invisible(DBI::dbExecute(db, "SET autoinstall_known_extensions=1"))
      invisible(DBI::dbExecute(db, "SET autoload_known_extensions=1"))
      invisible(DBI::dbExecute(db, glue("SET s3_region={q(sub('^gbif-open-data-', '', bucket))}")))
      invisible(DBI::dbExecute(db, glue("SET threads={as.integer(threads)}")))
      invisible(DBI::dbExecute(db, glue("SET memory_limit={q(memory_limit)}")))
      invisible(DBI::dbExecute(db, "SET preserve_insertion_order=false"))

      prefix <- glue("s3://{bucket}/occurrence/{date}/occurrence.parquet")
      parts <- sort(DBI::dbGetQuery(db, glue("SELECT file FROM glob('{prefix}/*')"))$file)
      if (length(parts) == 0) stop(glue("No parquet parts found under {prefix}."))
      if (!is.null(part_limit)) parts <- utils::head(parts, part_limit)

      needed <- unique(c(cols, "species", "specieskey", "year", "basisofrecord",
                         if (isTRUE(f$has_coordinate)) c("decimallatitude", "decimallongitude"),
                         if (!is.null(f$has_geo_issue)) "issue",
                         if (!is.null(cfg$region$country_code)) "countrycode",
                         if (!is.null(bb)) c("decimallatitude", "decimallongitude")))
      have <- DBI::dbGetQuery(
        db, glue("DESCRIBE SELECT * FROM read_parquet('{parts[[1]]}')"))$column_name
      if (length(setdiff(needed, have)) > 0) {
        stop(glue("Snapshot lacks columns: {paste(setdiff(needed, have), collapse = ', ')}"))
      }

      where <- c(
        glue("specieskey IN ({q(all_keys)})"),
        glue("year BETWEEN {as.integer(yrs[1])} AND {as.integer(yrs[2])}"),
        glue("basisofrecord IN ({q(f$basis_of_record)})"),
        if (isTRUE(f$has_coordinate))
          "decimallatitude IS NOT NULL AND decimallongitude IS NOT NULL",
        ### coalesce: a NULL issue array must count as "no issue", not drop the row.
        if (identical(f$has_geo_issue, FALSE))
          glue("NOT coalesce(list_has_any(issue, [{q(gbif_snapshot_geo_issues)}]), false)"),
        if (identical(f$has_geo_issue, TRUE))
          glue("coalesce(list_has_any(issue, [{q(gbif_snapshot_geo_issues)}]), false)"),
        if (!is.null(cfg$region$country_code))
          glue("countrycode = {q(cfg$region$country_code)}"),
        if (!is.null(bb))
          glue("decimallongitude BETWEEN {bb[['xmin']]} AND {bb[['xmax']]} ",
               "AND decimallatitude BETWEEN {bb[['ymin']]} AND {bb[['ymax']]}")
      )
      where_sql <- paste0("(", where, ")", collapse = " AND ")

      scan_hash <- rlang::hash(list(gbif_snapshot_hash(all_keys, date, cols %||% "*", cfg),
                                    part_limit))
      work_dir <- ensure_dir(file.path(cfg$paths$data_raw, "gbif", "_snapshot_work",
                                       glue("{date}_{substr(scan_hash, 1, 10)}")))
      if (!isTRUE(resume)) unlink(list.files(work_dir, full.names = TRUE))

      n_chunks <- min(as.integer(n_chunks), length(parts))
      chunks <- unname(split(parts, (seq_along(parts) - 1) %/% ceiling(length(parts) / n_chunks)))
      n_chunks <- length(chunks)
      gbif_log(glue("snapshot {date}: scanning {length(parts)} parts in {n_chunks} chunks ",
                    "for {length(all_keys)} specieskeys ({length(queued)} taxa); checkpoints in {work_dir}"))
      t0 <- Sys.time()
      for (i in seq_along(chunks)) {
        out <- file.path(work_dir, sprintf("chunk_%03d.parquet", i))
        if (file.exists(out)) {
          gbif_log(glue("chunk {i}/{n_chunks}: checkpoint exists; skipping."))
          next
        }
        tmp <- paste0(out, ".tmp")
        ### Retry transient S3 failures (503 / slow down / timeouts) with exponential
        ### backoff; anything else, or exhausting retries, fails the scan. The
        ### checkpoints already written are kept either way.
        max_tries <- cfg$gbif_snapshot_retries %||% 6
        for (attempt in seq_len(max_tries)) {
          unlink(tmp)
          res_try <- tryCatch({
            DBI::dbExecute(db, glue(
              "COPY (SELECT {cols_sql} FROM read_parquet([{paste(shQuote(chunks[[i]]), collapse = ', ')}]) ",
              "WHERE {where_sql}) TO '{tmp}' (FORMAT parquet)"))
            NULL
          }, error = function(e) e)
          if (is.null(res_try)) break
          transient <- grepl("503|SlowDown|Slow Down|timed? ?out|timeout|temporarily|connection|HTTP", 
                             conditionMessage(res_try), ignore.case = TRUE)
          if (!transient || attempt == max_tries) {
            unlink(tmp)
            stop(res_try)
          }
          wait <- min(30 * 2^(attempt - 1), 600)
          gbif_log(glue("chunk {i}/{n_chunks}: transient error (try {attempt}/{max_tries}): ",
                        "{substr(conditionMessage(res_try), 1, 150)} | retrying in {wait}s"))
          Sys.sleep(wait)
        }
        file.rename(tmp, out)
        el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
        gbif_log(glue("chunk {i}/{n_chunks} done | {round(el, 1)} min elapsed | ",
                      "ETA {round(el / i * (n_chunks - i), 1)} min"))
      }
      list(work_dir = work_dir)
    }, error = function(e) e)

    if (inherits(scan_res, "error")) {
      for (tx in queued) {
        jobs[[tx]]$state <- "error"
        jobs[[tx]]$error <- glue("snapshot scan failed: {conditionMessage(scan_res)} ",
                                 "(set config$gbif_source = 'download' to use the API instead)")
      }
      gbif_log(glue("snapshot scan FAILED: {conditionMessage(scan_res)}"))
    } else {
      ### 3. Split the scan into per-taxon datasets (staged, then swapped in).
      chunk_glob <- file.path(scan_res$work_dir, "chunk_*.parquet")
      for (tx in queued) {
        j <- jobs[[tx]]
        fin <- tryCatch({
          ensure_dir(j$out_dir)
          stage <- file.path(j$out_dir, ".snapshot_stage")
          unlink(stage, recursive = TRUE)
          dir.create(stage)
          DBI::dbExecute(db, glue(
            "COPY (SELECT * FROM read_parquet('{chunk_glob}') WHERE specieskey IN ({q(j$sp_keys)})) ",
            "TO '{file.path(stage, 'part-0.parquet')}' (FORMAT parquet)"))
          n_rec <- DBI::dbGetQuery(db, glue(
            "SELECT count(*) n FROM read_parquet('{file.path(stage, 'part-0.parquet')}')"))$n
          ### Crosswalk keys with zero records in the extract (no data, or a key that
          ### does not exist in this snapshot). Reported, not an error.
          got <- DBI::dbGetQuery(db, glue(
            "SELECT DISTINCT specieskey FROM read_parquet('{file.path(stage, 'part-0.parquet')}')"
          ))$specieskey
          zero_keys <- setdiff(j$sp_keys, got)
          if (length(zero_keys) > 0) {
            zero_nm <- unname(j$key_names[zero_keys])
            gbif_log(glue("{tx}: {length(zero_keys)}/{length(j$sp_keys)} crosswalk keys have ",
                          "ZERO records in the extract: ",
                          "{paste(utils::head(paste0(zero_keys, ' (', zero_nm, ')'), 15), collapse = '; ')}",
                          "{if (length(zero_keys) > 15) '; ...' else ''}"))
          } else {
            gbif_log(glue("{tx}: all {length(j$sp_keys)} crosswalk keys have records in the extract."))
          }
          unlink(j$dataset_path, recursive = TRUE)
          file.rename(stage, j$dataset_path)
          jsonlite::write_json(
            list(taxon = tx, source = "gbif_snapshot", snapshot_date = date,
                 snapshot_bucket = bucket, snapshot_citation = citation,
                 request_hash = j$req_hash, dataset = basename(j$dataset_path),
                 n_names = length(j$sp_names), n_keys = length(j$sp_keys),
                 keys_zero_records = as.list(zero_keys),
                 n_records = as.integer(n_rec),
                 downloaded_utc = format(Sys.time(), tz = "UTC", "%Y-%m-%d %H:%M:%S")),
            manifest_path(j$out_dir), auto_unbox = TRUE, pretty = TRUE
          )
          n_rec
        }, error = function(e) e)
        if (inherits(fin, "error")) {
          jobs[[tx]]$state <- "error"
          jobs[[tx]]$error <- conditionMessage(fin)
          gbif_log(glue("{tx}: extract FAILED: {conditionMessage(fin)}"))
        } else {
          jobs[[tx]]$state <- "done"
          jobs[[tx]]$n_records <- as.integer(fin)
          gbif_log(glue("{tx}: {fin} records written to {j$dataset_path}"))
        }
      }
      if (length(ids("error")) == 0 && !isTRUE(cfg$gbif_snapshot_keep_work %||% FALSE)) {
        unlink(scan_res$work_dir, recursive = TRUE)
      }
      gbif_log(glue("Cite the snapshot ({date}) and create a GBIF derived-dataset citation ",
                    "for the records used."))
    }
  }

  ### 4. One summary row per taxon.
  purrr::map_dfr(taxa, function(taxon) {
    j <- jobs[[taxon]]
    tibble::tibble(
      taxon = taxon,
      status = switch(j$state, skipped = "skipped", done = "ok", error = "error",
                      queued = "not_run"),
      n_species = j$n_species, n_flagged = j$n_flagged,
      n_names = length(j$sp_names), n_records = j$n_records %||% NA_integer_,
      source = "gbif_snapshot", snapshot_date = date,
      output_dir = j$out_dir, run_time_utc = as.character(Sys.time()),
      error_message = j$error %||% NA_character_
    )
  })
}
