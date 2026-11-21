### =============================================================================
### R/02_acquire/02_download_obis.R  —  download OBIS occurrence records
### =============================================================================
###
### Marine counterpart of 01_download_gbif.R. For large taxon downloads,
### queries the OBIS GeoParquet dataset hosted on AWS S3 (via DuckDB) instead
### of paging through the API. For small downloads or when S3 access is unavailable,
### falls back to the OBIS API.
###
### All pending taxa share ONE combined S3 query (names unioned, results split by
### taxon); taxa the S3 query does not cover fall back to the API.
###
### Large downloads (>100k records estimated) use:
###   - AWS S3 GeoParquet files at s3://obis-open-data/occurrence/*.parquet
###   - DuckDB for spatial + attribute filtering (faster than API paging)
###   - Direct parquet write to taxa_config$obis_dataset_path
###
### Small downloads use the OBIS API (https://api.obis.org/v3/occurrence) with
### species-name matching, paged through 10k-record batches.
###
### Filters apply uniformly: year range (config$years), region bounding box,
### and basis_of_record (config$obis_filters). Columns are renamed to GBIF-style
### names so parse/clean stages treat both sources alike.
###
### No credentials needed. Skipped entirely when config$obis_filters$enabled is FALSE.
###
### Requires config + taxa_config + R/00_utils.R + R/02_acquire/01_download_gbif.R
### (for load_crosswalk() and region_bbox()) sourced first.
###
### Usage:
###   res <- run_obis_downloads()                 # all configured taxa
###   res <- run_obis_downloads(taxa = "mammals_marine", force = TRUE)
###   res <- run_obis_downloads(use_aws = FALSE)  # API only, no S3
###   res <- run_obis_downloads(match_by = "id")   # S3: AphiaID only ("name", "both")

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(glue)
  library(httr2)
  library(arrow)
  library(here)
  library(DBI)
})

obis_api_url <- "https://api.obis.org/v3/occurrence"
obis_s3_prefix <- "s3://obis-open-data/occurrence"

### Fields requested from OBIS (these are also the columns written to parquet).
obis_fields <- c(
  "id", "datasetID", "aphiaID", "scientificName", "species",
  "decimalLatitude", "decimalLongitude", "coordinateUncertaintyInMeters",
  "date_year", "month", "day", "eventDate",
  "basisOfRecord", "individualCount", "occurrenceStatus"
)

### Accepted binomials to query. Flagged (needs_review) rows are included by
### default; set include_flagged = FALSE to drop them. Species with no WoRMS
### aphia_id are dropped: they are not in OBIS and are handled by GBIF only.
get_obis_names <- function(crosswalk, include_flagged = TRUE) {
  nms <- dplyr::mutate(crosswalk, query_name = dplyr::coalesce(name_gbif, name_accepted))
  nms <- dplyr::filter(nms, !is.na(query_name), !is.na(aphia_id))
  if (!isTRUE(include_flagged)) nms <- dplyr::filter(nms, !needs_review)
  sort(unique(nms$query_name))
}

### WKT polygon for the region, or NULL when no spatial restriction is needed.
### Reuses region_bbox() from 01_download_gbif.R.
obis_region_wkt <- function(cfg = config) {
  bb <- region_bbox(cfg)
  if (is.null(bb)) return(NULL)
  as.character(glue::glue(
    "POLYGON(({bb[['xmin']]} {bb[['ymin']]}, {bb[['xmax']]} {bb[['ymin']]}, ",
    "{bb[['xmax']]} {bb[['ymax']]}, {bb[['xmin']]} {bb[['ymax']]}, ",
    "{bb[['xmin']]} {bb[['ymin']]}))"
  ))
}

### Accepted WoRMS AphiaIDs to query on S3 (same row filter as get_obis_names()).
### OBIS stores the accepted id, which is what the crosswalk's aphia_id holds.
get_obis_aphia_ids <- function(crosswalk, include_flagged = TRUE) {
  ids <- dplyr::filter(crosswalk, !is.na(aphia_id))
  if (!isTRUE(include_flagged)) ids <- dplyr::filter(ids, !needs_review)
  sort(unique(as.integer(ids$aphia_id)))
}

### --- AWS S3 GeoParquet path (new, faster path) -------------------------------------------
### Query OBIS GeoParquet on AWS S3 via DuckDB. Much faster for large datasets than API paging.
### Returns a data frame with OBIS fields (or NULL if DuckDB unavailable or query fails).
### match_by: "both" (default) = aphiaid in `aphia_ids` OR species in `species_names`;
### "id" = aphiaid only; "name" = species name only.
### The AphiaID arm catches renamed taxa; the name arm catches subspecies/synonym
### records whose aphiaID differs from the species-level id but whose `species` is rolled up.
fetch_obis_aws <- function(aphia_ids, species_names, year_range, bbox = NULL, cfg = config,
                           match_by = c("both", "name", "id")) {
  match_by <- match.arg(match_by)
  # DuckDB support is optional; return NULL if not available
  if (!requireNamespace("duckdb", quietly = TRUE)) {
    return(NULL)
  }
  
  tryCatch({
    db <- duckdb::dbConnect(duckdb::duckdb(), read_only = TRUE)
    on.exit(duckdb::dbDisconnect(db), add = TRUE)
    
    # Enable S3 access (auto-install and auto-load httpfs extension)
    DBI::dbExecute(db, "SET autoinstall_known_extensions=1")
    DBI::dbExecute(db, "SET autoload_known_extensions=1")
    
    ### The GeoParquet nests QC'd fields in the `interpreted` struct (see
    ### update_docs/08_AWS_OBIS_SCHEMA_DEBUG.md). Alias to the names in
    ### obis_fields so standardize_obis() treats S3 and API results alike.
    bbox_filter <- if (!is.null(bbox)) {
      glue::glue(
        "AND interpreted.decimalLongitude >= {bbox[['xmin']]} ",
        "AND interpreted.decimalLongitude <= {bbox[['xmax']]} ",
        "AND interpreted.decimalLatitude >= {bbox[['ymin']]} ",
        "AND interpreted.decimalLatitude <= {bbox[['ymax']]}"
      )
    } else ""

    ### Register ids and names as tables, matched via IN (subquery) below rather than
    ### large literal lists. (IN subqueries, unlike an OR join, stay hash-based.)
    duckdb::duckdb_register(
      db, "obis_ids",
      data.frame(aphia_id = unique(as.integer(aphia_ids)))
    )
    duckdb::duckdb_register(
      db, "obis_names",
      data.frame(query_name = unique(species_names), stringsAsFactors = FALSE)
    )
    t_start <- Sys.time()
    ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z", tz = "America/Los_Angeles")
    message(glue("[obis] ⏱️ {ts()} S3 query start ({length(unique(aphia_ids))} AphiaIDs, {length(unique(species_names))} names)"))

    ### Build the match clause from the chosen strategy.
    id_clause   <- "o.interpreted.aphiaid IN (SELECT aphia_id FROM obis_ids)"
    name_clause <- "o.interpreted.species IN (SELECT query_name FROM obis_names)"
    match_clause <- switch(match_by,
      both = paste0("(", id_clause, " OR ", name_clause, ")"),
      id   = id_clause,
      name = name_clause
    )

    sql <- glue::glue("
      SELECT
        CAST(_id AS VARCHAR) AS id,
        CAST(dataset_id AS VARCHAR) AS datasetID,
        interpreted.aphiaid AS aphiaID,
        interpreted.scientificName AS scientificName,
        interpreted.species AS species,
        interpreted.decimalLatitude AS decimalLatitude,
        interpreted.decimalLongitude AS decimalLongitude,
        interpreted.coordinateUncertaintyInMeters AS coordinateUncertaintyInMeters,
        interpreted.date_year AS date_year,
        interpreted.month AS month,
        interpreted.day AS day,
        interpreted.eventDate AS eventDate,
        interpreted.basisOfRecord AS basisOfRecord,
        interpreted.individualCount AS individualCount,
        interpreted.occurrenceStatus AS occurrenceStatus
      FROM read_parquet('{obis_s3_prefix}/*.parquet') AS o
      WHERE
        {match_clause}
        AND interpreted.date_year >= {as.integer(year_range[1])}
        AND interpreted.date_year <= {as.integer(year_range[2])}
        {bbox_filter}
    ")
    
    # Use DBI interface: dbGetQuery or dbExecute + dbFetch
    recs <- DBI::dbGetQuery(db, sql)
    message(glue("[obis] ⏱️ {ts()} S3 query done: {nrow(recs)} rows in ",
                 "{round(as.numeric(difftime(Sys.time(), t_start, units = 'mins')), 1)} min"))
    dplyr::as_tibble(recs)
  },
  error = function(e) {
    warning(glue("AWS S3 query failed: {conditionMessage(e)}; falling back to API."))
    NULL
  })
}

### --- Download manifest -------------------------------------------------------
### Same idea as the GBIF manifest: records what was downloaded and from which
### request, so unchanged requests are not re-pulled.
obis_manifest_path <- function(out_dir) file.path(out_dir, "obis_download_manifest.json")

### Hash of everything that determines the download's content.
obis_request_hash <- function(obis_names, cfg = config, aphia_ids = NULL,
                              match_by = "both") {
  rlang::hash(list(
    match_by = match_by,
    names   = obis_names,
    aphia   = aphia_ids,
    filters = cfg$obis_filters$basis_of_record,
    years   = range(cfg$years),
    region  = region_bbox(cfg),
    fields  = obis_fields
  ))
}

read_obis_manifest <- function(out_dir) {
  p <- obis_manifest_path(out_dir)
  if (!file.exists(p)) return(NULL)
  tryCatch(jsonlite::read_json(p, simplifyVector = TRUE), error = function(e) NULL)
}

### Returns NULL if the existing download is reusable, else a short reason string.
obis_download_stale_reason <- function(manifest, request_hash, out_file,
                                       max_age_days = Inf) {
  if (is.null(manifest)) return("no manifest")
  if (!identical(manifest$request_hash, request_hash)) return("request changed")
  if (!file.exists(out_file)) return("downloaded file missing")
  age <- as.numeric(difftime(Sys.time(), as.POSIXct(manifest$downloaded_utc, tz = "UTC"),
                             units = "days"))
  if (is.na(age) || age > max_age_days) {
    return(glue("older than max_age_days ({round(age, 1)} d)"))
  }
  NULL
}

### Rename to GBIF-style columns and apply the shared basis-of-record filter.
### OBIS uses CamelCase ("HumanObservation"); GBIF uses "HUMAN_OBSERVATION".
standardize_obis <- function(df, cfg = config) {
  if (nrow(df) == 0) return(df)
  ### Ensure every expected column exists even if no record populated it.
  for (col in setdiff(obis_fields, names(df))) df[[col]] <- NA
  df <- df |>
    dplyr::select(dplyr::all_of(obis_fields)) |>
    dplyr::rename(year = date_year, obisID = id) |>
    dplyr::mutate(
      basisOfRecord = toupper(gsub("([a-z])([A-Z])", "\\1_\\2", basisOfRecord)),
      year = as.integer(year)
    )
  keep <- cfg$obis_filters$basis_of_record
  if (!is.null(keep)) df <- dplyr::filter(df, basisOfRecord %in% keep)
  dplyr::distinct(df, obisID, .keep_all = TRUE)
}

### Page through the API for one batch of names. OBIS pages by an `after` cursor
### (the last record id), so we loop until an empty page comes back.
fetch_obis_batch <- function(names_batch, year_range, wkt = NULL, page_size = 10000) {
  query <- list(
    scientificname = paste(names_batch, collapse = ","),
    startdate = sprintf("%d-01-01", as.integer(year_range[[1]])),
    enddate   = sprintf("%d-12-31", as.integer(year_range[[2]])),
    size      = page_size,
    fields    = paste(obis_fields, collapse = ",")
  )
  if (!is.null(wkt)) query$geometry <- wkt

  pages <- list()
  after <- NULL
  repeat {
    q <- query
    if (!is.null(after)) q$after <- after
    body <- httr2::request(obis_api_url) |>
      httr2::req_url_query(!!!q) |>
      httr2::req_retry(max_tries = 4) |>
      httr2::req_perform() |>
      httr2::resp_body_json()
    res <- body$results
    if (length(res) == 0) break
    pages[[length(pages) + 1]] <- dplyr::bind_rows(lapply(res, tibble::as_tibble_row))
    after <- res[[length(res)]]$id
    if (length(res) < page_size) break
  }
  dplyr::bind_rows(pages)
}

### --- Per-taxon helpers --------------------------------------------------------
### Resolve one taxon's crosswalk rows into a download plan (names, paths, hash).
### Species without an aphia_id are excluded in get_obis_names() (GBIF-only).
prepare_obis_taxon <- function(taxon, cfg, taxa_cfg, include_flagged, match_by = "both") {
  tcfg <- taxa_cfg[[taxon]]
  if (is.null(tcfg)) stop(glue("Taxon '{taxon}' not in taxa_config."))
  if (is.null(tcfg$obis_dataset_path)) {
    stop(glue("Taxon '{taxon}' has no obis_dataset_path in taxa_config."))
  }

  crosswalk <- load_crosswalk(taxon, cfg)
  obis_names <- get_obis_names(crosswalk, include_flagged = include_flagged)
  if (length(obis_names) == 0) {
    stop(glue("No species with an aphia_id resolved for '{taxon}'."))
  }

  aphia_ids <- get_obis_aphia_ids(crosswalk, include_flagged = include_flagged)

  list(
    taxon = taxon, crosswalk = crosswalk, obis_names = obis_names, aphia_ids = aphia_ids,
    n_no_aphia = sum(is.na(crosswalk$aphia_id)),
    out_file = tcfg$obis_dataset_path,
    out_dir = dirname(tcfg$obis_dataset_path),
    req_hash = obis_request_hash(obis_names, cfg, aphia_ids, match_by)
  )
}

### One-row result tibble (a run-log entry).
obis_summary_row <- function(prep, status, n_records, source = NA_character_) {
  tibble::tibble(
    taxon = prep$taxon, status = status,
    n_species = nrow(prep$crosswalk), n_no_aphia = prep$n_no_aphia,
    n_names = length(prep$obis_names),
    n_records = n_records, output_file = prep$out_file,
    source = source, run_time_utc = as.character(Sys.time())
  )
}

obis_error_row <- function(taxon, e) {
  tibble::tibble(
    taxon = taxon, status = "error",
    n_species = NA_integer_, n_no_aphia = NA_integer_, n_names = NA_integer_,
    n_records = NA_integer_, output_file = NA_character_, source = NA_character_,
    run_time_utc = as.character(Sys.time()),
    error_message = conditionMessage(e)
  )
}

### Standardise, write parquet, then the manifest last (an interrupted run is never
### treated as cached).
write_obis_taxon <- function(prep, recs, source, cfg) {
  recs <- standardize_obis(recs, cfg)
  ensure_dir(prep$out_dir)
  arrow::write_parquet(recs, prep$out_file)
  jsonlite::write_json(
    list(taxon = prep$taxon, request_hash = prep$req_hash, file = basename(prep$out_file),
         n_names = length(prep$obis_names), n_records = nrow(recs), source = source,
         downloaded_utc = format(Sys.time(), tz = "UTC", "%Y-%m-%d %H:%M:%S")),
    obis_manifest_path(prep$out_dir), auto_unbox = TRUE, pretty = TRUE
  )
  obis_summary_row(prep, "ok", nrow(recs), source)
}

### Driver. Three steps:
###   1. Plan each taxon and skip those whose cached download is current.
###   2. ONE combined S3 query for the union of names across all pending taxa; the
###      result is split back to each taxon by its own name list.
###   3. Taxa the S3 query did not cover (query failed, DuckDB missing, or zero rows
###      for that taxon) fall back to the API, then every taxon is written.
### Per-taxon error capture so one failure doesn't abort the rest. Taxa without an
### `obis_dataset_path` are skipped by default (GBIF-only); pass them explicitly via
### `taxa` to get the error instead.
### max_age_days: re-download if the stored file is older than this
### (default: config$obis_max_age_days, else 180; Inf = never).
### use_aws: if TRUE, tries AWS S3 GeoParquet first before falling back to API.
### match_by: how S3 records are matched to species: "both" (default; AphiaID OR
### name), "name" (species name only) or "id" (AphiaID only). The API fallback can
### only query by name, so it ignores match_by (with a warning unless "both"/"name").
run_obis_downloads <- function(taxa = Filter(function(t) !is.null(taxa_cfg[[t]]$obis_dataset_path),
                                             cfg$taxa),
                               cfg = config, taxa_cfg = taxa_config,
                               include_flagged = TRUE,
                               max_age_days = cfg$obis_max_age_days %||% 180,
                               force = config$force,
                               use_aws = config$obis_use_aws,
                               batch_size = 100,
                               match_by = c("both", "name", "id")) {
  match_by <- match.arg(match_by)
  if (!source_enabled("obis", cfg)) {
    message("[obis] disabled (config$obis_filters$enabled = FALSE); skipping.")
    return(tibble::tibble())
  }

  ### 1. Plan + cache check -------------------------------------------------------
  planned <- purrr::map(rlang::set_names(taxa), function(taxon) {
    tryCatch({
      prep <- prepare_obis_taxon(taxon, cfg, taxa_cfg, include_flagged, match_by)
      if (prep$n_no_aphia > 0) {
        message(glue("[obis] {taxon}: skipping {prep$n_no_aphia} species with no aphia_id ",
                     "(GBIF-only)."))
      }
      if (!isTRUE(force)) {
        manifest <- read_obis_manifest(prep$out_dir)
        reason <- obis_download_stale_reason(manifest, prep$req_hash, prep$out_file, max_age_days)
        if (is.null(reason)) {
          message(glue("[obis] {taxon}: using CACHED download ",
                       "(downloaded {manifest$downloaded_utc} UTC; nothing re-fetched)."))
          return(list(row = obis_summary_row(prep, "skipped", manifest$n_records,
                                             manifest$source %||% "api")))
        }
        message(glue("[obis] {taxon}: downloading ({reason})."))
      }
      list(prep = prep)
    }, error = function(e) list(row = obis_error_row(taxon, e)))
  })
  rows <- purrr::compact(purrr::map(planned, "row"))
  pending <- purrr::compact(purrr::map(planned, "prep"))

  ### 2. One combined S3 query ----------------------------------------------------
  recs_by_taxon <- list()
  if (isTRUE(use_aws) && length(pending) > 0) {
    all_ids <- sort(unique(unlist(purrr::map(pending, "aphia_ids"))))
    message(glue("[obis] combined AWS S3 query: {length(all_ids)} AphiaIDs across ",
                 "{length(pending)} taxa ({paste(names(pending), collapse = ', ')})..."))
    all_names <- sort(unique(unlist(purrr::map(pending, "obis_names"))))
    aws <- fetch_obis_aws(all_ids, all_names, range(cfg$years), region_bbox(cfg), cfg, match_by)
    if (is.null(aws) || nrow(aws) == 0) {
      message("[obis] S3 query failed or returned 0 rows; all pending taxa will use the API.")
    } else {
      ### Split by each taxon's own AphiaIDs and/or names (per match_by); a match shared
      ### by two taxa goes to both.
      for (taxon in names(pending)) {
        hit_id   <- aws$aphiaID %in% pending[[taxon]]$aphia_ids
        hit_name <- aws$species %in% pending[[taxon]]$obis_names
        keep <- switch(match_by, both = hit_id | hit_name, id = hit_id, name = hit_name)
        sub <- aws[keep, ]
        message(glue("[obis] {taxon}: {nrow(sub)} records from AWS S3 GeoParquet"))
        if (nrow(sub) > 0) recs_by_taxon[[taxon]] <- sub
      }
    }
  }

  ### 3. API fallback where needed, then write each taxon -------------------------
  wkt <- obis_region_wkt(cfg)
  for (taxon in names(pending)) {
    prep <- pending[[taxon]]
    rows[[taxon]] <- tryCatch({
      recs <- recs_by_taxon[[taxon]]
      source <- "aws_geoparquet"
      if (is.null(recs)) {
        source <- "api"
        if (match_by == "id") {
          warning(glue("[obis] {taxon}: API fallback queries by name only; match_by = 'id' not honoured."))
        }
        message(glue("[obis] {taxon}: querying API (batch_size = {batch_size})..."))
        batches <- split(prep$obis_names, ceiling(seq_along(prep$obis_names) / batch_size))
        recs <- purrr::imap(batches, function(b, i) {
          message(glue("[obis] {taxon}: batch {i}/{length(batches)} ({length(b)} names)"))
          fetch_obis_batch(b, range(cfg$years), wkt)
        }) |>
          dplyr::bind_rows()
      }
      write_obis_taxon(prep, recs, source, cfg)
    }, error = function(e) obis_error_row(taxon, e))
  }

  dplyr::bind_rows(rows[intersect(taxa, names(rows))])
}

### Single-taxon entry point (kept for benchmark_downloads.R and ad-hoc use).
download_one_taxon_obis <- function(taxon, cfg = config, taxa_cfg = taxa_config,
                                    include_flagged = TRUE, max_age_days = Inf,
                                    force = config$force, use_aws = TRUE, batch_size = 100,
                                    match_by = "both") {
  run_obis_downloads(taxa = taxon, cfg = cfg, taxa_cfg = taxa_cfg,
                     include_flagged = include_flagged, max_age_days = max_age_days,
                     force = force, use_aws = use_aws, batch_size = batch_size,
                     match_by = match_by)
}
