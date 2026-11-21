### =============================================================================
### R/03_process/01_parse_occurrence.R  —  subset raw GBIF and/or OBIS into per-species-year parquet
### =============================================================================
###
### Reads each taxon's raw GBIF occurrence dataset and writes one parquet file per
### species-year under data/processed/species_populations/<taxon>/occurrences/<year>/.
###
### Sources are chosen in config (gbif_filters$enabled / obis_filters$enabled):
### GBIF only, OBIS only, or both. If a taxon has `obis_dataset_path`
### (taxa_config) and OBIS is enabled, the OBIS parquet from 02_download_obis.R is
### merged in. An occurrence found in both OBIS and GBIF (same
### coordinates to ~11 m and same date) is kept once, as the GBIF record. Output
### gains `source` ("obis"/"gbif") and `record_id` columns.
###
### Species come from the master crosswalk's `name_gbif` (the GBIF accepted
### binomial) — the single canonical name column — so the filter always matches
### GBIF's own `species` field regardless of taxon.
###
### Caching: after a taxon is parsed, a manifest (parse_occurrence_manifest_<taxon>.json,
### next to the run log) records fingerprints of the GBIF and OBIS download manifests,
### the run settings (years, species, columns, enabled sources) and the parsing code.
### If all of them match and the written files still exist, the taxon is skipped and
### the previous run log is returned. `force = TRUE` ignores the cache, reparses, and
### removes per-species-year files (and their _clean outputs) that are no longer produced.
###
### Requires config + taxa_config + R/00_utils.R + a built crosswalk.
###
### Usage:  res <- run_occurrence_parse()
###         res <- run_occurrence_parse(force = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(tidyr)
  library(glue)
  library(arrow)
  library(here)
})

### Case-insensitive resolution of requested columns to the dataset's columns.
resolve_columns_to_keep <- function(dataset, columns_to_keep) {
  available <- names(dataset)
  lookup <- setNames(available, tolower(available))
  resolved <- unname(lookup[tolower(columns_to_keep)])
  missing <- columns_to_keep[is.na(resolved)]
  if (length(missing) > 0) {
    stop("Missing GBIF columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  resolved
}

### Case-insensitive column accessor (GBIF parquet is lower-case, OBIS camelCase).
.col_ci <- function(df, name) {
  hit <- names(df)[tolower(names(df)) == tolower(name)]
  if (length(hit) == 0) rep(NA, nrow(df)) else df[[hit[[1]]]]
}

### Key used to recognise the same occurrence in OBIS and GBIF: coordinates
### (rounded to `digits` decimals, ~11 m at 4) + calendar date (first 10 chars of
### eventDate, falling back to year when no date). Species and year are already
### fixed by the caller, so they are not part of the key. The same record is
### usually published to both networks under different IDs, so IDs can't be used.
occ_dedup_key <- function(df, digits = 4) {
  date <- substr(as.character(.col_ci(df, "eventDate")), 1, 10)
  date <- dplyr::coalesce(date, as.character(.col_ci(df, "year")))
  paste(round(as.numeric(.col_ci(df, "decimalLatitude")), digits),
        round(as.numeric(.col_ci(df, "decimalLongitude")), digits),
        date, sep = "|")
}

### Align OBIS records to the GBIF columns kept (names and types) and tag source.
### OBIS ids go in the same `record_id` column GBIF ids are copied to.
.align_obis <- function(obis, gbif_cols, template) {
  out <- tibble::tibble(.rows = nrow(obis))
  for (col in gbif_cols) {
    v <- if (tolower(col) == "gbifid") rep(NA, nrow(obis)) else .col_ci(obis, col)
    ### Match the GBIF column's type (OBIS returns e.g. individualCount as text).
    tmpl <- template[[col]]
    out[[col]] <- if (is.numeric(tmpl)) as.numeric(v)
                  else if (is.character(tmpl)) as.character(v)
                  else v
  }
  out$record_id <- as.character(obis$obisID)
  out$source <- "obis"
  out
}

### eventDate is only used for the dedup key (first 10 chars), and GBIF (datetime)
### and OBIS (text, may be an interval) disagree on type, so keep it as YYYY-MM-DD text.
.as_event_date_chr <- function(df) {
  for (n in names(df)[tolower(names(df)) == "eventdate"]) {
    v <- df[[n]]
    df[[n]] <- if (inherits(v, "POSIXt")) format(v, "%Y-%m-%d", tz = "UTC")
               else substr(as.character(v), 1, 10)
  }
  df
}

### Merge OBIS into GBIF for one species-year. Where an OBIS and a GBIF record
### share a key (see occ_dedup_key), only the GBIF record is kept. Duplicates
### within a single source are left for cleaning (cc_dupl) to handle.
### Returns list(data, n_gbif, n_obis, n_dup).
merge_obis_gbif <- function(gbif, obis, gbif_cols, digits = 4) {
  ### GBIF disabled in config: OBIS only.
  if (is.null(gbif)) {
    # out <- .align_obis(obis, gbif_cols, template = list())
    out <- .as_event_date_chr(.align_obis(obis, gbif_cols, template = list()))
    return(list(data = out, n_gbif = 0L, n_obis = nrow(out), n_dup = 0L))
  }
  gbif$record_id <- as.character(.col_ci(gbif, "gbifID"))
  gbif$source <- rep("gbif", nrow(gbif))
  gbif <- .as_event_date_chr(gbif)
  n_obis <- if (is.null(obis)) 0L else nrow(obis)
  if (n_obis == 0) {
    return(list(data = gbif, n_gbif = nrow(gbif), n_obis = 0L, n_dup = 0L))
  }
  # obis <- .align_obis(obis, gbif_cols, gbif)
  obis <- .as_event_date_chr(.align_obis(obis, gbif_cols, gbif))
  dup <- occ_dedup_key(obis, digits) %in% occ_dedup_key(gbif, digits)
  list(data = dplyr::bind_rows(gbif, obis[!dup, , drop = FALSE]),
       n_gbif = nrow(gbif), n_obis = n_obis, n_dup = sum(dup))
}

### Write one species-year subset from already-collected GBIF and/or OBIS tables
### (NULL means that source is off). `gbif`/`obis` are pre-split to the target year
### by the caller; species is matched here, in memory. Returns a one-row status tibble.
write_occurrence_subset <- function(gbif, obis, target_species, target_year,
                                    columns_to_keep, output_dir) {
  if (!is.null(gbif)) {
    gbif <- gbif[!is.na(gbif$species) & gbif$species == target_species,
                 columns_to_keep, drop = FALSE]
  }
  if (!is.null(obis)) {
    hit <- (!is.na(obis$species) & obis$species == target_species) |
           (!is.na(obis$scientificName) & obis$scientificName == target_species)
    obis <- obis[hit, , drop = FALSE]
  }

  m <- merge_obis_gbif(gbif, obis, columns_to_keep)
  res <- m$data

  if (nrow(res) == 0) {
    return(tibble::tibble(species = target_species, year = target_year,
                          rows = 0L, gbif_rows = 0L, obis_rows = 0L, dup_dropped = 0L,
                          file_path = NA_character_, status = "no_data"))
  }

  year_dir <- ensure_dir(file.path(output_dir, as.character(target_year)))
  out_path <- file.path(year_dir,
                        paste0(standardize_species_name(target_species), "_occs.parquet"))
  arrow::write_parquet(res, sink = out_path)

  tibble::tibble(species = target_species, 
                 year    = target_year,
                 rows    = nrow(res), 
                 gbif_rows   = m$n_gbif, 
                 obis_rows   = m$n_obis,
                 dup_dropped = m$n_dup, 
                 ### ensure path is relative
                 file_path   = fs::path_rel(out_path, start = here()),
                 status      = "written")
}

### --- Parse cache -------------------------------------------------------------
parse_manifest_path <- function(taxon, cfg = config) {
  file.path(taxon_output_dir(taxon, cfg), "logs", glue("parse_occurrence_manifest_{taxon}.json"))
}

parse_log_path <- function(taxon, cfg = config) {
  file.path(taxon_output_dir(taxon, cfg), "logs", glue("parse_occurrence_{taxon}.csv"))
}

### Hash of a download manifest JSON ("none" if absent, so a missing manifest never
### matches a recorded one by accident).
.download_manifest_hash <- function(dir, file) {
  p <- file.path(dir, file)
  if (!file.exists(p)) return("none")
  m <- tryCatch(jsonlite::read_json(p, simplifyVector = TRUE), error = function(e) NULL)
  if (is.null(m)) "none" else rlang::hash(m)
}

### Hash of the parsing code. Deparsed (not hashed as closures) so comments,
### whitespace and line numbers don't invalidate the cache.
parse_code_hash <- function() {
  fns <- c("resolve_columns_to_keep", ".col_ci", "occ_dedup_key", ".align_obis",
           "merge_obis_gbif", "write_occurrence_subset", "parse_one_taxon")
  rlang::hash(lapply(fns, function(f) deparse(get(f, envir = globalenv()))))
}

### Fingerprint of everything that determines a taxon's parsed output.
parse_fingerprint <- function(tcfg, use_gbif, use_obis, species_vec, cfg = config) {
  list(
    gbif_manifest = if (use_gbif) .download_manifest_hash(dirname(tcfg$gbif_dataset_path),
                                                          "gbif_download_manifest.json") else "off",
    obis_manifest = if (use_obis) .download_manifest_hash(dirname(tcfg$obis_dataset_path),
                                                          "obis_download_manifest.json") else "off",
    run = rlang::hash(list(years    = sort(as.integer(cfg$years)),
                           species  = sort(species_vec),
                           columns  = tcfg$gbif_columns_keep)),
    code = parse_code_hash()
  )
}

### Returns NULL if the previous parse is reusable, else a short reason string.
parse_stale_reason <- function(taxon, fingerprint, cfg = config) {
  mp <- parse_manifest_path(taxon, cfg)
  if (!file.exists(mp)) return("no manifest")
  old <- tryCatch(jsonlite::read_json(mp, simplifyVector = TRUE), error = function(e) NULL)
  if (is.null(old)) return("unreadable manifest")
  for (nm in names(fingerprint)) {
    if (!identical(old$fingerprint[[nm]], fingerprint[[nm]])) return(glue("{nm} changed"))
  }
  lp <- parse_log_path(taxon, cfg)
  if (!file.exists(lp)) return("run log missing")
  log <- readr::read_csv(lp, show_col_types = FALSE)
  written <- log$file_path[log$status == "written"]
  if (!all(file.exists(here::here(written)))) return("output file missing")
  NULL
}

### Remove per-species-year files under `base_dir` that the current parse did not
### write (species dropped from the crosswalk, years no longer configured, or
### species-years that now have no data), plus their derived `_clean.parquet`.
### Returns the removed paths (invisibly).
remove_orphaned_occurrences <- function(base_dir, results) {
  kept <- normalizePath(here::here(results$file_path[results$status == "written"]),
                        mustWork = FALSE)
  on_disk <- list.files(base_dir, pattern = "_occs\\.parquet$", recursive = TRUE,
                        full.names = TRUE)
  orphans <- on_disk[!normalizePath(on_disk, mustWork = FALSE) %in% kept]
  derived <- sub("_occs\\.parquet$", "_clean.parquet", orphans)
  gone <- c(orphans, derived[file.exists(derived)])
  if (length(gone) > 0) file.remove(gone)
  invisible(gone)
}

### Parse one taxon across all configured years. Skips when the cache is current
### (see parse_stale_reason) unless `force = TRUE`.
parse_one_taxon <- function(taxon, cfg = config, taxa_cfg = taxa_config, force = config$force) {
  tcfg <- taxa_cfg[[taxon]]
  dataset_path <- tcfg$gbif_dataset_path

  ### Sources are switched with config$gbif_filters$enabled / obis_filters$enabled.
  use_gbif <- source_enabled("gbif", cfg)
  use_obis <- source_enabled("obis", cfg) && !is.null(tcfg$obis_dataset_path)
  if (!use_gbif && !use_obis) {
    stop(glue("No occurrence source available for '{taxon}': GBIF is disabled and ",
              "OBIS is disabled or has no obis_dataset_path."))
  }

  dataset <- NULL
  if (use_gbif) {
    if (!file.exists(dataset_path) && !dir.exists(dataset_path)) {
      stop(glue("GBIF dataset not found for '{taxon}' at {dataset_path}. ",
                "Run download_gbif first, or point taxa_config at your data."))
    }
    ### GBIF partitions are extensionless (000000, 000001, ...), and one is zero-byte.
    remove_empty_files(dataset_path, pattern = "\\.parquet$|^[0-9]+$")
    dataset <- arrow::open_dataset(dataset_path)
  }

  ### Optional OBIS source (marine taxa). Missing file is an error, not a skip,
  ### so a configured-but-undownloaded OBIS dataset is never silently ignored.
  obis_dataset <- NULL
  if (use_obis) {
    if (!file.exists(tcfg$obis_dataset_path)) {
      stop(glue("OBIS dataset not found for '{taxon}' at {tcfg$obis_dataset_path}. ",
                "Run download_obis first, or remove obis_dataset_path from taxa_config."))
    }
    obis_dataset <- arrow::open_dataset(tcfg$obis_dataset_path)
  }

  species_vec <- readr::read_csv(cfg$paths$crosswalk, show_col_types = FALSE) |>
    dplyr::filter(taxon == !!taxon, !is.na(name_gbif)) |>
    dplyr::pull(name_gbif) |>
    stringr::str_squish() |>
    unique()

  fingerprint <- parse_fingerprint(tcfg, use_gbif, use_obis, species_vec, cfg)
  if (!force) {
    reason <- parse_stale_reason(taxon, fingerprint, cfg)
    if (is.null(reason)) {
      message(glue("[{taxon}] parse cache is current; skipping (force = TRUE to reparse)"))
      return(readr::read_csv(parse_log_path(taxon, cfg), show_col_types = FALSE) |>
               dplyr::select(-dplyr::any_of("taxon")))
    }
    message(glue("[{taxon}] reparsing: {reason}"))
  }

  cols <- if (use_gbif) resolve_columns_to_keep(dataset, tcfg$gbif_columns_keep) else tcfg$gbif_columns_keep
  base_dir <- file.path(taxon_processed_dir(taxon, cfg), "occurrences")

  years <- as.integer(cfg$years)

  ### Scan each source once, restricted to this taxon's species and years, then
  ### split by year so each species-year lookup is an in-memory subset.
  gbif_by_year <- NULL
  if (use_gbif) {
    gbif_by_year <- dataset |>
      dplyr::filter(species %in% species_vec, year %in% years) |>
      dplyr::select(dplyr::all_of(union(cols, c("species", "year")))) |>
      dplyr::collect() |>
      (\(d) split(d, factor(d$year, levels = years)))()
  }
  obis_by_year <- NULL
  if (use_obis) {
    obis_by_year <- obis_dataset |>
      dplyr::filter(species %in% species_vec | scientificName %in% species_vec,
                    year %in% years) |>
      dplyr::collect() |>
      (\(d) split(d, factor(d$year, levels = years)))()
  }

  grid <- tidyr::expand_grid(species = species_vec, year = years)
  results <- purrr::pmap_dfr(
    list(grid$species, grid$year),
    function(species, year) {
      yr <- as.character(year)
      write_occurrence_subset(
        gbif_by_year[[yr]], obis_by_year[[yr]], species, year,
        columns_to_keep = cols, output_dir = base_dir
      )
    }
  )

  if (force) {
    gone <- remove_orphaned_occurrences(base_dir, results)
    if (length(gone) > 0) message(glue("[{taxon}] removed {length(gone)} orphaned file(s)"))
  }

  ### Per-taxon run log, then the manifest last so an interrupted run is never cached.
  log_path <- parse_log_path(taxon, cfg)
  ensure_parent_dir(log_path)
  readr::write_csv(dplyr::mutate(results, taxon = taxon), log_path)
  jsonlite::write_json(
    list(taxon = taxon, fingerprint = fingerprint,
         parsed_utc = format(Sys.time(), tz = "UTC", usetz = FALSE)),
    parse_manifest_path(taxon, cfg), auto_unbox = TRUE, pretty = TRUE
  )
  results
}

run_occurrence_parse <- function(taxa = config$taxa, cfg = config, taxa_cfg = taxa_config,
                                 force = config$force) {
  purrr::map_dfr(taxa, function(taxon) {
    parse_one_taxon(taxon, cfg, taxa_cfg, force = force) |> dplyr::mutate(taxon = taxon, .before = 1)
  })
}
