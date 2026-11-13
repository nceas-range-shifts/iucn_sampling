### =============================================================================
### R/03_process/parse_occurrence.R  —  subset raw GBIF and/or OBIS into per-species-year parquet
### =============================================================================
###
### Reads each taxon's raw GBIF occurrence dataset and writes one parquet file per
### species-year under data/processed/species_populations/<taxon>/occurrences/<year>/.
###
### Sources are chosen in config (gbif_filters$enabled / obis_filters$enabled):
### GBIF only, OBIS only, or both. If a taxon has `obis_dataset_path`
### (taxa_config) and OBIS is enabled, the OBIS parquet from download_obis.R is
### merged in. An occurrence found in both OBIS and GBIF (same
### coordinates to ~11 m and same date) is kept once, as the OBIS record. Output
### gains `source` ("obis"/"gbif") and `record_id` columns.
###
### Species come from the master crosswalk's `name_gbif` (the GBIF accepted
### binomial) — the single canonical name column — so the filter always matches
### GBIF's own `species` field regardless of taxon.
###
### Requires config + taxa_config + R/00_utils.R + a built crosswalk.
###
### Usage:  res <- run_occurrence_parse()

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

### Merge OBIS into GBIF for one species-year. Where an OBIS and a GBIF record
### share a key (see occ_dedup_key), only the OBIS record is kept. Duplicates
### within a single source are left for cleaning (cc_dupl) to handle.
### Returns list(data, n_gbif, n_obis, n_dup).
merge_obis_gbif <- function(gbif, obis, gbif_cols, digits = 4) {
  ### GBIF disabled in config: OBIS only.
  if (is.null(gbif)) {
    out <- .align_obis(obis, gbif_cols, template = list())
    return(list(data = out, n_gbif = 0L, n_obis = nrow(out), n_dup = 0L))
  }
  gbif$record_id <- as.character(.col_ci(gbif, "gbifID"))
  gbif$source <- rep("gbif", nrow(gbif))
  n_obis <- if (is.null(obis)) 0L else nrow(obis)
  if (n_obis == 0) {
    return(list(data = gbif, n_gbif = nrow(gbif), n_obis = 0L, n_dup = 0L))
  }
  obis <- .align_obis(obis, gbif_cols, gbif)
  dup <- occ_dedup_key(gbif, digits) %in% occ_dedup_key(obis, digits)
  list(data = dplyr::bind_rows(obis, gbif[!dup, , drop = FALSE]),
       n_gbif = nrow(gbif), n_obis = n_obis, n_dup = sum(dup))
}

### Write one species-year subset from GBIF and/or OBIS (a NULL dataset means that
### source is off); returns a one-row status tibble.
write_occurrence_subset <- function(dataset, target_species, target_year,
                              species_col, year_col, columns_to_keep, output_dir,
                              obis_dataset = NULL) {
  gbif <- NULL
  if (!is.null(dataset)) {
    gbif <- dataset |>
      dplyr::filter(.data[[year_col]] == target_year,
                    .data[[species_col]] == target_species) |>
      dplyr::select(dplyr::all_of(columns_to_keep)) |>
      dplyr::collect()
  }

  obis <- NULL
  if (!is.null(obis_dataset)) {
    obis <- obis_dataset |>
      dplyr::filter(year == target_year,
                    species == target_species | scientificName == target_species) |>
      dplyr::collect()
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

  tibble::tibble(species = target_species, year = target_year,
                 rows = nrow(res), gbif_rows = m$n_gbif, obis_rows = m$n_obis,
                 dup_dropped = m$n_dup, file_path = out_path, status = "written")
}

### Parse one taxon across all configured years.
parse_one_taxon <- function(taxon, cfg = config, taxa_cfg = taxa_config) {
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
    remove_empty_files(dataset_path)
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

  cols <- if (use_gbif) resolve_columns_to_keep(dataset, tcfg$gbif_columns_keep) else tcfg$gbif_columns_keep
  base_dir <- file.path(taxon_processed_dir(taxon, cfg), "occurrences")

  grid <- tidyr::expand_grid(species = species_vec, year = as.integer(cfg$years))
  results <- purrr::pmap_dfr(
    list(grid$species, grid$year),
    function(species, year) write_occurrence_subset(
      dataset, species, year,
      species_col = "species", year_col = "year",
      columns_to_keep = cols, output_dir = base_dir,
      obis_dataset = obis_dataset
    )
  )

  ### Per-taxon run log.
  log_path <- file.path(taxon_output_dir(taxon, cfg), "logs",
                        glue("parse_occurrence_{taxon}.csv"))
  ensure_parent_dir(log_path)
  readr::write_csv(dplyr::mutate(results, taxon = taxon), log_path)
  results
}

run_occurrence_parse <- function(taxa = config$taxa, cfg = config, taxa_cfg = taxa_config) {
  purrr::map_dfr(taxa, function(taxon) {
    parse_one_taxon(taxon, cfg, taxa_cfg) |> dplyr::mutate(taxon = taxon, .before = 1)
  })
}
