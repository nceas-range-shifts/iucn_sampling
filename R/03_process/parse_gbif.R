### =============================================================================
### R/03_process/parse_gbif.R  —  subset raw GBIF into per-species-year parquet
### =============================================================================
###
### Reads each taxon's raw GBIF occurrence dataset and writes one parquet file per
### species-year under data/processed/species_populations/<taxon>/gbif/<year>/.
###
### Species come from the master crosswalk's `name_gbif` (the GBIF accepted
### binomial) — the single canonical name column — so the filter always matches
### GBIF's own `species` field regardless of taxon.
###
### Requires config + taxa_config + R/00_utils.R + a built crosswalk.
###
### Usage:  res <- run_gbif_parse()

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

### Write one species-year subset; returns a one-row status tibble.
write_gbif_subset <- function(dataset, target_species, target_year,
                              species_col, year_col, columns_to_keep, output_dir) {
  res <- dataset |>
    dplyr::filter(.data[[year_col]] == target_year,
                  .data[[species_col]] == target_species) |>
    dplyr::select(dplyr::all_of(columns_to_keep)) |>
    dplyr::collect()

  if (nrow(res) == 0) {
    return(tibble::tibble(species = target_species, year = target_year,
                          rows = 0L, file_path = NA_character_, status = "no_data"))
  }

  year_dir <- ensure_dir(file.path(output_dir, as.character(target_year)))
  out_path <- file.path(year_dir,
                        paste0(standardize_species_name(target_species), "_occs.parquet"))
  arrow::write_parquet(res, sink = out_path)

  tibble::tibble(species = target_species, year = target_year,
                 rows = nrow(res), file_path = out_path, status = "written")
}

### Parse one taxon across all configured years.
parse_one_taxon <- function(taxon, cfg = config, taxa_cfg = taxa_config) {
  tcfg <- taxa_cfg[[taxon]]
  dataset_path <- tcfg$gbif_dataset_path
  if (!file.exists(dataset_path) && !dir.exists(dataset_path)) {
    stop(glue("GBIF dataset not found for '{taxon}' at {dataset_path}. ",
              "Run download_gbif first, or point taxa_config at your data."))
  }

  remove_empty_files(dataset_path)
  dataset <- arrow::open_dataset(dataset_path)

  species_vec <- readr::read_csv(cfg$paths$crosswalk, show_col_types = FALSE) |>
    dplyr::filter(taxon == !!taxon, !is.na(name_gbif)) |>
    dplyr::pull(name_gbif) |>
    stringr::str_squish() |>
    unique()

  cols <- resolve_columns_to_keep(dataset, tcfg$gbif_columns_keep)
  base_dir <- file.path(taxon_processed_dir(taxon, cfg), "gbif")

  grid <- tidyr::expand_grid(species = species_vec, year = as.integer(cfg$years))
  results <- purrr::pmap_dfr(
    list(grid$species, grid$year),
    function(species, year) write_gbif_subset(
      dataset, species, year,
      species_col = "species", year_col = "year",
      columns_to_keep = cols, output_dir = base_dir
    )
  )

  ### Per-taxon run log.
  log_path <- file.path(taxon_output_dir(taxon, cfg), "logs",
                        glue("parse_gbif_{taxon}.csv"))
  ensure_parent_dir(log_path)
  readr::write_csv(dplyr::mutate(results, taxon = taxon), log_path)
  results
}

run_gbif_parse <- function(taxa = config$taxa, cfg = config, taxa_cfg = taxa_config) {
  purrr::map_dfr(taxa, function(taxon) {
    parse_one_taxon(taxon, cfg, taxa_cfg) |> dplyr::mutate(taxon = taxon, .before = 1)
  })
}
