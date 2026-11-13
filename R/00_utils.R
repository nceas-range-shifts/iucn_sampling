### =============================================================================
### R/00_utils.R  —  small shared helpers used across pipeline stages
### =============================================================================
#
### The original project re-implemented these in several scripts (ensure_dir /
### dir.create in gbif_download_core, gbif_parse_core, clean_gbif_occs,
### parse_gbif_occs; remove_empty_files and standardize_species_name twice each).
### They live here once.

suppressPackageStartupMessages({
  library(stringr)
})

### Create a directory (and parents) if it does not already exist. Returns `path`
### invisibly so it can be used inline.
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

### Ensure the *directory containing* a file path exists, then return the path.
ensure_parent_dir <- function(file_path) {
  ensure_dir(dirname(file_path))
  invisible(file_path)
}

### Turn a scientific name into a safe file-name token: lower-case, non-alnum ->
### underscore, collapsed and trimmed. "Branta canadensis" -> "branta_canadensis".
standardize_species_name <- function(x) {
  x |>
    stringr::str_to_lower() |>
    stringr::str_replace_all("[^a-z0-9]+", "_") |>
    stringr::str_replace_all("^_+|_+$", "")
}

### Remove zero-byte files left behind by an interrupted write. Returns the paths
### that were removed (invisibly).
remove_empty_files <- function(dir, pattern = "\\.parquet$") {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) return(invisible(character(0)))
  empty <- files[file.info(files)$size == 0]
  if (length(empty) > 0) file.remove(empty)
  invisible(empty)
}

### Canonical singular taxon slug used in every path (bird/reptile/mammal). Kept
### as a function so any future aliasing lives in one place.
taxon_slug <- function(taxon) tolower(taxon)
