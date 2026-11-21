### =============================================================================
### R/00_utils.R  —  small shared helpers used across pipeline stages
### =============================================================================
###
### The original project re-implemented these in several scripts (ensure_dir /
### dir.create in gbif_download_core, gbif_parse_core, clean_gbif_occs,
### parse_gbif_occs; remove_empty_files and standardize_species_name twice each).
### They live here once.

suppressPackageStartupMessages({
  library(stringr)
})

### Print `...` if `progress` is at least `level` ("terse" < "verbose").
.progress_msg <- function(progress, level, ...) {
  rank <- c(quiet = 0, terse = 1, verbose = 2)
  if (rank[[progress]] >= rank[[level]]) message(...)
  invisible()
}

### Stage-level cache check. TRUE when `force` is off, the manifest at `mpath`
### records the same `settings` hash, and every file in `outputs` still exists.
.stage_cached <- function(mpath, settings, outputs, force = FALSE) {
  if (isTRUE(force) || !file.exists(mpath)) return(FALSE)
  old <- tryCatch(jsonlite::read_json(mpath, simplifyVector = TRUE), error = function(e) NULL)
  !is.null(old) && identical(old$settings, settings) && all(file.exists(outputs))
}

### Write the manifest LAST, so an interrupted run is never treated as cached.
.write_stage_manifest <- function(mpath, settings, ...) {
  ensure_parent_dir(mpath)
  jsonlite::write_json(list(settings = settings, ...), mpath,
                       auto_unbox = TRUE, pretty = TRUE, na = "null")
}

### TRUE if a taxon's range_source denotes polygon (vector) range data.
### "polygon" is canonical; "shapefile" and "vector" are accepted aliases.
is_polygon_source <- function(range_source) {
  !is.null(range_source) && length(range_source) == 1 &&
    range_source %in% c("polygon", "shapefile", "vector")
}

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

### Path to a resolution-specific hex layer: inserts "_res<N>" before the
### extension, e.g. hex_grid.gpkg -> hex_grid_res4.gpkg, so grids at different
### H3 levels coexist.
hex_layer_path <- function(path, cfg = config) {
  sub("(\\.[^.]+)$", paste0("_res", cfg$hex_resolution, "\\1"), path)
}
hex_grid_path       <- function(cfg = config) hex_layer_path(cfg$paths$hex_grid, cfg)
cons_areas_hex_path <- function(cfg = config) hex_layer_path(cfg$paths$cons_areas_hex, cfg)

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

### Starting in version 3.9, rgbif defaults to Catalog of Life
### taxonomic key - alphanumeric speciesKey.  Ensure that user has at least 
### this version number.
check_package_version <- function(pkg = 'rgbif', min_ver = '3.9') {
  if (!requireNamespace(pkg, quietly = TRUE) || packageVersion(pkg) < min_ver) {
    stop('Requires {pkg} package version >= {min_ver}...')
  }
}