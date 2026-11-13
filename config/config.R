# =============================================================================
# config/config.R  —  GLOBAL project configuration
# =============================================================================
#
# This is the first place to edit when adapting the pipeline to your own study.
# Everything here is read by the pipeline stages via `config`; no analysis
# script hardcodes a region, year, CRS, or path of its own.
#
# Load order (see R/run_pipeline.R):
#   config/config.R  ->  config/taxa_config.R  ->  config/sampling_config.R
#
# Paths are all built with here::here() so the working directory only has to be
# the project root (open the .Rproj, or set here::i_am() once).

library(here)

config <- list(

  # --- Study region ----------------------------------------------------------
  # `admin_level`/`admin_name` select the focal area from a boundary layer
  # (e.g. a US state via tigris/rnaturalearth). `country_code` scopes GBIF/eBird
  # downloads. Swap these to move the study to a different place.
  ### NOTE: this region is used for spatial filtering in OBIS/GBIF.  Use global bbox
  region = list(
    country_code = NULL,
    admin_level  = NULL,
    admin_name   = "Global Oceans",  # Label only
    label        = "global_oceans",  # slug for output filenames
    # Region boundary: manually-defined bounding box (or path to polygon)
    boundary_source = "bbox",
    # Bounding box: c(xmin, ymin, xmax, ymax) in EPSG:4326
    bbox = c(-180, 180, -90, 90)
    # Alternative: point to a shapefile/geopackage
  ),

  # --- Conservation-areas layer ----------------------------------------------
  # A protected-areas polygon layer for the region. Field names are configurable
  # so you can point this at your own dataset without editing code. The defaults
  # match the California Protected Areas / GAP-status schema used in the study.
  conservation_areas = list(
    path          = here("data", "raw", "wdpa_2026", "WDPA_Nov2026_Public.gdb"),
    gap_field     = "IUCN_CAT",            # protection-status field
    gap_values    = c('Ia', 'Ib', 'II', 'III', 'IV'), # statuses to keep (strict protection)
    realm_field   = "MARINE",          # terrestrial/marine field (NULL to skip)
    realm_value   = c(1, 2),
    unit_field    = "NAME",        # protected-area name
    agency_field  = "MANG_AUTH",        # managing agency
    agency_level_field = "DESIG_TYPE"
  ),

  # --- Years -----------------------------------------------------------------
  # `years` is every year the pipeline processes; `focal_year` is the single
  # year used where one snapshot is needed (e.g. current-condition indicators).
  # Comparisons (range dynamics) use `compare_years`.
  years        = c(2014, 2015, 2021, 2022, 2023, 2024),
  focal_year   = 2024,
  compare_years = c(2014, 2024),

  # --- Coordinate reference systems ------------------------------------------
  # Named so scripts never sprinkle bare EPSG codes. `working` is lon/lat for
  # H3 indexing; `equal_area` for area/summary stats; `display` for basemaps.
  crs = list(
    working    = "EPSG:4326",   # H3 / geographic
    equal_area = "EPSG:6933",     # Cylindrical Equal Area (global, good for oceans)
    display    = "EPSG:3857"    # web mercator, plotting only
  ),

  # --- Spatial resolution ----------------------------------------------------
  ### level 4: ~1440 km^2; 5: ~250 km^2; 6: ~36 km^2
  hex_resolution = 4,           # H3 resolution for the analysis grid

  # --- GBIF download filters -------------------------------------------------
  # Shared occurrence-download predicates (region is added from `region` above).
  gbif_filters = list(
    has_geo_issue   = FALSE,
    has_coordinate  = TRUE,
    basis_of_record = c("HUMAN_OBSERVATION", "OBSERVATION",
                        "MACHINE_OBSERVATION", "PRESERVED_SPECIMEN"),
    download_format = "SIMPLE_PARQUET"
  ),

  # --- Paths (roots only; per-taxon/year paths are derived from these) -------
  paths = list(
    data_raw       = here("data", "raw"),
    data_processed = here("data", "processed"),
    output         = here("output"),
    # The single authoritative species crosswalk written by Stage 0.
    crosswalk      = here("data", "processed", "species_crosswalk.csv"),
    # Shared spatial layers produced by 03_process/build_hex_grid.R
    hex_grid       = here("data", "processed", "spatial_layers", "hex_grid.gpkg"),
    cons_areas_hex = here("data", "processed", "spatial_layers", "cons_areas_hex.gpkg")
  ),

  # --- Occurrence-cleaning thresholds ----------------------------------------
  # Used by 03_process/clean_gbif.R (CoordinateCleaner + quality filters).
  cleaning = list(
    max_coord_uncertainty_m = 750,   # drop records less precise than this
    min_coord_precision     = 0.01,  # drop coarser-than-this precision
    max_individual_count    = 99,    # drop implausibly large counts
    exclude_basis_of_record = c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN")
  ),

  # --- Default species-list source -------------------------------------------
  # By default the template runs on the committed sample list. Point this at
  # your own CSV (a `scientific_name` column is enough) to use real data. A
  # Google Sheet can be used per-taxon instead via config/taxa_config.R.
  species_list_default = here("data/raw/iucn_2026/iucn_species_list.csv"),

  ### --- Which taxa to run -------------------------------------------------
  ### Must match the names defined in config/taxa_config.R.
  taxa = c("mammals_marine")
)

# --- Output-path helpers (ONE canonical scheme) ------------------------------
# The original project mixed singular/plural taxon folders and an inconsistent
# `indicators/` segment. The template uses ONE scheme everywhere:
#   output/indicators/species_populations/<taxon>/...
# with the SINGULAR taxon slug (bird/reptile/mammal). Everything derives paths
# from these helpers so producer and consumer can never disagree.

taxon_output_dir <- function(taxon, cfg = config) {
  file.path(cfg$paths$output, "indicators", "species_populations", taxon)
}

taxon_processed_dir <- function(taxon, cfg = config) {
  file.path(cfg$paths$data_processed, "species_populations", taxon)
}

# Per-species-year parquet directory (parsed/cleaned GBIF occurrences).
taxon_year_occ_dir <- function(taxon, year, cfg = config) {
  file.path(taxon_processed_dir(taxon, cfg), "gbif", as.character(year))
}
