# =============================================================================
# data/sample/make_sample_data.R  —  synthesize the committed sample dataset
# =============================================================================
#
# Generates tiny, self-contained inputs so the pipeline runs end-to-end offline
# (no API keys, no large downloads) for reptiles + mammals:
#
#   species_list.csv                     scientific_name, taxon
#   species_crosswalk.csv                pre-built crosswalk (skips GBIF/eBird API)
#   spatial_layers/region.gpkg           small study-region polygon
#   spatial_layers/conservation_areas.gpkg  a couple of protected polygons
#   ranges/<taxon>/<species>.shp         one range polygon per species (SName field)
#   gbif/<taxon>/occurrence.parquet      synthetic occurrence points (2014, 2024)
#
# Run once from the project root:  Rscript data/sample/make_sample_data.R
# Committed outputs are force-added past .gitignore (see .gitignore tail).

suppressPackageStartupMessages({
  library(sf); library(terra); library(dplyr); library(arrow); library(here)
})
set.seed(42)

root <- here("data", "sample")
dir.create(file.path(root, "spatial_layers"), recursive = TRUE, showWarnings = FALSE)

# --- Species -----------------------------------------------------------------
species <- tibble::tribble(
  ~scientific_name,             ~taxon,
  "Sceloporus occidentalis",    "reptile",
  "Elgaria multicarinata",      "reptile",
  "Lampropeltis californiae",   "reptile",
  "Sylvilagus bachmani",        "mammal",
  "Otospermophilus beecheyi",   "mammal",
  "Odocoileus hemionus",        "mammal"
)
readr::write_csv(species, file.path(root, "species_list.csv"))

# --- Study region (a ~1x1 degree box in central California) ------------------
region_box <- function(xmin, ymin, xmax, ymax) {
  st_sf(name = "sample_region",
        geometry = st_sfc(st_polygon(list(rbind(
          c(xmin, ymin), c(xmax, ymin), c(xmax, ymax),
          c(xmin, ymax), c(xmin, ymin)))), crs = 4326))
}
region <- region_box(-120, 35, -119, 36)
st_write(region, file.path(root, "spatial_layers", "region.gpkg"), delete_dsn = TRUE, quiet = TRUE)

# --- Conservation areas (two protected polygons inside the region) -----------
cons <- st_sf(
  UNIT_NAME  = c("Sample Reserve A", "Sample Reserve B"),
  reGAP      = c(1L, 2L),
  TerrMar    = c("Terrestrial", "Terrestrial"),
  MNG_AGNCY  = c("Sample Agency", "Sample Agency"),
  MNG_AG_LEV = c("State", "Federal"),
  geometry = st_sfc(
    st_polygon(list(rbind(c(-119.9, 35.1), c(-119.6, 35.1),
                          c(-119.6, 35.4), c(-119.9, 35.4), c(-119.9, 35.1)))),
    st_polygon(list(rbind(c(-119.5, 35.5), c(-119.2, 35.5),
                          c(-119.2, 35.8), c(-119.5, 35.8), c(-119.5, 35.5)))),
    crs = 4326)
)
st_write(cons, file.path(root, "spatial_layers", "conservation_areas.gpkg"),
         delete_dsn = TRUE, quiet = TRUE)

# --- Per-species range shapefiles + synthetic occurrences --------------------
# Each species gets a random sub-box of the region as its range; occurrences are
# points sampled inside that range across two years.
first_word <- function(x) vapply(strsplit(x, "\\s+"), `[`, character(1), 1)

gbif_cols <- function(sp_range, sp_name, n, yr) {
  bb <- st_bbox(sp_range)
  lon <- runif(n, bb["xmin"], bb["xmax"]); lat <- runif(n, bb["ymin"], bb["ymax"])
  tibble::tibble(
    gbifID = paste0("g", sample(1e6:9e6, n)),
    datasetKey = "sample", occurrenceID = paste0("o", seq_len(n)),
    family = "Sample", genus = first_word(sp_name), species = sp_name,
    scientificName = sp_name, occurrenceStatus = "PRESENT", publishingOrgKey = "sample",
    decimalLatitude = lat, decimalLongitude = lon,
    coordinateUncertaintyInMeters = sample(10:200, n, replace = TRUE),
    coordinatePrecision = NA_real_, year = yr,   # usually absent in real GBIF
    month = sample(1:12, n, replace = TRUE), day = sample(1:28, n, replace = TRUE),
    eventDate = as.character(as.Date(paste0(yr, "-01-01"))),
    basisOfRecord = "HUMAN_OBSERVATION", individualCount = 1L
  )
}

crosswalk_rows <- list()
for (tx in c("reptile", "mammal")) {
  dir.create(file.path(root, "ranges", tx), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "gbif", tx), recursive = TRUE, showWarnings = FALSE)
  sp_tx <- species$scientific_name[species$taxon == tx]

  occ_all <- list()
  for (sp_name in sp_tx) {
    x0 <- runif(1, -119.95, -119.5); y0 <- runif(1, 35.05, 35.5)
    sp_range <- st_sf(SName = sp_name,
      geometry = st_sfc(st_polygon(list(rbind(
        c(x0, y0), c(x0 + 0.4, y0), c(x0 + 0.4, y0 + 0.4),
        c(x0, y0 + 0.4), c(x0, y0)))), crs = 4326))
    shp <- file.path(root, "ranges", tx, paste0(gsub(" ", "_", sp_name), ".shp"))
    st_write(sp_range, shp, delete_dsn = TRUE, quiet = TRUE)

    occ_all[[sp_name]] <- dplyr::bind_rows(
      gbif_cols(sp_range, sp_name, 40, 2014L),
      gbif_cols(sp_range, sp_name, 55, 2024L))

    crosswalk_rows[[sp_name]] <- tibble::tibble(
      taxon = tx, name_input = sp_name, name_accepted = sp_name, name_gbif = sp_name,
      gbif_usage_key = sample(1e6:9e6, 1), ebird_code = NA, ebird_name = NA,
      common_name = NA, class = NA, taxon_group = tx, gbif_status = "ACCEPTED",
      match_type = "EXACT", verifier_status = NA, verifier_source = NA,
      name_changed = FALSE, needs_review = FALSE, review_reason = "", override = FALSE,
      range_shapefile_name = sp_name, range_shapefile_path = shp)
  }
  arrow::write_parquet(dplyr::bind_rows(occ_all),
                       file.path(root, "gbif", tx, "occurrence.parquet"))
}

readr::write_csv(dplyr::bind_rows(crosswalk_rows),
                 file.path(root, "species_crosswalk.csv"))

message("Sample dataset written to data/sample/. Run the pipeline with:\n",
        "  source(here::here('R','run_pipeline.R')); run_pipeline(sample = TRUE)")
