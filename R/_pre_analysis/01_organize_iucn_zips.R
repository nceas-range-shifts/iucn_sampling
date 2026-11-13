### =============================================================================
### R/_pre_analysis/organize_iucn_zips.R
### =============================================================================
###
### Unzip IUCN downloads in data/raw/iucn/ and organize them into taxon folders:
###
###   data/raw/iucn/<taxon>/
###     spatial_data.gpkg
###     meta/*.csv
###
### Zips are either spatial (contain spatial_data.gpkg) or metadata (contain
### taxonomy.csv etc.). Each spatial gpkg is paired with the metadata zip whose
### taxonomy.csv `internalTaxonId` values best overlap the gpkg `id_no` values.
### Folders are provisionally named after the dominant `className` in the
### taxonomy table (suffixed _2, _3... on collision); rename as desired.
###
### Usage: source("R/_pre_analysis/organize_iucn_zips.R")

suppressPackageStartupMessages({
  library(sf)
  library(readr)
  library(dplyr)
})


### Unzip -------------------------------------------------------------------

### Unzip into a separate _unzipped folder

iucn_dir <- "data/raw/iucn"
stage_dir <- file.path(iucn_dir, "_unzipped")

zips <- list.files(iucn_dir, pattern = "\\.zip$", full.names = TRUE)
stopifnot(length(zips) > 0, nzchar(Sys.which("unzip")))

for (z in zips) {
  out <- file.path(stage_dir, tools::file_path_sans_ext(basename(z)))
  if (dir.exists(out)) next
  message("Unzipping ", basename(z))
  dir.create(out, recursive = TRUE)
  status <- system2("unzip", c("-q", shQuote(z), "-d", shQuote(out)))
  if (status != 0) {
    unlink(out, recursive = TRUE)
    stop("Failed to unzip ", z)
  }
}

### Identify whether each unzipped folder is spatial or meta based
### on typical file names
stage_dirs <- list.dirs(stage_dir, recursive = FALSE)
is_spatial <- file.exists(file.path(stage_dirs, "spatial_data.gpkg"))
is_meta <- file.exists(file.path(stage_dirs, "taxonomy.csv"))
spatial_dirs <- stage_dirs[is_spatial]
meta_dirs <- stage_dirs[is_meta]
message(length(spatial_dirs), " spatial and ", length(meta_dirs), " metadata zips")

### Species IDs per zip -----------------------------------------------------

### For each spatial file, extract the id_no field
spatial_ids <- lapply(spatial_dirs, function(d) {
  ### d <- spatial_dirs[1]
  gpkg <- file.path(d, "spatial_data.gpkg")
  layer <- st_layers(gpkg)$name[1]
  ### Attribute-only query avoids reading geometries
  q <- sprintf("SELECT DISTINCT id_no FROM \"%s\"", layer)
  result <- st_read(gpkg, query = q, quiet = TRUE)$id_no
})

meta_tax <- lapply(meta_dirs, function(d) {
  read_csv(file.path(d, "taxonomy.csv"), col_types = cols(.default = "c"),
           progress = FALSE)
})

### Match spatial to metadata -----------------------------------------------

### Fraction of each gpkg's species found in each metadata taxonomy
### (rows = spatial zips, columns = metadata zips)
overlap <- matrix(NA_real_, nrow = length(spatial_dirs), ncol = length(meta_dirs),
                  dimnames = list(basename(spatial_dirs), basename(meta_dirs)))

for (s in seq_along(spatial_dirs)) {
  for (m in seq_along(meta_dirs)) {
    meta_ids <- as.numeric(meta_tax[[m]]$internalTaxonId)
    ### populate matrix cell with the percent of match
    overlap[s, m] <- mean(spatial_ids[[s]] %in% meta_ids)
  }
}

overlap_df <- as.data.frame(overlap)
overlap_df$spatial <- rownames(overlap_df)
meta_names <- 
overlap_long <- reshape(
  data = overlap_df,
  direction = 'long',
  varying = basename(meta_dirs),
  v.names = 'pct_match',
  timevar = 'meta',
  times = basename(meta_dirs),
  idvar = 'spatial'
)
rownames(overlap_long) <- NULL

### each spatial matches exactly for one meta and vice versa
best <- overlap_long[overlap_long$pct_match == 1, ]
if(nrow(best) != length(spatial_dirs) | nrow(best) != length(meta_dirs)) {
  stop('Mismatch between spatial and meta directories!')
}

### Organize ----------------------------------------------------------------

### Dominant className of each matched metadata set
best$taxon <- NA_character_
for (i in 1:nrow(best)) {
  tax <- meta_tax[[match(best$meta[i], basename(meta_dirs))]]
  nm <- names(sort(table(tax$className), decreasing = TRUE))[1]
  best$taxon[i] <- gsub("[^a-z0-9]+", "_", tolower(nm))
}
### if any duplicated class names, append a number to differentiate
best$taxon <- ave(best$taxon, best$taxon, FUN = function(x) {
  if (length(x) == 1) x else paste0(x, "_", seq_along(x))
})
print(best)

### repeat for ambiguous matches
amb <- overlap_long[overlap_long$pct_match < 1 & overlap_long$pct_match > 0, ]
amb$taxon <- NA_character_
for (i in 1:nrow(amb)) {
  tax <- meta_tax[[match(amb$meta[i], basename(meta_dirs))]]
  nm <- names(sort(table(tax$className), decreasing = TRUE))[1]
  amb$taxon[i] <- gsub("[^a-z0-9]+", "_", tolower(nm))
}
### if any duplicated class names, append a number to differentiate
amb$taxon <- ave(amb$taxon, amb$taxon, FUN = function(x) {
  if (length(x) == 1) x else paste0(x, "_", seq_along(x))
})
print(amb)

### Now using best, organize all the files

for (i in seq_len(nrow(best))) {
  spatial_dir <- file.path(stage_dir, best$spatial[i])
  meta_dir    <- file.path(stage_dir, best$meta[i])
  target      <- file.path(iucn_dir, best$taxon[i])
  if (dir.exists(target)) {
    warning(target, " already exists; skipping")
    next()
  }
  dir.create(file.path(target, "meta"), recursive = TRUE)
  file.rename(file.path(spatial_dir, "spatial_data.gpkg"),
              file.path(target, "spatial_data.gpkg"))
  csvs <- list.files(meta_dir, pattern = "\\.csv$", full.names = TRUE)
  file.rename(csvs, file.path(target, "meta", basename(csvs)))
  message(best$spatial[i], " + ", best$meta[i], " -> ", target)
}

### Clean up staging dir only if everything has been moved out
leftover <- list.files(stage_dir, pattern = "\\.(gpkg|csv)$", recursive = TRUE)
if (length(leftover) == 0) {
  unlink(stage_dir, recursive = TRUE)
} else {
  message("Unmoved files remain in ", stage_dir, "; not deleting it.")
}
