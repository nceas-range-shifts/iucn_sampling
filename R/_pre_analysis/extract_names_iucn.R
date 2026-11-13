### =============================================================================
### R/_pre_analysis/extract_names_iucn.R
### =============================================================================
###
### Extract scientific names (binomial nomenclature) from IUCN spatial data.
###
### - Finds all geopackages/shapefiles/geodatabases in data/raw/iucn_2026/
### - Reads the spatial data attributes (no geometry needed)
### - Extracts the scientific name field (typically `sci_name`) and metadata
### - Adds a `taxon` column based on the folder name (e.g., "mammals_marine")
### - Returns a deduplicated species list per taxon with IUCN metadata
###
### Output: A tibble with columns:
###   - taxon: folder name (e.g., "mammals_marine", "mammals_mar_terr")
###   - scientific_name: binomial nomenclature from `sci_name` attribute
###   - presence: presence status code (1=Extant, 3=Probably Extant, etc.)
###   - origin: origin code (1=Native, 2=Reintroduced, 3=Introduced, etc.)
###   - seasonal: seasonal status (1=Resident, 2=Breeding, 4=Passage, etc.)
###   - legend: human-readable status (e.g., "Extant (resident)")
###   - yrcompiled: year the range was compiled
###   - assessment_id: IUCN assessment ID
###   - id_no: IUCN species ID number
###
### Usage:
###   species_from_iucn <- extract_iucn_species()
###   head(species_from_iucn)
###   write_csv(species_from_iucn, "data/raw/iucn_species_list.csv")

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(readr)
  library(purrr)
  library(glue)
  library(here)
})

### Extract scientific names and metadata from one spatial file (gpkg, shp, gdb, etc.)
### Returns a tibble with columns: taxon, scientific_name, presence, origin, seasonal, 
###                               legend, yrcompiled, assessment_id, id_no
extract_from_spatial_file <- function(file_path, taxon_name, sci_name_field = "sci_name") {
  tryCatch({
    ### Read spatial data, drop geometry to speed up
    data <- sf::st_read(file_path, quiet = TRUE) |>
      sf::st_drop_geometry()
    
    ### Check if the science name field exists
    if (!sci_name_field %in% names(data)) {
      warning(glue("File '{basename(file_path)}' has no '{sci_name_field}' field. ",
                   "Available fields: {paste(names(data), collapse=', ')}"))
      return(tibble::tibble(
        taxon = character(), scientific_name = character(),
        presence = integer(), origin = integer(), seasonal = integer(),
        legend = character(), yrcompiled = integer(),
        assessment_id = numeric(), id_no = numeric()
      ))
    }
    
    ### Select relevant columns (use only those that exist in the file)
    cols_to_keep <- c(sci_name_field, "presence", "origin", "seasonal", 
                      "legend", "yrcompiled", "assessment_id", "id_no")
    cols_available <- cols_to_keep[cols_to_keep %in% names(data)]
    
    ### Extract and clean scientific names with metadata
    species <- data |>
      dplyr::select(all_of(cols_available)) |>
      dplyr::rename(scientific_name = all_of(sci_name_field)) |>
      ### Remove empty/NA scientific names
      dplyr::filter(!is.na(scientific_name) & scientific_name != "") |>
      ### Deduplicate by scientific name, keeping first occurrence
      dplyr::distinct(scientific_name, .keep_all = TRUE) |>
      ### Add taxon column
      dplyr::mutate(taxon = taxon_name, .before = 1)
    
    species
  }, error = function(e) {
    warning(glue("Error reading '{basename(file_path)}': {conditionMessage(e)}"))
    tibble::tibble(
      taxon = character(), scientific_name = character(),
      presence = integer(), origin = integer(), seasonal = integer(),
      legend = character(), yrcompiled = integer(),
      assessment_id = numeric(), id_no = numeric()
    )
  })
}

### Main function: scan iucn_2026 folder for all spatial data and extract species
extract_iucn_species <- function(iucn_root = here("data", "raw", "iucn_2026"),
                                  sci_name_field = "sci_name") {
  
  ### Verify root exists
  if (!dir.exists(iucn_root)) {
    stop(glue("IUCN directory not found: {iucn_root}"))
  }
  
  ### List all immediate subdirectories (taxa folders)
  taxa_dirs <- list.dirs(iucn_root, recursive = FALSE, full.names = FALSE)
  
  if (length(taxa_dirs) == 0) {
    warning(glue("No subdirectories found in {iucn_root}"))
    return(tibble::tibble(taxon = character(), scientific_name = character()))
  }
  
  message(glue("Found {length(taxa_dirs)} taxon folder(s): {paste(taxa_dirs, collapse=', ')}"))
  
  ### For each taxon folder, find and read spatial files
  all_species <- purrr::map_dfr(taxa_dirs, function(taxon_name) {
    taxon_path <- file.path(iucn_root, taxon_name)
    
    message(glue("  Processing '{taxon_name}'..."))
    
    ### Find spatial files (geopackages first, then shapefiles, then geodatabases)
    gpkg_files <- list.files(taxon_path, pattern = "\\.gpkg$", 
                             full.names = TRUE, recursive = FALSE)
    shp_files  <- list.files(taxon_path, pattern = "\\.shp$", 
                             full.names = TRUE, recursive = FALSE)
    gdb_files  <- list.files(taxon_path, pattern = "\\.gdb$", 
                             include.dirs = TRUE,
                             full.names = TRUE, recursive = FALSE)
    
    spatial_files <- c(gpkg_files, shp_files, gdb_files)
    
    if (length(spatial_files) == 0) {
      warning(glue("  No spatial files found in {taxon_path}"))
      return(tibble::tibble(taxon = character(), scientific_name = character()))
    }
    
    message(glue("    Found {length(spatial_files)} spatial file(s)"))
    
    ### Read each spatial file and extract species
    taxon_species <- purrr::map_dfr(spatial_files, function(f) {
      message(glue("      Reading {basename(f)}..."))
      extract_from_spatial_file(f, taxon_name, sci_name_field)
    })
    
    ### Deduplicate within this taxon (keep all columns)
    if (nrow(taxon_species) > 0) {
      taxon_species <- taxon_species |>
        dplyr::distinct(scientific_name, .keep_all = TRUE)
      message(glue("    → {nrow(taxon_species)} unique species"))
    } else {
      message(glue("    → No species extracted"))
    }
    
    taxon_species
  })
  
  ### Overall summary
  if (nrow(all_species) > 0) {
    by_taxon <- all_species |> dplyr::count(taxon)
    message(glue("\nSummary:"))
    for (i in seq_len(nrow(by_taxon))) {
      message(glue("  {by_taxon$taxon[i]}: {by_taxon$n[i]} species"))
    }
    message(glue("Total: {nrow(all_species)} species across {dplyr::n_distinct(all_species$taxon)} taxa"))
  }
  
  all_species
}

### Run extraction
species_from_iucn <- extract_iucn_species()

### Write to CSV
output_path <- here("data", "raw", "iucn_species_list.csv")
readr::write_csv(species_from_iucn, output_path)
message(glue("\nSpecies list written to: {output_path}"))

