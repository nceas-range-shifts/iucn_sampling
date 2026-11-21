### =============================================================================
### R/01_taxonomy/01_build_species_crosswalk.R  —  STAGE 0 driver
### =============================================================================
###
### Turns each taxon's raw input species list into ONE authoritative crosswalk
### that every downstream stage joins on. Runs FIRST in the pipeline.
###
### What it adds beyond resolve_species_names():
###   - range_polygon_name / range_polygon_path / range_polygon_layer for polygon-range taxa.
###     Range polygons carry the vendor's own species name (the polygon_name_
###     field, e.g. "SName"), which may not match the GBIF canonical name. We
###     resolve those polygon names through the GBIF backbone and join them to
###     the crosswalk on the shared gbif_usage_key -- the same key-bridge the
###     resolver already uses for eBird. This is the join key that lets range
###     cells (keyed by polygon name) line up with occurrence/richness data
###     (keyed by GBIF canonical name).
###
###   Polygon sources may be multi-species (e.g. IUCN geopackages holding many
###   species as rows of one layer) or single-species (one file per species), so
###   the link records the file path, the layer within it, and the species name.
###   Rows of a layer are later selected by matching the name field.
###
###   - aphia_id (+ aphia_id_matched, worms_name, worms_status, worms_n_matches)
###     for taxa with an obis_dataset_path: the accepted WoRMS AphiaID, OBIS's
###     taxon key, matched exactly by name. NA for non-marine taxa or no match.
###
### Output: data/processed/species_crosswalk.csv (one row per species, all taxa).
###
### Requires config/config.R + config/taxa_config.R sourced first.
###
### Usage:
###   source("config/config.R"); source("config/taxa_config.R")
###   source("R/00_utils.R")
###   source("R/01_taxonomy/01_build_species_crosswalk.R")
###   crosswalk <- build_species_crosswalk()

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(glue)
  library(janitor)
  library(sf)
  library(here)
})

source(here::here("R", "01_taxonomy", "01a_resolve_species_names.R"))

### Load one taxon's raw species list from its configured CSV or Google Sheet.
### When reading a shared CSV that carries a `taxon` column, filter to this taxon.
load_species_list <- function(taxon, tcfg) {
  src <- tcfg$species_source
  has_sheet <- !is.null(src$sheet_url) && !is.na(src$sheet_url) && nzchar(src$sheet_url)

  if (has_sheet) {
    googlesheets4::gs4_deauth()
    df <- googlesheets4::read_sheet(ss = src$sheet_url, sheet = src$sheet_name)
  } else if (!is.null(src$csv_path) && file.exists(src$csv_path)) {
    df <- readr::read_csv(src$csv_path, show_col_types = FALSE)
  } else {
    stop(glue("No species source found for taxon '{taxon}' (csv_path missing)."))
  }

  ### Optional filter when a shared list holds multiple taxa.
  if (!is.null(src$taxon_filter)) {
    df_clean <- janitor::clean_names(df)
    if ("taxon" %in% names(df_clean)) {
      df <- df[tolower(df_clean$taxon) == tolower(src$taxon_filter), , drop = FALSE]
    }
  }
  df
}

### Count features per species name in one layer of a vector file, attributes
### only. Tries an OGR SQL query first (avoids loading geometry for big layers),
### falling back to a full read. Returns tibble(name, n_features) or NULL.
.layer_polygon_names <- function(path, layer, name_field) {
  q <- glue('SELECT "{name_field}" AS nm, COUNT(*) AS n FROM "{layer}" GROUP BY "{name_field}"')
  counts <- tryCatch({
    res <- suppressWarnings(sf::st_read(path, query = q, quiet = TRUE))
    tibble::tibble(name = as.character(res$nm), n_features = as.integer(res$n))
  }, error = function(e) NULL)

  if (is.null(counts)) {
    atts <- sf::st_read(path, layer = layer, quiet = TRUE) |> sf::st_drop_geometry()
    if (!name_field %in% names(atts)) {
      warning(glue("'{basename(path)}' layer '{layer}' has no field '{name_field}'; skipped"))
      return(NULL)
    }
    counts <- dplyr::count(tibble::tibble(name = as.character(atts[[name_field]])),
                           name, name = "n_features")
  }
  dplyr::filter(counts, !is.na(name), name != "")
}

### Read the species-name attribute out of every vector file in a directory.
### Handles multi-species files (many species per layer, e.g. IUCN geopackages)
### and single-species files alike: every layer of every .gpkg, plus each .shp
### (whose layer is its file stem). Returns one row per species name per layer:
### range_polygon_name, range_polygon_path, range_polygon_layer, range_polygon_n_features.
harvest_polygon_names <- function(polygon_dir, name_field) {
  empty <- tibble::tibble(range_polygon_name = character(),
                          range_polygon_path = character(),
                          range_polygon_layer = character(),
                          range_polygon_n_features = integer())
  if (is.null(polygon_dir) || !dir.exists(polygon_dir)) return(empty)

  files <- list.files(polygon_dir, pattern = "\\.(shp|gpkg)$",
                      full.names = TRUE, recursive = TRUE)

  ### Paths are stored in the crosswalk relative to the project root so the file
  ### is portable across machines; consumers resolve them with here::here().
  root_prefix <- paste0(normalizePath(here::here(), mustWork = TRUE), "/")
  to_project_relative <- function(p) {
    p_abs <- normalizePath(p, mustWork = TRUE)
    if (!startsWith(p_abs, root_prefix)) {
      warning(glue("'{p}' is outside the project root; storing absolute path"))
      return(p_abs)
    }
    substring(p_abs, nchar(root_prefix) + 1)
  }

  purrr::map_dfr(files, function(f) {
    tryCatch({
      layers <- if (grepl("\\.gpkg$", f)) sf::st_layers(f)$name
                else tools::file_path_sans_ext(basename(f))
      purrr::map_dfr(layers, function(lyr) {
        cnt <- .layer_polygon_names(f, lyr, name_field)
        if (is.null(cnt) || nrow(cnt) == 0) return(NULL)
        tibble::tibble(range_polygon_name = cnt$name, range_polygon_path = to_project_relative(f),
                       range_polygon_layer = lyr, range_polygon_n_features = cnt$n_features)
      })
    }, error = function(e) {
      warning(glue("could not read '{basename(f)}': {conditionMessage(e)}"))
      NULL
    })
  }) |>
    dplyr::bind_rows(empty)
}

### Resolve polygon names to GBIF keys and attach them to the crosswalk on
### gbif_usage_key. Rows in the crosswalk with no matching polygon get NA.
attach_polygon_link <- function(crosswalk, polygon_names) {
  if (nrow(polygon_names) == 0) {
    return(dplyr::mutate(crosswalk,
                         range_polygon_name = NA_character_,
                         range_polygon_path = NA_character_,
                         range_polygon_layer = NA_character_,
                         range_polygon_n_features = NA_integer_))
  }
  keys <- rgbif::name_backbone_checklist(polygon_names$range_polygon_name) |>
    dplyr::select(range_polygon_name = verbatim_name,
                  gbif_usage_key = speciesKey)

  ### One crosswalk row per species. If a species has polygons in several
  ### layers/files (or under several vendor names), the sources are collapsed into
  ### parallel "|"-delimited name/path/layer strings that the range step expands
  ### and combines. n_features totals polygon features across all sources.
  link <- polygon_names |>
    dplyr::left_join(keys, by = "range_polygon_name") |>
    dplyr::filter(!is.na(gbif_usage_key)) |>
    dplyr::summarise(
      range_polygon_name = paste(range_polygon_name, collapse = "|"),
      range_polygon_path = paste(range_polygon_path, collapse = "|"),
      range_polygon_layer = paste(range_polygon_layer, collapse = "|"),
      range_polygon_n_features = sum(range_polygon_n_features),
      .by = gbif_usage_key
    )

  multi <- dplyr::filter(link, range_polygon_n_features > 1)
  if (nrow(multi) > 0) {
    warning(glue("{nrow(multi)} species have >1 polygon feature (e.g. differing ",
                 "presence/seasonality); the range step will combine them: ",
                 "{paste(utils::head(multi$range_polygon_name, 5), collapse = ', ')}",
                 "{if (nrow(multi) > 5) ', ...' else ''}"), call. = FALSE)
  }

  crosswalk |>
    dplyr::left_join(link, by = "gbif_usage_key")
}

### Add WoRMS AphiaIDs (OBIS's taxon key) to a crosswalk. Names are matched
### exactly (no fuzzy matching) on name_gbif, falling back to name_accepted.
### aphia_id is the ACCEPTED (valid) AphiaID, so synonyms resolve to the same id
### OBIS stores; aphia_id_matched keeps the id of the name as matched, and
### worms_status records the match outcome ("no match" -> aphia_id is NA).
attach_worms_aphia <- function(crosswalk, taxon, chunk_size = 50,
                               progress = c("terse", "verbose", "quiet")) {
  progress <- match.arg(progress)
  query <- dplyr::coalesce(crosswalk$name_gbif, crosswalk$name_accepted)
  uniq <- unique(stats::na.omit(query))

  chunks <- split(uniq, ceiling(seq_along(uniq) / chunk_size))
  worms <- purrr::imap_dfr(chunks, function(nms, i) {
    .progress_msg(progress, "verbose",
                  glue("    [{taxon}] WoRMS chunk {i}/{length(chunks)} ({length(nms)} names)"))
    res <- tryCatch(
      worrms::wm_records_names(nms, marine_only = FALSE),
      error = function(e) {
        warning(glue("WoRMS lookup failed for {length(nms)} names: {conditionMessage(e)}"),
                call. = FALSE)
        NULL
      }
    )
    if (is.null(res)) return(NULL)
    purrr::map2_dfr(nms, res, function(nm, df) {
      if (nrow(df) == 0) return(NULL)
      ### Several records can share a name (homonyms): prefer accepted, then first.
      df |>
        dplyr::arrange(status != "accepted") |>
        dplyr::slice(1) |>
        dplyr::transmute(
          query_name = nm,
          aphia_id_matched = as.integer(AphiaID),
          aphia_id = as.integer(valid_AphiaID),
          worms_name = valid_name,
          worms_status = status,
          worms_n_matches = nrow(df)
        )
    })
  })

  out <- crosswalk |>
    dplyr::mutate(.query_name = query) |>
    dplyr::left_join(worms, by = c(".query_name" = "query_name")) |>
    dplyr::select(-.query_name) |>
    dplyr::mutate(worms_status = dplyr::coalesce(worms_status, "no match"))

  n_miss <- sum(is.na(out$aphia_id))
  .progress_msg(progress, "terse",
                glue("[crosswalk] WoRMS: {nrow(out) - n_miss}/{nrow(out)} {taxon} names matched ",
                     "to an AphiaID ({n_miss} NA)"))
  out
}

### Build one taxon's crosswalk (names resolved + polygon linkage where used).
build_taxon_crosswalk <- function(taxon, tcfg, cfg = config, force = config$force,
                                  progress = c("terse", "verbose", "quiet"),
                                  raw = NULL) {
  progress <- match.arg(progress)
  ### for new alphanumeric speciesKey taxonomic keys from Catalog of Life
  check_package_version('rgbif', '3.9')

  .progress_msg(progress, "terse", glue("[crosswalk] {taxon}: resolving species names"))
  ### `raw` is passed in when the caller already loaded the list (for the cache hash).
  if (is.null(raw)) raw <- load_species_list(taxon, tcfg)
  .progress_msg(progress, "verbose", glue("    [{taxon}] {nrow(raw)} species in input list"))

  cw <- resolve_species_names(
    spp_df = raw,
    sci_name_col = tcfg$species_source$sci_name_col,
    taxon_group = taxon,          ### stamp group so reptiles aren't left NA
    rebuild_ebird_cache = isTRUE(force)
  )

  ### Report GBIF backbone match rate (mirrors the WoRMS message).
  n_gbif_miss <- sum(is.na(cw$gbif_usage_key))
  .progress_msg(progress, "terse",
                glue("[crosswalk] GBIF: {nrow(cw) - n_gbif_miss}/{nrow(cw)} {taxon} names matched ",
                     "to a GBIF usage key ({n_gbif_miss} NA)"))

  if (is_polygon_source(tcfg$range_source)) {
    .progress_msg(progress, "verbose", glue("    [{taxon}] harvesting range polygon names"))
    polygon_names <- harvest_polygon_names(tcfg$polygon_dir, tcfg$polygon_name_field)
    .progress_msg(progress, "verbose",
                  glue("    [{taxon}] {nrow(polygon_names)} polygon name/layer records found"))
    cw <- attach_polygon_link(cw, polygon_names)
    n_poly_miss <- sum(is.na(cw$range_polygon_name))
    .progress_msg(progress, "terse",
                  glue("[crosswalk] Polygons: {nrow(cw) - n_poly_miss}/{nrow(cw)} {taxon} species ",
                       "linked to a range polygon ({n_poly_miss} NA)"))

    ### Flag species whose range polygon we could not match -> review.
    cw <- cw |>
      dplyr::mutate(
        needs_review = needs_review | is.na(range_polygon_name),
        review_reason = dplyr::if_else(
          is.na(range_polygon_name) & !override,
          stringr::str_trim(paste(review_reason, "no range polygon matched", sep = "; ")),
          review_reason
        )
      )
  } else {
    cw <- dplyr::mutate(cw,
                        range_polygon_name = NA_character_,
                        range_polygon_path = NA_character_,
                        range_polygon_layer = NA_character_,
                        range_polygon_n_features = NA_integer_)
  }

  ### Marine taxa (those with an obis_dataset_path) also get WoRMS AphiaIDs. Other
  ### taxa carry NA columns (not applicable) so the crosswalk stays rectangular.
  if (!is.null(tcfg$obis_dataset_path)) {
    cw <- attach_worms_aphia(cw, taxon, progress = progress)
  } else {
    cw <- dplyr::mutate(cw,
                        aphia_id_matched = NA_integer_, aphia_id = NA_integer_,
                        worms_name = NA_character_, worms_status = NA_character_,
                        worms_n_matches = NA_integer_)
  }

  dplyr::mutate(cw, taxon = taxon, .before = 1)
}

### Build the master crosswalk across all configured taxa and write it once.
### `force = TRUE` also refreshes the cached eBird-to-GBIF key map.
### `progress`: "terse" (default; per-taxon GBIF/polygon/WoRMS match counts),
### "verbose" (adds input sizes and WoRMS chunk progress), or "quiet".
build_species_crosswalk <- function(taxa = config$taxa,
                                    taxa_cfg = taxa_config,
                                    cfg = config,
                                    write = TRUE,
                                    force = config$force,
                                    progress = c("terse", "verbose", "quiet")) {
  progress <- match.arg(progress)
  ### test to warn if config and taxa_config don't match for all taxa!
  all_missed <- c(setdiff(taxa, names(taxa_config)), setdiff(names(taxa_config), taxa))
  if(length(all_missed) != 0) {
    warning('Unmatched taxa!',
            '\n  Taxa in config:      ', paste0(taxa, collapse = ', '),
            '\n  Taxa in taxa_config: ', paste0(names(taxa_cfg), collapse = ', '))
  }
  
  ### Load every raw species list once; they feed both the cache hash and the build.
  raw_lists <- purrr::map(rlang::set_names(taxa), function(taxon) {
    tcfg <- taxa_cfg[[taxon]]
    if (is.null(tcfg)) stop(glue("Taxon '{taxon}' not found in taxa_config."))
    load_species_list(taxon, tcfg)
  })

  ### --- Cache check ---------------------------------------------------------
  ### Settings hash covers the species list CONTENT (sheets can change in place),
  ### each taxon's config, polygon file checksums, package versions and the code of
  ### the functions that build the crosswalk. It does NOT cover upstream changes
  ### (GBIF/WoRMS/eBird backbones); use force = TRUE to re-query those.
  out <- cfg$paths$crosswalk
  mpath <- file.path(dirname(out), "species_crosswalk_manifest.json")
  polygon_files <- sort(unlist(purrr::map(taxa, function(taxon) {
    tcfg <- taxa_cfg[[taxon]]
    if (is_polygon_source(tcfg$range_source) && !is.null(tcfg$polygon_dir) &&
        dir.exists(tcfg$polygon_dir)) {
      list.files(tcfg$polygon_dir, pattern = "\\.(shp|gpkg)$",
                 full.names = TRUE, recursive = TRUE)
    }
  })))
  fns <- c("build_species_crosswalk", "build_taxon_crosswalk", "load_species_list",
           "harvest_polygon_names", ".layer_polygon_names", "attach_polygon_link",
           "attach_worms_aphia", "resolve_species_names")
  settings <- rlang::hash(list(
    taxa     = taxa,
    species  = raw_lists,
    taxa_cfg = taxa_cfg[taxa],
    polygons = unname(tools::md5sum(polygon_files)),
    packages = lapply(c("rgbif", "worrms"), function(p) as.character(utils::packageVersion(p))),
    code     = lapply(fns, function(f) deparse(get(f, envir = globalenv())))
  ))

  if (isTRUE(write) && .stage_cached(mpath, settings, out, force)) {
    .progress_msg(progress, "terse",
                  glue("[crosswalk] CACHED (inputs unchanged; force = TRUE to rebuild): {out}"))
    ### gbif_usage_key is alphanumeric (e.g. "5384V"); keep it character.
    return(invisible(readr::read_csv(
      out, show_col_types = FALSE,
      col_types = readr::cols(gbif_usage_key = readr::col_character(),
                              .default = readr::col_guess()))))
  }

  crosswalk <- purrr::map_dfr(taxa, function(taxon) {
    build_taxon_crosswalk(taxon, taxa_cfg[[taxon]], cfg, force = force,
                          progress = progress, raw = raw_lists[[taxon]])
  })

  if (isTRUE(write)) {
    ensure_parent_dir(out)
    readr::write_csv(crosswalk, out)
    ### Manifest last, so an interrupted run is never treated as cached.
    .write_stage_manifest(mpath, settings, n_species = nrow(crosswalk))
    .progress_msg(progress, "terse", glue(
      "[crosswalk] {nrow(crosswalk)} species across {length(taxa)} taxa -> {out} ",
      "({sum(crosswalk$needs_review)} flagged for review)"
    ))
  }

  return(invisible(crosswalk))
}
