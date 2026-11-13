### =============================================================================
### R/02_acquire/get_redlist_status.R  —  IUCN Red List status (OPTIONAL)
### =============================================================================
###
### Look up each species' IUCN Red List category. Optional: only needed if you use
### conservation status in downstream summaries.
###
### API-gated: needs config$credentials$iucn_redlist_key (default: IUCN_REDLIST_KEY
### from .Renviron). Reads species from the master
### crosswalk (name_gbif), keyed off the canonical name — not a hand-edited path.
###
### Requires config + R/00_utils.R + a built crosswalk.  Usage: run_redlist_status()

suppressPackageStartupMessages({
  library(rredlist)
  library(dplyr)
  library(purrr)
  library(readr)
  library(glue)
  library(here)
})

### One species -> one-row status tibble (never stops the run).
get_status_one <- function(sci_name, key, pause = 2) {
  res <- tryCatch({
    parts <- strsplit(trimws(sci_name), "\\s+")[[1]]
    if (length(parts) < 2) {
      tibble::tibble(scientific_name = sci_name, category = NA_character_,
                     note = "not a binomial")
    } else {
      r <- rredlist::rl_species_latest(genus = parts[1], species = parts[2],
                                       key = key)
      code <- r$red_list_category$code
      tibble::tibble(scientific_name = sci_name,
                     category = if (length(code) == 0) "Not Found" else code,
                     note = if (length(code) == 0) "no record" else "ok")
    }
  }, error = function(e) tibble::tibble(scientific_name = sci_name,
                                        category = "Error", note = conditionMessage(e)))
  Sys.sleep(pause)
  res
}

### Manifest next to the output records the request hash and the Red List
### version the statuses came from. Cached results are reused until either the
### species request or the Red List version (e.g. "2026-1") changes.
redlist_manifest_path <- function(cfg) {
  file.path(cfg$paths$data_processed, "iucn_redlist_status_manifest.json")
}

run_redlist_status <- function(taxa = config$taxa, cfg = config, force = FALSE,
                               pause = 2) {
  key <- cfg$credentials$iucn_redlist_key
  if (is.null(key) || !nzchar(key)) {
    stop("IUCN Red List key not set: fill config$credentials$iucn_redlist_key ",
         "(default reads IUCN_REDLIST_KEY from .Renviron).")
  }

  crosswalk <- readr::read_csv(cfg$paths$crosswalk, show_col_types = FALSE) |>
    dplyr::filter(taxon %in% taxa, !is.na(name_gbif))

  out <- file.path(cfg$paths$data_processed, "iucn_redlist_status.csv")
  mpath <- redlist_manifest_path(cfg)
  req_hash <- rlang::hash(list(
    taxa = sort(unique(taxa)),
    species = sort(unique(crosswalk$name_gbif))
  ))
  rl_version <- rredlist::rl_version(key = key)

  if (!isTRUE(force) && file.exists(mpath) && file.exists(out)) {
    m <- tryCatch(jsonlite::read_json(mpath, simplifyVector = TRUE),
                  error = function(e) NULL)
    if (!is.null(m) && identical(m$request_hash, req_hash) &&
        identical(m$redlist_version, rl_version)) {
      message(glue("[redlist] using CACHED statuses from {out} (Red List {rl_version}, ",
                   "retrieved {m$retrieved_utc} UTC; no API queries made)."))
      return(readr::read_csv(out, show_col_types = FALSE))
    }
    message(glue("[redlist] refreshing (cached: {m$redlist_version %||% 'none'}, ",
                 "current: {rl_version}; request ",
                 "{if (identical(m$request_hash, req_hash)) 'unchanged' else 'changed'})."))
  }

  results <- purrr::map_dfr(unique(crosswalk$name_gbif), get_status_one,
                            key = key, pause = pause) |>
    dplyr::left_join(dplyr::select(crosswalk, taxon, name_gbif),
                     by = c("scientific_name" = "name_gbif"))

  ensure_parent_dir(out)
  readr::write_csv(results, out)
  jsonlite::write_json(
    list(request_hash = req_hash, redlist_version = rl_version,
         file = basename(out), n_species = nrow(results),
         retrieved_utc = format(Sys.time(), tz = "UTC", "%Y-%m-%d %H:%M:%S")),
    mpath, auto_unbox = TRUE, pretty = TRUE
  )
  message(glue("[redlist] {nrow(results)} species (Red List {rl_version}) -> {out}"))
  results
}
