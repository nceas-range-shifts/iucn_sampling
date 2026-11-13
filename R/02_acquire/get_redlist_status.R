### =============================================================================
### R/02_acquire/get_redlist_status.R  —  IUCN Red List status (OPTIONAL)
### =============================================================================
#
### Look up each species' IUCN Red List category. Optional: only needed if you use
### conservation status in downstream summaries.
#
### API-gated: needs IUCN_REDLIST_KEY in .Renviron. Reads species from the master
### crosswalk (name_gbif), keyed off the canonical name — not a hand-edited path.
#
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

run_redlist_status <- function(taxa = config$taxa, cfg = config) {
  key <- Sys.getenv("IUCN_REDLIST_KEY")
  if (!nzchar(key)) stop("IUCN_REDLIST_KEY not set in .Renviron.")

  crosswalk <- readr::read_csv(cfg$paths$crosswalk, show_col_types = FALSE) |>
    dplyr::filter(taxon %in% taxa, !is.na(name_gbif))

  purrr::map_dfr(unique(crosswalk$name_gbif), get_status_one, key = key) |>
    dplyr::left_join(dplyr::select(crosswalk, taxon, name_gbif),
                     by = c("scientific_name" = "name_gbif")) -> results

  out <- file.path(cfg$paths$data_processed, "iucn_redlist_status.csv")
  ensure_parent_dir(out)
  readr::write_csv(results, out)
  message(glue("[redlist] {nrow(results)} species -> {out}"))
  results
}
