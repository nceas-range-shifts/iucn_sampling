# =============================================================================
# R/04_indicators/greedy_setcover_fn.R  —  greedy weighted set-cover (functions)
# =============================================================================
#
# Pick a minimal set of hex cells that brings every undersampled species up to
# its sampling target, with a conservation-area balance adjustment. Consumes the
# undersamp_cells / undersamp_spp outputs of find_optimal_sampling_locs.
#
# Reworked from run_greedy_randomized_fn.R (renamed: the algorithm is
# deterministic greedy, not randomized). The three tuning knobs now come from
# sampling_params; the conservation flag is derived inline from the hex grid +
# conservation cells (no separate ca_h3_grid_w_cons file). The un-parameterized
# third copy (99_utilities/run_greedy_randomized.R) is dropped.
#
# Requires config + sampling_config + R/00_utils.R.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(sf)
  library(arrow)
  library(readr)
  library(glue)
  library(here)
})

# Resolve the taxa a run covers ("combined" spans all configured taxa).
resolve_greedy_taxa <- function(taxon, cfg = config) {
  if (identical(taxon, "combined")) return(cfg$taxa)
  if (!taxon %in% cfg$taxa) stop(glue("taxon must be one of {toString(cfg$taxa)}, or 'combined'."))
  taxon
}

build_greedy_output_dir <- function(year, taxon, cfg = config, taxa_resolved = NULL) {
  year <- as.character(year)
  if (identical(taxon, "combined")) {
    out_dir <- file.path(cfg$paths$output, "indicators", "species_populations",
                         "combined", "sampling_optimized",
                         paste(taxa_resolved, collapse = "_"), year)
  } else {
    out_dir <- file.path(taxon_output_dir(taxon, cfg), "sampling_optimized", year)
  }
  ensure_dir(out_dir)
}

# Region hex grid annotated with a 0/1 cons_area flag.
.hex_grid_with_cons <- function(cfg = config) {
  grid <- sf::st_read(cfg$paths$hex_grid, quiet = TRUE)
  cons_ids <- sf::st_read(cfg$paths$cons_areas_hex, quiet = TRUE)$h3_index
  dplyr::mutate(grid, cons_area = as.integer(h3_index %in% cons_ids))
}

greedy_setcover <- function(year, taxon, cfg = config, params = sampling_params) {
  year <- as.character(year)
  taxa <- resolve_greedy_taxa(taxon, cfg)
  samp <- make_sampling_config(year, cfg = cfg)

  cons_balance_weight <- params$cons_balance_weight
  batch_fraction      <- params$batch_fraction
  deficit_percentile  <- params$deficit_percentile

  output_dir <- build_greedy_output_dir(year, taxon, cfg, taxa)

  # Stack candidate cells, species deficits, and range cells across taxa.
  updated_ranges <- rbindlist(lapply(taxa, function(t) {
    dt <- as.data.table(read_parquet(samp[[t]]$undersamp_cells_out)); dt[, taxon := t]; dt
  }))
  undersamp_species <- rbindlist(lapply(taxa, function(t) {
    df <- read_csv(samp[[t]]$undersamp_spp_out, show_col_types = FALSE); df$taxon <- t; df
  }), fill = TRUE)
  range_cells_long <- rbindlist(lapply(taxa, function(t) {
    rc <- as.data.table(read_parquet(samp[[t]]$range_cells_path))
    melt(rc, id.vars = "h3_index", variable.name = "species",
         value.name = "in_range")[in_range == 1]
  }))

  CA_hex_grid <- .hex_grid_with_cons(cfg)
  cons_lookup <- unique(as.data.table(sf::st_drop_geometry(CA_hex_grid))[, .(h3_index, cons_area)])

  dup_spp <- undersamp_species[, .N, by = species][N > 1, species]
  if (length(dup_spp) > 0) {
    warning("Species in multiple taxa: ", paste(dup_spp, collapse = ", "))
  }

  # Sparse candidate edge list (cell-species pairs), with cons flag.
  edges <- as.data.table(updated_ranges)[pres_abs == 1, .(h3_index, species)]
  rm(updated_ranges)
  edges <- merge(edges, cons_lookup, by = "h3_index", all.x = TRUE)
  edges[is.na(cons_area), cons_area := 0L]

  # Per-species conservation fraction of full range (for the balance target).
  zm_long <- merge(range_cells_long, cons_lookup, by = "h3_index", all.x = TRUE)
  zm_long[is.na(cons_area), cons_area := 0L]
  cons_target <- zm_long[, .(cons_frac_range = sum(cons_area) / .N), by = species]

  spp_dt <- as.data.table(undersamp_species)[, .(species, taxon, cells_needed)]
  spp_dt <- merge(spp_dt, cons_target, by = "species", all.x = TRUE)
  spp_dt[is.na(cons_frac_range), cons_frac_range := 0]
  rm(zm_long, range_cells_long); gc()

  remaining <- copy(spp_dt); setkey(remaining, species)
  selected <- vector("list", 5000L)
  n_batches <- 0L; n_selected <- 0L
  total_cons_selected <- 0L; total_noncons_selected <- 0L
  setkey(edges, species)
  edges <- edges[species %in% remaining$species]

  while (nrow(remaining[cells_needed > 0]) > 0) {
    active_spp <- remaining[cells_needed > 0, species]
    edges <- edges[species %in% active_spp]
    cell_scores <- edges[, .(n_spp = .N, cons_area = cons_area[1]), by = h3_index]
    if (nrow(cell_scores) == 0) break

    median_target <- median(remaining[cells_needed > 0, cons_frac_range])
    current_cons_frac <- if ((total_cons_selected + total_noncons_selected) == 0) 0.5 else
      total_cons_selected / (total_cons_selected + total_noncons_selected)
    cons_gap <- median_target - current_cons_frac
    cell_scores[, adj_score := n_spp + cons_balance_weight * n_spp *
                  fifelse(cons_area == 1, cons_gap, -cons_gap)]

    deficit_p <- remaining[cells_needed > 0, quantile(cells_needed, probs = deficit_percentile)]
    batch_size <- max(1L, floor(deficit_p * batch_fraction))
    batch_size <- min(batch_size, nrow(cell_scores))

    setorder(cell_scores, -adj_score)
    batch <- cell_scores[1:batch_size]
    batch_cells <- batch$h3_index

    batch_edges <- edges[h3_index %in% batch_cells]
    spp_hits <- batch_edges[, .(hits = .N), by = species]
    remaining[spp_hits, cells_needed := cells_needed - i.hits, on = "species"]

    batch_cons <- sum(batch$cons_area == 1L)
    total_cons_selected <- total_cons_selected + batch_cons
    total_noncons_selected <- total_noncons_selected + (batch_size - batch_cons)

    n_batches <- n_batches + 1L
    selected[[n_batches]] <- data.table(
      order = seq(n_selected + 1L, n_selected + batch_size),
      h3_index = batch$h3_index, cons_area = batch$cons_area,
      n_spp_covered = batch$n_spp, batch_id = n_batches,
      spp_remaining = nrow(remaining[cells_needed > 0])
    )
    n_selected <- n_selected + batch_size

    edges <- edges[!h3_index %in% batch_cells]
    satisfied <- remaining[cells_needed <= 0, species]
    if (length(satisfied) > 0) edges <- edges[!species %in% satisfied]
  }

  selected_cells <- rbindlist(selected[seq_len(n_batches)])
  species_summary <- merge(spp_dt, remaining[, .(species, final_deficit = cells_needed)],
                           by = "species")
  species_summary[, cells_selected := cells_needed - final_deficit]
  species_summary[, target_met := final_deficit <= 0]
  taxon_summary <- species_summary[, .(n_species = .N, n_at_target = sum(target_met),
    pct_at_target = round(100 * sum(target_met) / .N, 1)), by = taxon]

  fwrite(selected_cells, file.path(output_dir, "greedy_selected_cells.csv"))
  fwrite(species_summary, file.path(output_dir, "greedy_species_summary.csv"))
  selected_sf <- CA_hex_grid |>
    dplyr::filter(h3_index %in% selected_cells$h3_index) |>
    dplyr::left_join(as.data.frame(selected_cells), by = "h3_index")
  sf::st_write(selected_sf, file.path(output_dir, "greedy_selected_cells.gpkg"),
               delete_dsn = TRUE, quiet = TRUE)

  message(glue("[greedy {taxon} {year}] {n_selected} cells, ",
               "{sum(species_summary$target_met)}/{nrow(species_summary)} species at target ",
               "-> {output_dir}"))
  list(selected_cells = selected_cells, species_summary = species_summary,
       taxon_summary = taxon_summary, output_dir = output_dir)
}
