# =============================================================================
# R/01_taxonomy/resolve_species_names.R  —  the taxonomy engine
# =============================================================================
#
# One function to turn a raw input species list into a clean, non-destructive
# crosswalk that every downstream stage joins on.
#
# Engine (in order):
#   1. gna_parse    - flag malformed / non-binomial names before we trust them.
#   2. GBIF backbone (name_backbone_checklist) - AUTHORITATIVE. Resolves synonyms
#      itself: for a synonym `status == "SYNONYM"` but `species` / `speciesKey`
#      point at the ACCEPTED taxon, which is what we want for occ_download().
#   3. gna_verifier - cross-source sanity check on the non-bird names, to surface
#      disagreements for manual review.
#   4. eBird        - attach the current speciesCode by joining on the GBIF
#      ACCEPTED key (not the verbatim name). Rescues renamed birds (Chen ->
#      Anser, Picoides -> Dryobates, Anas -> Spatula/Mareca, ...).
#
# Output: ONE tidy tibble, one row per submitted name, birds and non-birds
# together, with a stable schema (see the return docs on resolve_species_names).
#
# The shapefile linkage (range_shapefile_name) is added on TOP of this output by
# build_species_crosswalk.R, reusing the same GBIF-key bridge.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(janitor)
  library(taxize)
  library(rgbif)
  library(rebird)
})

# Map GBIF class -> the taxon_group labels used across this project.
# NOTE: the GBIF backbone reports reptiles at the "class" of Squamata /
# Testudines / Crocodylia rather than "Reptilia", so those are mapped here too.
.class_to_group <- c(
  AVES = "bird",
  MAMMALIA = "mammal",
  REPTILIA = "reptile",
  SQUAMATA = "reptile",
  TESTUDINES = "reptile",
  CROCODYLIA = "reptile",
  AMPHIBIA = "amphibian"
)

# Ranks we accept from gna_parse as "a real species-or-finer name".
.accepted_parse_ranks <- c("sp.", "var.", "ssp.", "subsp.", "f.")

# Default cache/override locations (portable, under data/processed).
.default_ebird_cache <- here::here("data", "processed", "ebird_gbif_key_map.csv")
.default_overrides   <- here::here("data", "processed", "name_overrides.csv")

# ---------------------------------------------------------------------------
# ebird_gbif_map(): eBird current name/code keyed to GBIF usageKey
# ---------------------------------------------------------------------------
# eBird tracks a more current taxonomy than the GBIF backbone (e.g. Ardenna,
# Hydrobates, Dryobates), so a name-to-name join between the two misses renamed
# birds. Instead we resolve every eBird species name through the GBIF backbone
# once and key it to the GBIF ACCEPTED speciesKey. Any GBIF-resolved name
# (Puffinus griseus, key 5739277) then joins to eBird's current name (Ardenna
# grisea, same key) on that shared key.
#
# Resolving the full eBird taxonomy (~11k species) takes ~40s, so the result is
# cached. Pass rebuild = TRUE to refresh, or cache_path = NULL to skip caching.
ebird_gbif_map <- function(cache_path = .default_ebird_cache, rebuild = FALSE) {
  if (!is.null(cache_path) && !isTRUE(rebuild) && file.exists(cache_path)) {
    return(readr::read_csv(cache_path, show_col_types = FALSE))
  }

  eb <- rebird::ebirdtaxonomy() |>
    dplyr::filter(category == "species") |>
    dplyr::select(
      ebird_name = sciName,
      common_name = comName,
      ebird_code = speciesCode
    )

  keys <- rgbif::name_backbone_checklist(eb$ebird_name) |>
    dplyr::transmute(
      ebird_name = verbatim_name,
      gbif_usage_key = speciesKey
    )

  map <- eb |>
    dplyr::left_join(keys, by = "ebird_name") |>
    dplyr::filter(!is.na(gbif_usage_key)) |>
    dplyr::distinct(gbif_usage_key, .keep_all = TRUE)

  if (!is.null(cache_path)) {
    dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(map, cache_path)
  }

  map
}

# ---------------------------------------------------------------------------
# apply_name_overrides(): fold in hand-curated fixes
# ---------------------------------------------------------------------------
# Some names legitimately can't be auto-resolved (usually the source list uses a
# newer taxonomy than the GBIF backbone). Rather than hand-edit the generated
# crosswalk (a rebuild would clobber it), record fixes ONCE in a version-
# controlled CSV and re-apply them every run.
#
# The overrides CSV must have a `name_input` column (the join key) plus any of:
#   name_accepted, name_gbif, gbif_usage_key, ebird_code, ebird_name, common_name
# Only the columns you fill are overridden. Any listed row is treated as
# human-reviewed: needs_review is set FALSE and `override` is set TRUE.
apply_name_overrides <- function(crosswalk, overrides_path = .default_overrides) {
  finish <- function(df) {
    if (!"override" %in% names(df)) df$override <- FALSE
    dplyr::relocate(df, override, .after = name_changed)
  }

  if (is.null(overrides_path) || !file.exists(overrides_path)) {
    return(finish(dplyr::mutate(crosswalk, override = FALSE)))
  }

  ov <- readr::read_csv(overrides_path, show_col_types = FALSE)
  if (!("name_input" %in% names(ov)) || nrow(ov) == 0) {
    return(finish(dplyr::mutate(crosswalk, override = FALSE)))
  }

  overridable <- intersect(
    c("name_accepted", "name_gbif", "gbif_usage_key",
      "ebird_code", "ebird_name", "common_name"),
    names(ov)
  )

  ov <- ov |>
    dplyr::mutate(name_input = stringr::str_trim(name_input)) |>
    dplyr::distinct(name_input, .keep_all = TRUE)
  if ("gbif_usage_key" %in% overridable) {
    ov$gbif_usage_key <- as.numeric(ov$gbif_usage_key)
  }

  crosswalk <- crosswalk |>
    dplyr::left_join(
      ov |> dplyr::select(name_input, dplyr::all_of(overridable)),
      by = "name_input",
      suffix = c("", "__ov")
    )

  for (col in overridable) {
    ovcol <- paste0(col, "__ov")
    crosswalk[[col]] <- dplyr::coalesce(crosswalk[[ovcol]], crosswalk[[col]])
    crosswalk[[ovcol]] <- NULL
  }

  crosswalk <- crosswalk |>
    dplyr::mutate(
      override = name_input %in% ov$name_input,
      name_changed = !is.na(name_accepted) & name_accepted != name_input,
      needs_review = dplyr::if_else(override, FALSE, needs_review),
      review_reason = dplyr::if_else(override, "manual override", review_reason)
    )

  finish(crosswalk)
}

# ---------------------------------------------------------------------------
# resolve_species_names()
# ---------------------------------------------------------------------------
# Args:
#   spp_df        data frame containing a scientific-name column.
#   sci_name_col  name of that column (cleaned with janitor internally, so
#                 "Scientific Name" and "scientific_name" both work).
#   taxon_group   optional label ("bird"/"mammal"/"reptile"/...) stamped on every
#                 row. Pass it when you already know the group; leave NULL to
#                 derive from GBIF class. Recommended for reptiles (GBIF class
#                 Squamata/Testudines rather than "Reptilia").
#   verify        run gna_verifier cross-check on non-birds (TRUE).
#   batch_size    gna_verifier batch size (the API dislikes huge calls).
#
# Returns a tibble leading with the name columns:
#   name_input, name_accepted, name_gbif, gbif_usage_key, ebird_code, ebird_name,
#   common_name, class, taxon_group, gbif_status, match_type, verifier_status,
#   verifier_source, name_changed, needs_review, review_reason
resolve_species_names <- function(
  spp_df,
  sci_name_col = "scientific_name",
  taxon_group = NULL,
  verify = TRUE,
  batch_size = 50,
  ebird_cache = .default_ebird_cache,
  rebuild_ebird_cache = FALSE,
  overrides = .default_overrides
) {
  spp_df <- janitor::clean_names(spp_df)
  sci_name_col <- janitor::make_clean_names(sci_name_col)
  if (!sci_name_col %in% names(spp_df)) {
    stop(glue::glue("Column '{sci_name_col}' not found in input data."))
  }

  names_in <- spp_df[[sci_name_col]] |>
    as.character() |>
    stringr::str_trim() |>
    (\(x) x[!is.na(x) & x != ""])() |>
    unique()

  if (length(names_in) == 0) {
    stop("No non-empty scientific names found in input.")
  }

  # 1. Parse: flag malformed / non-binomial names ------------------------------
  parsed <- taxize::gna_parse(names_in)
  parse_ok <- parsed |>
    dplyr::transmute(
      sci_name_input = verbatim,
      parse_ok = (rank %in% .accepted_parse_ranks) & cardinality >= 2
    ) |>
    dplyr::distinct(sci_name_input, .keep_all = TRUE)

  # 2. GBIF backbone: authoritative accepted name + key ------------------------
  gbif <- rgbif::name_backbone_checklist(names_in) |>
    dplyr::transmute(
      sci_name_input = verbatim_name,
      sci_name_gbif = species,             # ACCEPTED binomial
      gbif_scientific_name = scientificName,
      # speciesKey ONLY -- never fall back to usageKey. A HIGHERRANK match (genus
      # only) has no speciesKey; falling back to the genus usageKey would let the
      # eBird key-join attach an arbitrary congener. NA keeps such rows honest.
      gbif_usage_key = speciesKey,
      gbif_status = status,
      match_type = matchType,
      gbif_rank = rank,
      class = class,
      confidence = confidence
    ) |>
    dplyr::mutate(
      taxon_group = if (!is.null(taxon_group)) {
        taxon_group
      } else {
        unname(.class_to_group[toupper(class)])
      }
    ) |>
    dplyr::distinct(sci_name_input, .keep_all = TRUE)

  is_bird <- gbif$class == "Aves" & !is.na(gbif$class)

  # 3. gna_verifier cross-check on NON-BIRDS (batched) -------------------------
  verifier <- tibble::tibble(
    sci_name_input = character(),
    verifier_status = character(),
    verifier_match = character(),
    verifier_source = character()
  )
  nonbird_names <- gbif$sci_name_input[!is_bird]
  if (isTRUE(verify) && length(nonbird_names) > 0) {
    batches <- split(nonbird_names, ceiling(seq_along(nonbird_names) / batch_size))
    results_list <- vector("list", length(batches))
    for (i in seq_along(batches)) {
      tryCatch(
        {
          results_list[[i]] <- taxize::gna_verifier(batches[[i]])
        },
        error = function(e) {
          warning(paste("gna_verifier batch", i, "failed:", conditionMessage(e)))
        }
      )
    }
    verifier_raw <- dplyr::bind_rows(results_list)
    if (nrow(verifier_raw) > 0) {
      verifier <- verifier_raw |>
        dplyr::transmute(
          sci_name_input = submittedName,
          verifier_status = taxonomicStatus,
          verifier_match = matchType,
          verifier_source = dataSourceTitleShort
        ) |>
        dplyr::distinct(sci_name_input, .keep_all = TRUE)
    }
  }

  # 4. eBird: attach the CURRENT name/code by joining on the GBIF key ----------
  is_bird <- gbif$class == "Aves" & !is.na(gbif$class)
  if (any(is_bird)) {
    emap <- ebird_gbif_map(cache_path = ebird_cache, rebuild = rebuild_ebird_cache) |>
      dplyr::mutate(gbif_usage_key = as.numeric(gbif_usage_key))
  } else {
    emap <- tibble::tibble(
      gbif_usage_key = numeric(),
      ebird_name = character(),
      common_name = character(),
      ebird_code = character()
    )
  }

  # 5. Assemble the unified crosswalk ------------------------------------------
  crosswalk <- gbif |>
    dplyr::mutate(gbif_usage_key = as.numeric(gbif_usage_key)) |>
    dplyr::left_join(parse_ok, by = "sci_name_input") |>
    dplyr::left_join(verifier, by = "sci_name_input") |>
    dplyr::left_join(emap, by = "gbif_usage_key") |>
    dplyr::mutate(
      is_bird = class == "Aves" & !is.na(class),
      ebird_name = dplyr::if_else(is_bird, ebird_name, NA_character_),
      ebird_code = dplyr::if_else(is_bird, ebird_code, NA_character_),
      common_name = dplyr::if_else(is_bird, common_name, NA_character_),
      name_accepted = dplyr::coalesce(ebird_name, sci_name_gbif),
      parse_ok = dplyr::coalesce(parse_ok, FALSE),
      name_changed = !is.na(name_accepted) &
        (name_accepted != sci_name_input | gbif_status == "SYNONYM")
    )

  # 6. Flag rows a human should check before trusting --------------------------
  crosswalk <- crosswalk |>
    dplyr::mutate(
      review_reason = purrr::pmap_chr(
        list(parse_ok, match_type, sci_name_gbif, is_bird, ebird_code, verifier_status),
        function(ok, mt, acc, bird, ecode, vstat) {
          reasons <- c(
            if (isFALSE(ok)) "unparseable name",
            if (isTRUE(mt == "NONE") || is.na(acc)) "no GBIF species match",
            if (isTRUE(mt == "FUZZY")) "fuzzy GBIF match",
            if (isTRUE(bird) && is.na(ecode)) "no eBird match (GBIF/eBird taxonomy differ)",
            if (!is.na(vstat) && !(vstat %in% c("Accepted", "Synonym", "N/A", ""))) {
              paste0("verifier: ", vstat)
            }
          )
          paste(reasons, collapse = "; ")
        }
      ),
      needs_review = review_reason != ""
    ) |>
    dplyr::select(
      name_input = sci_name_input,
      name_accepted,
      name_gbif = sci_name_gbif,
      gbif_usage_key,
      ebird_code,
      ebird_name,
      common_name,
      class,
      taxon_group,
      gbif_status,
      match_type,
      verifier_status,
      verifier_source,
      name_changed,
      needs_review,
      review_reason
    )

  # 7. Fold in any hand-curated fixes (re-applied every run, never lost) -------
  apply_name_overrides(crosswalk, overrides)
}
