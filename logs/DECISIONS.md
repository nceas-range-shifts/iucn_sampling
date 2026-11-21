# DECISIONS

## 2026-11-17 18:30 UTC — OBIS source: AWS S3 GeoParquet primary, API fallback

**Decision:** Configure OBIS downloads to use AWS S3 GeoParquet as the primary source (via DuckDB), with automatic fallback to the API.

**Rationale:**
- Both `mammals_marine` and `abalones` downloads this morning (04:24–04:48 UTC) used the API, which is slower and requires pagination through 10k-record batches.
- S3 GeoParquet query is a single fast scan instead of repeated API calls; for `mammals_marine` (1.8M records), this saves significant time.
- AWS access verified working (360M total OBIS records accessible).
- API fallback ensures robustness if S3 is unavailable or returns nothing.

**Changes:**
- Added `config$obis_use_aws = TRUE` to `config/config.R` with documentation.
- Added `config$obis_max_age_days = 180` to match GBIF re-download age threshold.
- Updated `run_pipeline.R` to pass `use_aws = config$obis_use_aws` to `run_obis_downloads()`.
- Script already supported this via `download_one_taxon_obis(use_aws = TRUE)` default; just needed config to wire it in.

**Next run:** Expect OBIS downloads to show `source: "aws_geoparquet"` in manifests instead of `"api"`.


## 2026-11-18 PST — Use h3o instead of h3jsr for polygon-to-H3 conversion

**Decision:** Production range step uses `h3o::sfc_to_cells(containment = "centroid")`.

**Evidence:** Abalones, res 4, one pass: h3jsr 76.2 s vs h3o 5.0 s (15.1x overall; median 11.9x per species, IQR 9.9-14.8x). Cell sets identical for every species benchmarked. Only one run, no variance estimate; only abalones tested.

**Rationale:** Same cells, an order of magnitude faster. `"centroid"` containment chosen because h3jsr keeps cells whose centre lies in the polygon; h3o's default `"intersect"` would change range sizes.

**Rejected / ruled out:** The prototype assumed `h3o::polygon_to_cells()`; it does not exist in h3o 0.3.0. Batching all species in one call (grouping cells by species afterwards) was discussed but not implemented; multi-feature species are still combined per species before flattening.

**Open:** Not yet verified on other polygon taxa, especially antimeridian-crossing ranges that triggered the strip-cutting fallback with h3jsr. Benchmark script has a many-to-many join inflating 43 species to 55 rows. `renv.lock` not yet re-snapshotted for h3o.

## 2026-11-18 PST — Use h3o in build_hex_grid(); polygons via st_as_sfc, not h3_to_vertexes

**Decision:** Hex-grid ids and polygons are generated with h3o.

**Evidence (global grid, level 3, 41,162 cells, 3 runs each):** h3jsr 4.93 s mean vs h3o 0.29 s (17x). Cell ids identical; all 41,162 polygons have identical vertices (rounded to 1e-6). A regional bbox fill (test extent, level 3) gave identical ids (678).

**Rejected:** `h3o::h3_to_vertexes()` for polygon construction — it returns only the 6 (5) corner vertices and drops the extra boundary vertices of cells crossing icosahedron edges (2,082 of 41,162 cells differed; areas off by up to ~11%). `sf::st_as_sfc()` on an H3 vector reproduces h3jsr's boundaries exactly.

**Notes:** res-0 ids are generated with `sprintf("8%03xfffffffffff", 2*(0:121)+1)`, checked equal to `h3jsr::get_res0()` in the benchmark script. Timed only the H3 steps; the s2 conservation-area intersection (same in both) dominates total build time at finer levels. Level 4 (current config) and the full `build_hex_grid()` not re-run with h3o yet.

## 2026-11-18 PST — Crop range polygons to <= 88 N and strip-cut features wider than 180 deg before h3o

**Decision:** `.polygon_cells()` (R/03_process/04a_range_from_polygons_fn.R) now (1) crops polygons to latitude <= `h3_max_lat` (88, set at the top of the file; north only), then (2) sends any feature whose longitude span exceeds 180 deg straight to the 60-degree strips; narrower features still try the direct h3o conversion first, with strips as the error fallback.

**Reason:** `run_build_ranges()` hung on mammals_marine *Ursus maritimus* (one valid 160-vertex ring, lon -179.999..179.999, top edge along 89.9 N spanning 360 deg). H3 treats edges as great circles, so that edge collapses to a sliver at the pole. A hang is not an error, so the old try-direct-then-strips fallback never ran. The polar cap is not needed for the analysis.

**Verified:** crop + strip geometry on the polar bear (no h3o call): bbox ymax 88, six valid POLYGON strips. **Not yet verified:** h3o conversion of the strips, and the resulting cell count; full `run_build_ranges()` not rerun.

**Notes:** The earlier session stalls (even `1 + 1`) were a wedged Positron kernel session, not the gpkg. The layer in spatial_data.gpkg is `exported`, not `spatial_data`. Cells above 88 N are now dropped for any species; the south is unchanged.

## 2026-11-18 PST — Follow-up: keep polygon parts only after st_intersection (WKB round trip) in .polygon_cells()

**Decision:** After each `st_intersection` (crop and strips), `polys_only()` extracts POLYGON parts, casts to MULTIPOLYGON and rebuilds the sfc via a WKB round trip before calling h3o.

**Reason:** The first crop version broke *Trichechus manatus* (59 features): `st_make_valid`/`st_intersection` return mixed `sfc_GEOMETRY`, which h3o rejects; `st_collection_extract` fixes the type but leaves nested ALTREP `wrap_list` objects (R 4.4) that h3o cannot take a writable pointer to. Manatus gives 761 cells with and without the crop (matches the raw direct conversion).

**Result:** `run_build_ranges(taxa = "mammals_marine", force = TRUE, progress = "verbose")` completed: 126 species, 1,540 s total (1,172 s in `.polygon_cells`), no warnings or errors in the log. *Ursus maritimus*: 11,750 cells in 0.1 s. Slowest: *Delphinus delphis* 174 s, *Phocoena phocoena* 144 s, *Orcinus orca* 69 s.

**Open:** ~370 s of the total is outside `.polygon_cells` (likely the verbose-mode `.area_km2`, which returns NA for several multi-feature species); not investigated.

## 2026-11-18 PST — Speed up .polygon_cells() (no up-front make_valid, per-polygon strips) and .area_km2() (equal-area planar)

**Investigation:** *Delphinus delphis* (3 features, 5,393 polygons, 3.4M vertices) and *Phocoena phocoena* (2,204 polygons, 2.1M vertices) took 174 s / 144 s. `st_make_valid` with s2 off (GEOS) was ~148 s / 124 s of that, although GEOS `st_is_valid` says all features are valid; the crop intersection took ~10-14 s and h3o on the exploded polygons 2-6 s. The `NA` area came from s2: `st_make_valid` (s2 on) output is rejected by `st_area` ("Loop ... crosses loop ..."), swallowed by `tryCatch`. Geodesic area with s2 off took ~150 s.

**Decision:** (1) `.polygon_cells()` crops the raw geometry and calls `st_make_valid` only if the crop errors; then explodes to single polygons and sends only polygons wider than 180 deg of longitude (circumpolar) to the 60-degree strips; narrower polygons go straight to h3o. (2) `.area_km2()` = planar area after `st_transform` to `config$crs$equal_area` (EPSG:6933): ~1 s, within 0.02% of geodesic on Delphinus/Phocoena.

**Verified:** cell counts identical to the previous run for Ursus maritimus (11,750), Trichechus manatus (761), Delphinus (46,396), Phocoena (6,691); `spp_cells_in_region.csv` and `spp_cells_in_cons_area.csv` byte-identical to the previous full run. Full mammals_marine build: 1,540 s -> ~305 s; slowest species now ~22 s (Delphinus).

**Notes:** Verbose areas from the old code were wrong for several whales (s2 area of make_valid output: e.g. Megaptera 163.8M km2 vs 346M km2 now; ~350M km2 expected from 197,574 cells x ~1,770 km2). Multi-feature areas are summed over features, so overlapping features can exceed the ocean area (B. musculus 858M km2); the figure is display only. A vctrs "structure(NULL, ...) is deprecated" warning (11 per run, from h3o with empty cell vectors) is harmless and not suppressed.

## 2026-11-19 19:19 PST — GBIF snapshot path now matches on alphanumeric specieskey, not species name
Snapshot (2026-11-01) has VARCHAR `specieskey`/`taxonkey` (e.g. 59PMB). Crosswalk rebuilt with matching alphanumeric `gbif_usage_key`. `run_gbif_snapshot()` filters on `specieskey IN (...)`; no name fallback. Crosswalk rows with no key (e.g. Lutra lutra, HIGHERRANK) are logged and not extracted. `load_crosswalk()` reads `gbif_usage_key` as character. Request hash now keyed on keys, so old snapshot checkpoints/caches are invalidated. Rejected: name matching (fragile to synonymy/spelling).

## 2026-11-19 19:50 PST — Trimmed gbif_columns_keep to the 11 columns actually used
Dropped occurrenceID, family, genus, scientificName, occurrenceStatus, publishingOrgKey, month, day (not referenced by parse, clean or any later stage). Kept datasetKey for the derived-dataset citation. Applies to config/taxa_config.R only; config/sample/taxa_config.R unchanged. Changes the snapshot request hash and parse fingerprint, so the snapshot rescans and parse reparses. Note: occurrenceStatus is not filtered anywhere, so ABSENT records are kept (open question). The snapshot code still adds `license` to the extract columns.

## 2026-11-20 11:01 PST — Manifest cache for build_species_crosswalk()
Stage 0 now skips when a settings hash matches `species_crosswalk_manifest.json` (next to the CSV) and the CSV exists. Hash covers: raw species list CONTENT (hashed after load, so Google Sheet edits are caught), each taxon's taxa_config entry, md5 of polygon files, rgbif/worrms versions, and the deparsed code of the crosswalk functions. Does NOT cover upstream backbone changes (GBIF/WoRMS/eBird); `force = TRUE` re-queries them. Reuses `.stage_cached()` / `.write_stage_manifest()`; manifest written after the CSV. Cached return reads the CSV, so all-NA columns come back as logical, not integer/character. Open: `run_pipeline()` calls it with `force = TRUE`, so the cache never applies there.

## 2026-11-20 11:01 PST — One combined OBIS S3 query across marine taxa
`run_obis_downloads()` now plans each taxon (cache check first), runs ONE S3 query for the union of ids/names of all pending taxa, splits results back by taxon, and falls back to the API per taxon if S3 fails, DuckDB is missing, or a taxon gets 0 rows. Species with no crosswalk `aphia_id` are excluded from OBIS (GBIF-only); logged per taxon (`n_no_aphia` in the result). Test (abalones + mammals_marine, force = TRUE, written to scratch): whole run 6.5 min -> 3.6 min (combined S3 query ~3.3 min vs ~3.2 + 3.0 min separate). `download_one_taxon_obis()` kept as a wrapper for `benchmark_downloads.R`.

## 2026-11-20 11:01 PST — OBIS S3 match on AphiaID OR species name
Tried AphiaID-only first (`interpreted.aphiaid` = crosswalk `aphia_id`): lost 54,338 mammals_marine and 281 abalone records vs the name-based download (OBIS `aphiaID` is often a subspecies/synonym id while `species` is rolled up; e.g. Phoca vitulina, Trichechus manatus, Enhydra lutris), but gained 7,661 under other names (Lobodon carcinophagus, Leucopleurus acutus). Rejected AphiaID-only. Adopted: keep a record if aphiaid IN ids OR species IN names (IN-subqueries, not an OR join). Result vs name-based: abalones 1,597,865 = 1,597,865; mammals_marine 2,003,077 -> 2,010,738 (0 lost, 7,661 gained). Query time ~3.3 min, no slowdown from the OR. Request hash now includes AphiaIDs, so existing OBIS caches are stale once. Real `data/raw/obis/` files not yet regenerated; needs `run_obis_downloads(force = TRUE)` then reparse/re-clean.
