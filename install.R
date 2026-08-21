# =============================================================================
# install.R  —  install the packages this pipeline needs
# =============================================================================
#
# Preferred: use renv to restore the exact pinned versions:
#   renv::restore()
#
# This script is the fallback / bootstrap when you are not using renv, or to
# seed the library before snapshotting. Run once:  Rscript install.R
#
# Note: h3jsr needs the V8 package (JavaScript engine). arrow, sf and terra
# pull in system libraries (GDAL/GEOS/PROJ, libarrow) — install those via your
# OS package manager if the R packages fail to build from source.

cran_pkgs <- c(
  # core / data
  "here", "dplyr", "tidyr", "readr", "purrr", "stringr", "glue", "tibble",
  "janitor", "data.table", "arrow",
  # spatial
  "sf", "terra", "h3jsr", "V8", "rnaturalearth",
  # biodiversity data APIs
  "rgbif", "rebird", "ebirdst", "rredlist", "taxize",
  # occurrence cleaning
  "CoordinateCleaner",
  # optional input source
  "googlesheets4"
)

installed <- rownames(installed.packages())
missing <- setdiff(cran_pkgs, installed)

if (length(missing) == 0) {
  message("All required packages already installed.")
} else {
  message("Installing: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}

# ebirdst is on CRAN but its data access also needs an access key in .Renviron
# (EBIRDST_KEY). See README.md -> Prerequisites.
message("\nDone. Next: copy .Renviron.example to .Renviron and add your keys,\n",
        "then see README.md for how to run the pipeline.")
