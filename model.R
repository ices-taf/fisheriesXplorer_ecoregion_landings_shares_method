# Run analysis, write model results

library(icesTAF)

mkdir("model")

# Fixed production settings shared by both modelling stages.
supported_ecoregions <- c(
  "Baltic Sea",
  "Bay of Biscay and the Iberian Coast",
  "Celtic Seas",
  "Greater North Sea",
  "Norwegian Sea",
  "Icelandic Waters",
  "Barents Sea",
  "Greenland Sea",
  "Faroes",
  "Oceanic Northeast Atlantic",
  "Azores"
)

reference_year <- 2026L
analysis_years <- NULL
rebuild_landings <- FALSE
refresh_reference <- FALSE

source("model_1_global_stock_matching_conflict_refined.R", local = FALSE)
source("model_2_ecoregion_shares.R", local = FALSE)
