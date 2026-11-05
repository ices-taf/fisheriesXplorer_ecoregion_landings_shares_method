# Extract production result, write FisheriesXplorer-ready output

library(icesTAF)
library(dplyr)

mkdir("output")

model_file <- "./model/stock_ecoregion_landings_shares.csv"
output_file <- "./output/stock_ecoregion_landings_shares.csv"

if (!file.exists(model_file)) {
  stop("Model output not found: ", model_file)
}

shares <- read.csv(
  model_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required <- c(
  "CLyear",
  "StockKeyLabel",
  "Ecoregion",
  "ShareLandingsInEcoregion",
  "Assigned_kg",
  "Inside_kg",
  "Missing_rectangle_kg",
  "share_method",
  "share_source",
  "special_rule",
  "ReferenceYear"
)
stopifnot(all(required %in% names(shares)))
stopifnot(!anyDuplicated(shares[c("CLyear", "StockKeyLabel", "Ecoregion")]))

shares <- shares |>
  arrange(Ecoregion, StockKeyLabel, CLyear)

write.csv(shares, output_file, row.names = FALSE, na = "")
message("Saved FisheriesXplorer-ready stock/ecoregion shares: ", output_file)
