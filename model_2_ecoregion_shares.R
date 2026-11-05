# Production stock-level landings shares for all FisheriesXplorer ecoregions
#
# Input:
#   ./model/rdbes_stock_matches.RData created by model_1_global_stock_matching.R
#
# Output:
#   ./model/stock_ecoregion_landings_shares.csv
#
# Method:
#   - Fixed stock matching from the global checkpoint.
#   - For each ecoregion, a landing is inside when its ICES statistical rectangle
#     is returned by icesVMS::get_csquare(ecoregion = ..., convert2sf = FALSE).
#   - Share = 100 * Inside_kg / Assigned_kg.
#   - Missing rectangles count outside, matching the accepted rectangle convention.
#   - Sardine 25E4/25E5 is already applied upstream as a production stock rule.
#   - The pooled anchovy Division 9.a proxy is applied only for Biscay/Iberian waters.

library(dplyr)
library(tidyr)

# 1. Settings ---------------------------------------------------------------
if (!exists("supported_ecoregions", inherits = FALSE)) {
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
}
if (!exists("analysis_years", inherits = FALSE)) analysis_years <- NULL

checkpoint_file <- "./model/rdbes_stock_matches.RData"
cache_dir <- "./data/ices_cache_simple"
model_output_dir <- "./model/ecoregion_shares"
master_file <- "./model/stock_ecoregion_landings_shares.csv"

dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_output_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Load completed stock matching ------------------------------------------
message("[1/5] Reading global stock-matching checkpoint")
if (!file.exists(checkpoint_file)) stop("Checkpoint not found: ", checkpoint_file)

e <- new.env(parent = baseenv())
load(checkpoint_file, envir = e)
stopifnot(all(c(
  "dat", "candidates", "stock_list", "reference_year",
  "matching_scope", "supported_ecoregions"
) %in% ls(e)))

checkpoint_regions <- e$supported_ecoregions
if (!setequal(checkpoint_regions, supported_ecoregions)) {
  stop(
    "Checkpoint ecoregion scope differs from the requested production scope. ",
    "Rerun model_1_global_stock_matching.R."
  )
}

required_dat <- c(
  "row_id", "CLyear", "CLarea", "CLstatisticalRectangle",
  "CLspeciesFaoCode", "CLscientificWeight", "fishstock",
  "match_status", "special_rule"
)
stopifnot(all(required_dat %in% names(e$dat)))

dat <- select(e$dat, all_of(required_dat))
candidates <- select(e$candidates, row_id, fishstock) |> distinct()
catalogue <- e$stock_list
matching_reference_year <- as.integer(e$reference_year)
matching_scope <- e$matching_scope
rm(e)
invisible(gc())

if (!is.null(analysis_years)) {
  if (!all(analysis_years %in% dat$CLyear)) {
    stop("Requested analysis_years are absent from the checkpoint.")
  }
  dat <- filter(dat, CLyear %in% analysis_years)
}
years <- sort(unique(dat$CLyear))
stopifnot(
  nrow(dat) > 0L,
  !anyDuplicated(dat$row_id),
  !anyDuplicated(catalogue$fishstock)
)

# 3. Fixed 2026 stock-to-ecoregion catalogue --------------------------------
message("[2/5] Preparing fixed stock-to-ecoregion catalogue")

make_stock_ecoregion_catalogue <- function(stock_table, sid_year) {
  base <- stock_table |>
    transmute(
      fishstock = as.character(fishstock),
      description = as.character(description),
      SID_ecoregions = as.character(SID_ecoregions)
    ) |>
    filter(!is.na(SID_ecoregions), nzchar(trimws(SID_ecoregions))) |>
    separate_rows(SID_ecoregions, sep = ",") |>
    transmute(
      fishstock,
      description,
      Ecoregion = trimws(SID_ecoregions),
      SID_selection = "SID ecoregion",
      SIDYear = as.integer(sid_year)
    ) |>
    filter(Ecoregion %in% supported_ecoregions)

  # Accepted documented additions where the checked SID tag omitted Biscay.
  additions <- tibble(
    fishstock = c("mon.27.78abd", "rjb.27.89a"),
    Ecoregion = "Bay of Biscay and the Iberian Coast"
  ) |>
    semi_join(stock_table, by = "fishstock") |>
    left_join(
      select(stock_table, fishstock, description),
      by = "fishstock",
      relationship = "many-to-one"
    ) |>
    mutate(
      SID_selection = "documented addition",
      SIDYear = as.integer(sid_year)
    ) |>
    select(fishstock, description, Ecoregion, SID_selection, SIDYear)

  bind_rows(base, additions) |>
    arrange(fishstock, Ecoregion, SID_selection) |>
    distinct(fishstock, Ecoregion, .keep_all = TRUE)
}

fixed_catalogue_long <- make_stock_ecoregion_catalogue(
  catalogue,
  matching_reference_year
)

write.csv(
  fixed_catalogue_long,
  file.path(model_output_dir, paste0(
    "stock_ecoregion_catalogue_", matching_reference_year, ".csv"
  )),
  row.names = FALSE
)

missing_sid_regions <- catalogue |>
  filter(is.na(SID_ecoregions) | !nzchar(trimws(SID_ecoregions))) |>
  select(fishstock, description)
write.csv(
  missing_sid_regions,
  file.path(model_output_dir, "stocks_missing_SID_ecoregion.csv"),
  row.names = FALSE
)

# 4. Prepare global stock totals once ---------------------------------------
assigned <- dat |>
  filter(match_status == "unique", !is.na(fishstock)) |>
  mutate(
    missing_rectangle = is.na(CLstatisticalRectangle) |
      !nzchar(trimws(CLstatisticalRectangle))
  )

# If the sardine production rule was used in a year, it affects both the stock
# receiving the records and the stock from which those records were removed.
sardine_years <- sort(unique(
  assigned$CLyear[assigned$special_rule == "sardine_25E4_25E5"]
))
sardine_impacts <- if (length(sardine_years)) {
  tidyr::crossing(
    CLyear = sardine_years,
    fishstock = c("pil.27.7", "pil.27.8abd")
  ) |>
    mutate(rule_impact = "sardine_25E4_25E5")
} else {
  tibble(CLyear = integer(), fishstock = character(), rule_impact = character())
}

stock_totals <- assigned |>
  summarise(
    Assigned_kg = sum(CLscientificWeight),
    Missing_rectangle_kg = sum(CLscientificWeight[missing_rectangle]),
    .by = c(CLyear, fishstock)
  ) |>
  left_join(
    sardine_impacts,
    by = c("CLyear", "fishstock"),
    relationship = "one-to-one"
  ) |>
  mutate(special_rule = coalesce(rule_impact, "none")) |>
  select(-rule_impact)

# Retain unresolved candidate exposure as a separate audit. It is overlapping
# candidate exposure and must never be added across stocks.
unresolved_candidate_exposure <- candidates |>
  inner_join(
    dat |>
      filter(match_status != "unique") |>
      select(row_id, CLyear, CLscientificWeight),
    by = "row_id",
    relationship = "many-to-many"
  ) |>
  summarise(
    Unresolved_candidate_kg = sum(CLscientificWeight),
    .by = c(CLyear, fishstock)
  )
write.csv(
  unresolved_candidate_exposure,
  file.path(model_output_dir, "unresolved_candidate_exposure.csv"),
  row.names = FALSE
)

# Fetch/cache one rectangle vector per ecoregion.
get_ecoregion_rectangles <- function(ecoregion) {
  path <- file.path(
    cache_dir,
    paste0("vms_rectangle_ids_", gsub("[^A-Za-z0-9]", "_", ecoregion), ".rds")
  )
  if (file.exists(path)) {
    area <- readRDS(path)
  } else {
    cells <- icesVMS::get_csquare(ecoregion = ecoregion, convert2sf = FALSE)
    stopifnot(is.data.frame(cells), "stat_rec" %in% names(cells))
    area <- unique(as.character(cells$stat_rec))
    area <- area[!is.na(area) & nzchar(trimws(area))]
    if (!length(area)) stop("No statistical rectangles returned for: ", ecoregion)
    saveRDS(area, path)
  }
  stopifnot(is.character(area), length(area) > 0L, !anyNA(area))
  unique(area)
}

calculate_region <- function(ecoregion) {
  message("  Ecoregion: ", ecoregion)
  area <- get_ecoregion_rectangles(ecoregion)
  regional_catalogue <- fixed_catalogue_long |>
    filter(Ecoregion == ecoregion)

  if (!nrow(regional_catalogue)) {
    warning("No stocks in fixed catalogue for ecoregion: ", ecoregion)
    return(tibble())
  }

  inside_totals <- assigned |>
    filter(CLstatisticalRectangle %in% area) |>
    summarise(
      Inside_kg = sum(CLscientificWeight),
      .by = c(CLyear, fishstock)
    )

  out <- crossing(
    CLyear = years,
    fishstock = regional_catalogue$fishstock
  ) |>
    left_join(stock_totals, by = c("CLyear", "fishstock")) |>
    left_join(inside_totals, by = c("CLyear", "fishstock")) |>
    mutate(
      Assigned_kg = coalesce(Assigned_kg, 0),
      Inside_kg = coalesce(Inside_kg, 0),
      Missing_rectangle_kg = coalesce(Missing_rectangle_kg, 0),
      special_rule = coalesce(special_rule, "none"),
      ShareLandingsInEcoregion = if_else(
        Assigned_kg > 0,
        100 * Inside_kg / Assigned_kg,
        NA_real_
      ),
      share_method = if_else(
        is.na(ShareLandingsInEcoregion),
        "unavailable",
        "uniquely_assigned_stock"
      ),
      share_source = fishstock,
      Ecoregion = ecoregion,
      ReferenceYear = matching_reference_year
    ) |>
    left_join(
      select(regional_catalogue, fishstock, description, SID_selection),
      by = "fishstock",
      relationship = "many-to-one"
    )

  # Accepted production proxy: the two split Division 9.a anchovy stocks cannot
  # be distinguished by the available broad area relationship. For Biscay only,
  # fill a missing stock-specific share with the annual pooled Division 9.a share.
  if (ecoregion == "Bay of Biscay and the Iberian Coast") {
    targets <- c("ane.27.9aS", "ane.27.9aW")

    anchovy_records <- dat |>
      filter(
        toupper(trimws(CLspeciesFaoCode)) == "ANE",
        !is.na(CLarea),
        tolower(trimws(CLarea)) == "27.9.a" |
          startsWith(tolower(trimws(CLarea)), "27.9.a.")
      ) |>
      mutate(
        pool_eligible = !is.na(match_status) &
          match_status != "unresolved_taxon",
        inside_region = CLstatisticalRectangle %in% area,
        missing_rectangle = is.na(CLstatisticalRectangle) |
          !nzchar(trimws(CLstatisticalRectangle))
      )

    anchovy_pool <- anchovy_records |>
      summarise(
        Pool_records = sum(pool_eligible),
        Pool_kg = sum(CLscientificWeight[pool_eligible]),
        Pool_inside_kg = sum(CLscientificWeight[pool_eligible & inside_region]),
        Pool_missing_rectangle_kg = sum(
          CLscientificWeight[pool_eligible & missing_rectangle]
        ),
        Excluded_taxon_or_status_kg = sum(
          CLscientificWeight[!pool_eligible]
        ),
        .by = CLyear
      ) |>
      mutate(
        PooledShare_pct = if_else(
          Pool_kg > 0,
          100 * Pool_inside_kg / Pool_kg,
          NA_real_
        ),
        Ecoregion = ecoregion,
        share_source = "ANE_Division_9a"
      )

    write.csv(
      anchovy_pool,
      file.path(model_output_dir, "anchovy_9a_proxy_audit.csv"),
      row.names = FALSE
    )

    out <- out |>
      left_join(
        select(anchovy_pool, CLyear, PooledShare_pct),
        by = "CLyear",
        relationship = "many-to-one"
      ) |>
      mutate(
        use_proxy = fishstock %in% targets &
          is.na(ShareLandingsInEcoregion) &
          !is.na(PooledShare_pct),
        ShareLandingsInEcoregion = if_else(
          use_proxy,
          PooledShare_pct,
          ShareLandingsInEcoregion
        ),
        share_method = if_else(
          use_proxy,
          "pooled_9a_proxy",
          share_method
        ),
        share_source = if_else(
          use_proxy,
          "ANE_Division_9a",
          share_source
        ),
        special_rule = if_else(
          use_proxy,
          "anchovy_9a_proxy",
          special_rule
        ),
        # These pooled kilograms are not stock-specific and therefore are not
        # duplicated into the two split-stock rows. See anchovy_9a_proxy_audit.csv.
        Assigned_kg = if_else(use_proxy, NA_real_, Assigned_kg),
        Inside_kg = if_else(use_proxy, NA_real_, Inside_kg),
        Missing_rectangle_kg = if_else(
          use_proxy,
          NA_real_,
          Missing_rectangle_kg
        )
      ) |>
      select(-PooledShare_pct, -use_proxy)
  }

  out
}

# 5. Calculate all regions and validate -------------------------------------
message("[3/5] Calculating production shares across all ecoregions")
master <- bind_rows(lapply(supported_ecoregions, calculate_region))
if (!nrow(master)) stop("No regional share rows were produced.")

# Validate final percentage values with a small tolerance for machine rounding.
p <- master$ShareLandingsInEcoregion
tolerance_pct <- 1e-8
invalid <- is.nan(p) | (!is.na(p) & (
  !is.finite(p) | p < -tolerance_pct | p > 100 + tolerance_pct
))
if (any(invalid)) {
  invalid_rows <- master[invalid, , drop = FALSE]
  invalid_rows$Share_full_precision <- sprintf("%.17g", p[invalid])
  write.csv(
    invalid_rows,
    file.path(model_output_dir, "invalid_share_values.csv"),
    row.names = FALSE
  )
  stop("Invalid regional share values remain; inspect invalid_share_values.csv")
}

roundoff <- !is.na(p) & (p < 0 | p > 100)
if (any(roundoff)) {
  adjustments <- master[roundoff, , drop = FALSE]
  adjustments$Share_before_full_precision <- sprintf("%.17g", p[roundoff])
  adjustments$Share_after_pct <- pmin(100, pmax(0, p[roundoff]))
  write.csv(
    adjustments,
    file.path(model_output_dir, "share_rounding_adjustments.csv"),
    row.names = FALSE
  )
  p[roundoff] <- adjustments$Share_after_pct
}
master$ShareLandingsInEcoregion <- p

stopifnot(!anyDuplicated(master[c("CLyear", "fishstock", "Ecoregion")]))

# Production master table. Scenario is intentionally omitted: special-rule
# information is row-specific through share_method, share_source and special_rule.
master_output <- master |>
  transmute(
    CLyear,
    StockKeyLabel = fishstock,
    Ecoregion,
    ShareLandingsInEcoregion,
    Assigned_kg,
    Inside_kg,
    Missing_rectangle_kg,
    share_method,
    share_source,
    special_rule,
    ReferenceYear
  ) |>
  arrange(Ecoregion, StockKeyLabel, CLyear)

write.csv(master_output, master_file, row.names = FALSE, na = "")

# 6. Annual SID audit only ---------------------------------------------------
message("[4/5] Auditing annual SID stock codes and ecoregion membership")
read_sid <- function(y) {
  path <- file.path(cache_dir, paste0("stocks_", y, ".rds"))
  if (file.exists(path)) {
    x <- readRDS(path)
  } else {
    x <- icesSD::getSD(year = y)
    stopifnot(is.data.frame(x), nrow(x) > 0L)
    saveRDS(x, path)
  }
  stopifnot(
    nrow(x) > 0L,
    all(c("StockKeyLabel", "StockKeyDescription", "EcoRegion") %in% names(x))
  )
  x
}

annual_sid_list <- lapply(years, function(y) {
  tryCatch({
    x <- read_sid(y) |>
      transmute(
        fishstock = as.character(StockKeyLabel),
        description = as.character(StockKeyDescription),
        SID_ecoregions = na_if(trimws(as.character(EcoRegion)), "")
      ) |>
      distinct(fishstock, .keep_all = TRUE)
    make_stock_ecoregion_catalogue(x, y)
  }, error = function(err) {
    warning("SID ", y, " audit unavailable: ", conditionMessage(err))
    NULL
  })
})

annual_sid <- bind_rows(annual_sid_list)
if (nrow(annual_sid)) {
  write.csv(
    annual_sid,
    file.path(model_output_dir, "SID_catalogue_by_year_long.csv"),
    row.names = FALSE
  )

  sid_counts <- annual_sid |>
    summarise(
      Regional_codes = n_distinct(fishstock),
      .by = c(SIDYear, Ecoregion)
    )
  write.csv(
    sid_counts,
    file.path(model_output_dir, "SID_regional_counts_by_year.csv"),
    row.names = FALSE
  )

  # Exact stock-code audit: fixed-reference matching is not changed by this.
  sid_codes_by_year <- annual_sid |>
    distinct(CLyear = SIDYear, fishstock)
  absent_codes <- assigned |>
    distinct(CLyear, fishstock) |>
    anti_join(sid_codes_by_year, by = c("CLyear", "fishstock"))
  write.csv(
    absent_codes,
    file.path(model_output_dir, "assigned_codes_absent_from_yearly_SID.csv"),
    row.names = FALSE
  )

  # Fixed 2026 stock-region membership versus each historical year's SID.
  membership_audit <- master_output |>
    distinct(CLyear, StockKeyLabel, Ecoregion) |>
    anti_join(
      annual_sid |>
        transmute(
          CLyear = SIDYear,
          StockKeyLabel = fishstock,
          Ecoregion
        ),
      by = c("CLyear", "StockKeyLabel", "Ecoregion")
    )
  write.csv(
    membership_audit,
    file.path(model_output_dir, "fixed_membership_absent_from_yearly_SID.csv"),
    row.names = FALSE
  )
}

writeLines(
  c(
    paste("Checkpoint:", normalizePath(checkpoint_file)),
    paste("Fixed matching SID year:", matching_reference_year),
    paste("RDBES years:", paste(years, collapse = ", ")),
    paste("Ecoregions:", paste(supported_ecoregions, collapse = "; ")),
    "Production share basis: ICES statistical-rectangle membership from icesVMS.",
    "Denominator: all uniquely assigned RDBES landings for the stock/year across all areas.",
    "Missing statistical rectangles count outside the ecoregion and are reported separately.",
    "Sardine 25E4/25E5 reassignment is an accepted upstream production stock-assignment rule.",
    "Anchovy ane.27.9aS / ane.27.9aW use a pooled Division 9.a proxy only for Bay of Biscay and the Iberian Coast when stock-specific shares are unavailable.",
    "Pooled anchovy kilograms are not duplicated into split-stock Assigned_kg/Inside_kg fields; those fields are NA on proxy rows.",
    "Annual SID is an audit only; stock assignments and final stock-region membership use the fixed reference catalogue.",
    "Ecoregion rectangle lists can overlap, so shares for a stock are not required to sum to 100 across ecoregions.",
    paste("Master output:", normalizePath(master_file))
  ),
  file.path(model_output_dir, "RUN_NOTES.txt")
)

message("[5/5] Saved production master table: ", master_file)
message("  Rows: ", nrow(master_output))
message("  Unique stocks: ", n_distinct(master_output$StockKeyLabel))
message("  Ecoregions: ", n_distinct(master_output$Ecoregion))
