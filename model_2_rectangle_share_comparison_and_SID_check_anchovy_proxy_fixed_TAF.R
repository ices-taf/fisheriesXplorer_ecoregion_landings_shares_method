# Re-export completed RDBES stock matches using the ORIGINAL rectangle rule.
# Also inspect yearly SID lists. This does NOT perform annual stock rematching.
# It does not run GIS overlays or read extendedRDBES.RData.
# Source this file on its own; do not source v5 or its old companion first.
# Additional export: a share-only lookup with a labelled annual pooled 9.a
# anchovy proxy. It does NOT allocate pooled RDBES kilograms to the split stocks.
# Existing stock summaries, coverage comparisons and checkpoint are unchanged.

library(dplyr)
library(tidyr)

# ---- Settings: edit paths, not the saved matching reference year ---------
checkpoint_file <- "./data/rdbes_stock_matches.RData"
original_file <- "./output/data_03_shareLandingsInEcoregion.csv"  # optional legacy comparison; skipped if absent
Ecoregion <- "Bay of Biscay and the Iberian Coast"
analysis_years <- NULL             # NULL = every year in the saved checkpoint
cache_dir <- "./data/ices_cache_simple"
original_rectangle_file <- NULL    # Optional RDS containing the ORIGINAL 'area' vector
# Otherwise fetch/cache the current VMS list using the original query convention.
# A newly fetched list is not proof of an identical historical service response.
output_root <- "./output/stock_matching_rectangle_check"
use_anchovy_9a_proxy <- TRUE  # Configured for Bay of Biscay and the Iberian Coast only.

run_rectangle_check <- function() {
  if (use_anchovy_9a_proxy && Ecoregion != "Bay of Biscay and the Iberian Coast") {
    stop("The pooled 9.a proxy is configured for Biscay/Iberian waters only. ",
         "Set use_anchovy_9a_proxy <- FALSE for a different ecoregion.")
  }
  if (!file.exists(checkpoint_file)) stop("Checkpoint not found: ", checkpoint_file)
  comparison_available <- !is.null(original_file) && length(original_file) == 1L &&
    !is.na(original_file) && nzchar(original_file) && file.exists(original_file)
  if (!comparison_available)
    message("  Legacy IFREMER comparison file not found; continuing without comparison audit.")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- file.path(output_root, gsub("[^A-Za-z0-9]", "_", Ecoregion),
                          format(Sys.time(), "%Y%m%d_%H%M%S"))
  if (dir.exists(output_dir)) stop("Output directory already exists; rerun with a new timestamp.")
  dir.create(output_dir, recursive = TRUE)
  write_table <- function(x, name) write.csv(x, file.path(output_dir, name), row.names = FALSE, na = "")

  message("[1/5] Reading completed stock matching; no rematching")
  e <- new.env(parent = baseenv())
  load(checkpoint_file, envir = e)
  stopifnot(all(c("dat", "candidates", "stock_list", "reference_year") %in% ls(e)))
  required <- c("row_id", "CLyear", "CLarea", "CLstatisticalRectangle", "CLspeciesFaoCode",
                "CLscientificWeight", "fishstock", "match_status")
  stopifnot(all(required %in% names(e$dat)))
  dat <- select(e$dat, all_of(required))
  candidates <- select(e$candidates, row_id, fishstock) |> distinct()
  catalogue <- e$stock_list
  matching_reference_year <- as.integer(e$reference_year)
  scope <- if (exists("matching_scope", envir = e, inherits = FALSE)) e$matching_scope else NULL
  rm(e)
  invisible(gc())
  stopifnot(length(matching_reference_year) == 1L, !is.na(matching_reference_year))
  if (!is.null(scope$ecoregions) && !Ecoregion %in% scope$ecoregions)
    stop("Checkpoint does not cover this ecoregion. Use the completed global checkpoint.")
  if (!is.null(analysis_years)) {
    if (!all(analysis_years %in% dat$CLyear)) stop("Requested years are absent from this checkpoint.")
    dat <- filter(dat, CLyear %in% analysis_years)
  }
  stopifnot(nrow(dat) > 0L, !anyNA(dat$CLyear), !anyNA(dat$CLscientificWeight),
            all(is.finite(dat$CLscientificWeight)), all(dat$CLscientificWeight >= 0),
            !anyDuplicated(catalogue$fishstock))
  years <- sort(unique(dat$CLyear))
  message("  Saved matching uses SID ", matching_reference_year,
          "; catch years: ", paste(years, collapse = ", "))

  read_sid <- function(y) {
    p <- file.path(cache_dir, paste0("stocks_", y, ".rds"))
    if (file.exists(p)) x <- readRDS(p) else {
      x <- icesSD::getSD(year = y)
      stopifnot(is.data.frame(x), nrow(x) > 0L)
      saveRDS(x, p)
    }
    stopifnot(nrow(x) > 0L, all(c("StockKeyLabel", "StockKeyDescription", "EcoRegion") %in% names(x)))
    if ("ActiveYear" %in% names(x)) {
      returned <- as.integer(as.character(x$ActiveYear))
      if (anyNA(returned) || any(returned != y)) stop("Wrong ActiveYear in SID cache: ", p)
    }
    x <- transmute(x, SIDYear = as.integer(y), fishstock = as.character(StockKeyLabel),
                   description = as.character(StockKeyDescription),
                   SID_ecoregions = na_if(trimws(as.character(EcoRegion)), "")) |> distinct()
    if (anyNA(x$fishstock) || anyDuplicated(x$fishstock))
      stop("Missing codes or conflicting duplicate stock metadata in SID ", y)
    x
  }
  # Backwards compatibility with v4 checkpoints; do not change their reference year.
  if (!"SID_ecoregions" %in% names(catalogue)) {
    catalogue <- left_join(catalogue,
      select(read_sid(matching_reference_year), fishstock, SID_ecoregions),
      by = "fishstock", relationship = "one-to-one")
  }
  select_region <- function(x) {
    tags <- strsplit(coalesce(x$SID_ecoregions, ""), ",", fixed = TRUE)
    # Retain the same explicit additions used in the supplied v5, not new overrides.
    additions <- if (Ecoregion == "Bay of Biscay and the Iberian Coast")
      c("mon.27.78abd", "rjb.27.89a") else character()
    x |> mutate(SID_selection = case_when(
      vapply(tags, function(z) Ecoregion %in% trimws(z), logical(1)) ~ "SID ecoregion",
      fishstock %in% additions ~ "v5 documented addition",
      is.na(SID_ecoregions) | SID_ecoregions == "" ~ "SID ecoregion missing",
      TRUE ~ NA_character_))
  }
  catalogue <- select_region(catalogue)
  selected_codes <- catalogue$fishstock[!is.na(catalogue$SID_selection)]

  message("[2/5] Reading the ecoregion rectangle list; no polygons")
  rectangle_source <- original_rectangle_file
  if (is.null(rectangle_source)) rectangle_source <- file.path(cache_dir,
    paste0("vms_rectangle_ids_", gsub("[^A-Za-z0-9]", "_", Ecoregion), ".rds"))
  if (file.exists(rectangle_source)) {
    area <- readRDS(rectangle_source)
  } else {
    if (!is.null(original_rectangle_file)) stop("Original rectangle vector not found: ", rectangle_source)
    cells <- icesVMS::get_csquare(ecoregion = Ecoregion, convert2sf = FALSE)
    stopifnot(is.data.frame(cells), "stat_rec" %in% names(cells))
    area <- unique(as.character(cells$stat_rec))
    area <- area[!is.na(area) & nzchar(trimws(area))]
    if (!length(area)) stop("VMS returned no usable rectangles.")
    saveRDS(area, rectangle_source)
    rm(cells)
  }
  stopifnot(is.character(area), length(area) > 0L, !anyNA(area), all(nzchar(trimws(area))))
  # Do not replace records' coordinates or stock identities to determine location.
  dat$inside_original_rule <- dat$CLstatisticalRectangle %in% area
  # Missing rectangles are outside under the original convention, but their weight is visible.
  dat$missing_rectangle <- is.na(dat$CLstatisticalRectangle) | !nzchar(trimws(dat$CLstatisticalRectangle))
  assigned <- filter(dat, match_status == "unique", !is.na(fishstock))
  stopifnot(!anyDuplicated(dat$row_id), all(assigned$fishstock %in% catalogue$fishstock))
  # Candidate exposure is not allocated landings and must NOT be summed across stocks.
  exposure <- candidates |>
    inner_join(dat |> filter(match_status != "unique") |>
      select(row_id, CLyear, CLscientificWeight), by = "row_id", relationship = "many-to-many") |>
    summarise(Unresolved_candidate_kg = sum(CLscientificWeight), .by = c(CLyear, fishstock))
  rm(candidates)
  make_summary <- function(x, scenario) {
    totals <- x |> summarise(Assigned_kg = sum(CLscientificWeight),
      Inside_kg = sum(CLscientificWeight[inside_original_rule]),
      Missing_rectangle_kg = sum(CLscientificWeight[missing_rectangle]), .by = c(CLyear, fishstock))
    out <- crossing(CLyear = years, fishstock = catalogue$fishstock) |>
      left_join(totals, by = c("CLyear", "fishstock")) |>
      left_join(exposure, by = c("CLyear", "fishstock")) |>
      mutate(across(ends_with("_kg"), ~ coalesce(.x, 0)), Outside_kg = Assigned_kg - Inside_kg,
        ShareLandingsInEcoregion = if_else(Assigned_kg > 0, 100 * Inside_kg / Assigned_kg, NA_real_),
        Scenario = scenario, ReferenceYear = matching_reference_year,
        StockReferenceMode = "fixed_completed_checkpoint", Ecoregion = Ecoregion,
        ShareBasis = "original_rectangle_rule_all_assigned_weight") |>
      left_join(catalogue, by = "fishstock", relationship = "many-to-one")
    stopifnot(isTRUE(all.equal(sum(out$Assigned_kg), sum(x$CLscientificWeight))))
    filter(out, !is.na(SID_selection))
  }

  message("[3/5] Recalculating shares and a SEPARATE sardine sensitivity")
  baseline <- make_summary(assigned, "unchanged_stock_assignments")
  write_table(baseline, "stock_summary_rectangles.csv")
  write_table(dat |> summarise(Landings_kg = sum(CLscientificWeight),
                              .by = c(CLyear, match_status)), "matching_summary.csv")
  summaries <- list(baseline)
  # Country was not retained in the supplied compact table. This is the IFREMER
  # rectangle convention, NOT a claim that the French-catch condition is verified.
  transfer <- assigned$CLspeciesFaoCode %in% "PIL" & assigned$fishstock == "pil.27.7" &
              assigned$CLstatisticalRectangle %in% c("25E4", "25E5")
  if (any(transfer) && "pil.27.8abd" %in% catalogue$fishstock) {
    transfer_log <- assigned[transfer, ] |>
      summarise(Landings_kg = sum(CLscientificWeight),
        .by = c(CLyear, CLarea, CLstatisticalRectangle, fishstock, inside_original_rule)) |>
      mutate(Proposed_stock = "pil.27.8abd", Country_condition_verified = FALSE)
    write_table(transfer_log, "sardine_IFREMER_sensitivity_transfers.csv")
    sensitivity <- assigned
    sensitivity$fishstock[transfer] <- "pil.27.8abd"
    stopifnot(isTRUE(all.equal(sum(sensitivity$CLscientificWeight), sum(assigned$CLscientificWeight))))
    summaries[[2L]] <- make_summary(sensitivity, "IFREMER_sardine_rectangle_sensitivity")
    write_table(summaries[[2L]], "stock_summary_rectangles_IFREMER_sensitivity.csv")
    rm(sensitivity)
  }

  # 3b. Share-only lookup for downstream regional estimates -----------------
  # Prefer the separate sardine sensitivity when it was produced above;
  # otherwise use the unchanged-assignment baseline. Scenario records this.
  lookup_base <- summaries[[length(summaries)]]
  share_lookup_export <- lookup_base |>
    select(CLyear, fishstock, ShareLandingsInEcoregion,
           Ecoregion, ReferenceYear, Scenario) |>
    mutate(
      share_method = if_else(is.na(ShareLandingsInEcoregion),
                             "unavailable", "uniquely_assigned_stock"),
      share_source = fishstock
    )
  anchovy_9a_pool <- tibble()

  if (use_anchovy_9a_proxy) {
    message("  Exporting annual pooled 9.a anchovy proxy; stock assignments stay unchanged")
    targets <- c("ane.27.9aS", "ane.27.9aW")
    stopifnot(all(targets %in% share_lookup_export$fishstock))

    # Use dat, NOT assigned or the candidate table. Ambiguous stock identity
    # is allowed here; unresolved/conflicting taxonomy is not. Each row_id
    # represents one aggregate and was checked for uniqueness above.
    anchovy_records <- dat |>
      filter(
        toupper(trimws(CLspeciesFaoCode)) == "ANE",
        !is.na(CLarea),
        tolower(trimws(CLarea)) == "27.9.a" |
          startsWith(tolower(trimws(CLarea)), "27.9.a.")
      ) |>
      mutate(pool_eligible = !is.na(match_status) &
               match_status != "unresolved_taxon")

    anchovy_9a_pool <- anchovy_records |>
      summarise(
        Pool_records = sum(pool_eligible),
        Pool_kg = sum(CLscientificWeight[pool_eligible]),
        Pool_inside_kg = sum(CLscientificWeight[pool_eligible & inside_original_rule]),
        Pool_missing_rectangle_kg = sum(CLscientificWeight[pool_eligible & missing_rectangle]),
        Excluded_taxon_or_status_kg = sum(CLscientificWeight[!pool_eligible]),
        .by = CLyear
      ) |>
      mutate(
        PooledShare_pct = if_else(Pool_kg > 0,
          100 * Pool_inside_kg / Pool_kg, NA_real_),
        Ecoregion = Ecoregion,
        share_source = "ANE_Division_9a"
      )
    # Save pooled weights ONCE per year, separately from the stock share lookup.
    write_table(anchovy_9a_pool, "anchovy_9a_proxy_audit.csv")

    share_lookup_export <- share_lookup_export |>
      left_join(select(anchovy_9a_pool, CLyear, PooledShare_pct),
                by = "CLyear", relationship = "many-to-one") |>
      mutate(
        # Fill missing shares only; preserve any existing stock-specific result.
        use_proxy = fishstock %in% targets &
                    is.na(ShareLandingsInEcoregion) & !is.na(PooledShare_pct),
        ShareLandingsInEcoregion = if_else(use_proxy,
          PooledShare_pct, ShareLandingsInEcoregion),
        share_method = if_else(use_proxy, "pooled_9a_proxy", share_method),
        share_source = if_else(use_proxy, "ANE_Division_9a", share_source)
      ) |>
      select(-PooledShare_pct, -use_proxy)
    if (!any(share_lookup_export$share_method == "pooled_9a_proxy"))
      warning("No missing split-stock shares were filled. Inspect anchovy_9a_proxy_audit.csv.")
    rm(anchovy_records)
  }

  # Validate percentages without rejecting machine-level boundary round-off.
  # Units remain PERCENT (0-100), not proportions (0-1).
  stopifnot(!anyDuplicated(share_lookup_export[c("CLyear", "fishstock")]))
  p <- share_lookup_export$ShareLandingsInEcoregion
  stopifnot(is.numeric(p))
  tolerance_pct <- 1e-8
  invalid <- is.nan(p) | (!is.na(p) &
    (!is.finite(p) | p < -tolerance_pct | p > 100 + tolerance_pct))

  # Fail with the actual stock-year values; never hide substantive errors.
  if (any(invalid)) {
    invalid_rows <- share_lookup_export[invalid, , drop = FALSE]
    invalid_rows$Share_full_precision <- sprintf("%.17g", p[invalid])
    invalid_file <- file.path(output_dir, "invalid_share_lookup_values.csv")
    write_table(invalid_rows, "invalid_share_lookup_values.csv")
    print(as.data.frame(invalid_rows), row.names = FALSE, digits = 17)
    stop("Invalid percentages remain; no share lookup was saved. Inspect: ",
         invalid_file, call. = FALSE)
  }

  # Only snap values JUST outside the endpoints, after the validation above.
  # Interior values and unavailable (NA) shares remain unchanged.
  roundoff <- !is.na(p) & (p < 0 | p > 100)
  if (any(roundoff)) {
    rounding_rows <- share_lookup_export[roundoff, , drop = FALSE]
    rounding_rows$Share_before_full_precision <- sprintf("%.17g", p[roundoff])
    rounding_rows$Share_after_pct <- pmin(100, pmax(0, p[roundoff]))
    write_table(rounding_rows, "share_rounding_adjustments.csv")
    p[roundoff] <- rounding_rows$Share_after_pct
    message("  Normalised ", sum(roundoff),
            " percentage(s) within ", tolerance_pct,
            " percentage points of a boundary; see share_rounding_adjustments.csv")
  }
  share_lookup_export$ShareLandingsInEcoregion <- p
  stopifnot(all(is.na(p) | (is.finite(p) & p >= 0 & p <= 100)))
  lookup_filename <- if (use_anchovy_9a_proxy)
    "stock_share_lookup_with_anchovy_proxy.csv" else "stock_share_lookup.csv"
  write_table(share_lookup_export, lookup_filename)
  message("  Share lookup base scenario: ", unique(share_lookup_export$Scenario))
  # Do NOT add proxy rows to summaries: the existing comparisons describe
  # stock-resolved results, not the separate proxy-assisted plotting lookup.

  comparison <- tibble()
  coverage <- tibble()

  if (comparison_available) {
    message("[4/5] Comparing exact stock codes with the optional legacy CSV")
    original <- read.csv(original_file, stringsAsFactors = FALSE, check.names = FALSE)
    stopifnot(all(c("CLyear", "fishstock", "CS", "Landings") %in% names(original)))
    original <- filter(original, CLyear %in% years, !is.na(fishstock), nzchar(fishstock))
    stopifnot(!anyNA(original$Landings), all(original$CS %in% c("inside", "outside")))
    original <- original |> summarise(OriginalAssigned_kg = sum(Landings),
      OriginalInside_kg = sum(Landings[CS == "inside"]), .by = c(CLyear, fishstock)) |>
      mutate(OriginalShare_pct = if_else(OriginalAssigned_kg > 0,
                                       100 * OriginalInside_kg / OriginalAssigned_kg, NA_real_))
    comparison <- bind_rows(lapply(summaries, function(s) {
      scenario <- unique(s$Scenario)
      full_join(original, select(s, CLyear, fishstock, NewAssigned_kg = Assigned_kg,
        NewInside_kg = Inside_kg, NewShare_pct = ShareLandingsInEcoregion),
        by = c("CLyear", "fishstock"), relationship = "one-to-one") |>
        mutate(Scenario = scenario, In_fixed_regional_catalogue = fishstock %in% selected_codes,
          Original_positive = coalesce(OriginalShare_pct > 0, FALSE),
          New_positive = coalesce(NewShare_pct > 0, FALSE),
          Difference_pp = NewShare_pct - OriginalShare_pct,
          Coverage = case_when(Original_positive & New_positive ~ "positive_in_both",
            New_positive ~ "positive_only_in_new", Original_positive ~ "positive_only_in_original",
            TRUE ~ "no_positive_share_in_either"))
    }))
    coverage <- comparison |> summarise(
      Original_positive_all_codes = sum(Original_positive), New_positive_selected_codes = sum(New_positive),
      Original_positive_same_catalogue = sum(Original_positive & In_fixed_regional_catalogue),
      New_numeric_shares = sum(!is.na(NewShare_pct)),
      Positive_both = sum(Original_positive & New_positive),
      Positive_only_original = sum(Original_positive & !New_positive),
      Positive_only_new = sum(!Original_positive & New_positive), .by = c(Scenario, CLyear))
    write_table(comparison, "stock_coverage_comparison.csv")
    write_table(coverage, "coverage_by_year.csv")
  } else {
    message("[4/5] Legacy comparison skipped")
  }

  message("[5/5] Auditing annual SID lists; NOT changing saved stock assignments")
  sid_results <- lapply(years, function(y) {
    tryCatch(select_region(read_sid(y)), error = function(err) {
      warning("SID ", y, " unavailable: ", conditionMessage(err),
              ". Shares were exported, but that year's SID audit is incomplete.")
      NULL
    })
  })
  downloaded <- years[!vapply(sid_results, is.null, logical(1))]
  write_table(tibble(SIDYear = years, Retrieved = years %in% downloaded), "SID_download_status.csv")
  annual_sid <- bind_rows(sid_results)
  if (nrow(annual_sid)) {
    annual_sid <- annual_sid |> mutate(In_fixed_matching_catalogue = fishstock %in% catalogue$fishstock)
    write_table(annual_sid, "SID_catalogue_by_year.csv")
    write_table(annual_sid |> filter(!is.na(SID_selection)) |>
      summarise(Regional_codes = n_distinct(fishstock),
        Codes_absent_from_fixed_catalogue = sum(!In_fixed_matching_catalogue), .by = SIDYear),
      "SID_regional_counts_by_year.csv")
    mismatches <- assigned |> filter(CLyear %in% downloaded) |>
      anti_join(transmute(annual_sid, CLyear = SIDYear, fishstock), by = c("CLyear", "fishstock")) |>
      summarise(Assigned_kg = sum(CLscientificWeight), .by = c(CLyear, fishstock))
    write_table(mismatches, "assigned_codes_absent_from_yearly_SID.csv")
  }
  writeLines(c(paste("Checkpoint:", normalizePath(checkpoint_file)),
    paste("Fixed matching SID year:", matching_reference_year),
    paste("Catch years:", paste(years, collapse = ", ")), paste("Ecoregion:", Ecoregion),
    paste("Rectangle vector:", normalizePath(rectangle_source)),
    paste("Legacy comparison performed:", comparison_available),
    "Stock assignments were NOT rematched using annual SID.",
    "No code renames, splits or mergers were automatically reconciled.",
    "Unresolved_candidate_kg is overlapping candidate exposure, not stock landings.",
    "Missing rectangles count outside, as in the original convention; their weight is reported.",
    "Sardine sensitivity reproduces the IFREMER rectangle convention; country condition is unverified.",
    paste("Additional share lookup:", lookup_filename),
    paste("Lookup base scenario:", paste(unique(share_lookup_export$Scenario), collapse = "; ")),
    paste("Pooled 9.a anchovy proxy enabled:", use_anchovy_9a_proxy),
    "Proxy eligibility: FAO ANE, CLarea 27.9.a or a child code; exclude unresolved taxonomy or missing status.",
    "Proxy uses each eligible landing aggregate once, including unresolved stock assignments.",
    "The proxy assumes the same ecoregion proportion for the two split stocks; it does not resolve their stock identity.",
    "Missing split-stock shares only are filled; existing stock-specific shares are preserved.",
    "Pooled kilograms are exported once per year in the proxy audit, not duplicated into stock totals.",
    "Share lookup percentages remain on the 0-100 scale. Divide by 100 when weighting assessment quantities.",
    "Existing stock-summary and coverage exports do not include the proxy; share_method identifies it in the new lookup.",
    "No raw RDBES data, original outputs or matching checkpoints were overwritten."),
    file.path(output_dir, "RUN_NOTES.txt"))
  if (nrow(coverage)) print(tibble::as_tibble(coverage), n = Inf, width = Inf)
  message("Results saved in: ", output_dir)
  invisible(list(output_directory = output_dir, coverage = coverage, comparison = comparison,
                 stock_summary = baseline, annual_sid = annual_sid,
                 share_lookup = share_lookup_export,
                 share_lookup_file = file.path(output_dir, lookup_filename),
                 anchovy_proxy_audit = anchovy_9a_pool))
}

rectangle_check_results <- run_rectangle_check()
