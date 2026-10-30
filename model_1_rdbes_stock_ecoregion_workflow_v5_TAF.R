# RDBES landings -> ICES stocks -> percentage in an ecoregion
# SID scope revision, 25 October 2026. Run the numbered sections in order.
# Weights are assumed to be kilograms. Start in a clean R session after saving
# dat to the six-column compact_file (see the accompanying recovery instructions).
# Existing data/ices_cache_simple downloads are reused with refresh_reference=FALSE.
# Shares describe the uniquely assigned RDBES subset, not all stock landings.
# Uses one reference stock list for all catch years; boundaries are current.
# Requires dplyr >= 1.1.1, tidyr, sf, httr, jsonlite, icesVocab and icesSD.
# install.packages(c("dplyr", "tidyr", "sf", "httr", "jsonlite", "icesVocab", "icesSD"),
#   repos = c(ICES = "https://ices-tools-prod.r-universe.dev", CRAN = "https://cloud.r-project.org"))

library(dplyr)
library(tidyr)
library(sf)

# 1. Settings ---------------------------------------------------------------
Ecoregion <- "Bay of Biscay and the Iberian Coast"
scope_ecoregions <- Ecoregion  # vector of regions to prepare; NULL keeps all taxa
stopifnot(is.null(scope_ecoregions) || Ecoregion %in% scope_ecoregions)
reference_year <- 2026L  # fixed SID stock catalogue used for matching
input_file <- "./data/extendedRDBES.RData"  # creates cl.all2
compact_file <- "./data/rdbes_landings_compact.rds"  # small, ALL areas
checkpoint_file <- "./data/rdbes_stock_matches.RData"  # reused by the rectangle-share workflow
rebuild_landings <- FALSE  # TRUE when the original RDBES extract changes
analysis_years <- NULL  # e.g. 2024L for a first run; NULL keeps all saved years
fu_file <- "./boot/initial/data/ICES_nephrops_FU/Nephrops Functional Units.shp"  # keep sidecar files together
output_dir <- file.path("./output/stock_matching_sid", gsub("[^A-Za-z0-9]", "_", Ecoregion))
cache_dir <- "./data/ices_cache_simple"
refresh_reference <- FALSE
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Small helpers ----------------------------------------------------------
clean <- function(x) na_if(trimws(as.character(x)), "")
# SID is a selection aid, not a stock boundary. Keep unknown tags visible.
# Documented Biscay additions (SID tags omit this region in the checked 2026 list):
# mon.27.78abd: https://asd.ices.dk/viewAdvice/7176
# rjb.27.89a: https://doi.org/10.17895/ices.advice.25019513
regional_stock_list <- function(stock_list, ecoregions) {
  if (is.null(ecoregions)) return(mutate(stock_list, SID_selection = "all stocks"))
  tags <- strsplit(coalesce(stock_list$SID_ecoregions, ""), ",", fixed = TRUE)
  if (!all(ecoregions %in% trimws(unlist(tags)))) stop("Unknown SID ecoregion name")
  additions <- if ("Bay of Biscay and the Iberian Coast" %in% ecoregions)
    c("mon.27.78abd", "rjb.27.89a") else character()
  stock_list |> mutate(SID_selection = case_when(
    vapply(tags, function(x) any(trimws(x) %in% ecoregions), logical(1)) ~ "SID ecoregion",
    fishstock %in% additions ~ "documented addition",
    is.na(SID_ecoregions) | SID_ecoregions == "" ~ "SID ecoregion missing",
    TRUE ~ NA_character_)) |> filter(!is.na(SID_selection))
}
cached <- function(name, expression) {
  path <- file.path(cache_dir, paste0(name, ".rds"))
  if (file.exists(path) && !refresh_reference) {
    message("  Reading cache: ", name)
    return(readRDS(path))
  }
  message("  Preparing: ", name)
  started <- proc.time()[[3]]
  x <- force(expression)
  saveRDS(x, path)
  message("  Finished ", name, " in ", round(proc.time()[[3]] - started, 1), " seconds")
  x
}
normalise <- function(x) st_make_valid(st_set_precision(st_transform(x, 3857), 1000))
empty <- st_sfc(st_geometrycollection(), crs = 3857)
union_geometry <- function(x) {
  g <- st_geometry(x)
  if (!length(g) || all(st_is_empty(g))) return(empty)
  if (length(g) == 1L) return(g)  # inputs already use the shared precision
  normalise(st_union(g))
}

# Query explicit IDs so a service limit cannot silently truncate the polygons.
read_gis <- function(layer, where = "1=1") {
  url <- paste0("https://gis.ices.dk/gis/rest/services/ICES_reference_layers/", layer, "/query")
  response <- httr::GET(url, query = list(where = where, returnIdsOnly = "true", f = "json"), httr::timeout(180))
  httr::stop_for_status(response)
  info <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))
  if (!is.null(info$error) || is.null(info$objectIdFieldName)) stop("GIS ID query failed: ", layer)
  ids <- info$objectIds
  if (!length(ids)) return(st_sf(geometry = st_sfc(crs = 3857)))
  message("    ", layer, ": ", length(ids), " features to download")
  pages <- lapply(split(ids, ceiling(seq_along(ids) / 100)), function(wanted) {
    path <- tempfile(fileext = ".geojson")
    on.exit(unlink(path))
    message("    Downloading ", length(wanted), " features; first ID ", wanted[1])
    r <- httr::GET(url, query = list(objectIds = paste(wanted, collapse = ","),
      outFields = "*", outSR = 3857, f = "geojson"), httr::timeout(180), httr::write_disk(path))
    httr::stop_for_status(r)
    message("    Download complete; reading geometry (", round(file.info(path)$size / 1024^2, 1), " MB)")
    x <- st_read(path, quiet = TRUE)
    stopifnot(nrow(x) == length(wanted), !is.na(st_crs(x)),
      setequal(as.character(x[[info$objectIdFieldName]]), as.character(wanted)))
    message("    Normalising geometry")
    normalise(x)
  })
  do.call(rbind, pages)
}

# Classify many supports against ONE unioned target, all already normalised
# to EPSG:3857 with 1 mm precision. Areas are numerical checks, not weights.
relation_vec <- function(support, target) {
  support <- sf::st_geometry(support)
  target <- sf::st_geometry(target)
  answer <- rep("unknown", length(support))
  if (!length(support) || !length(target) || any(sf::st_is_empty(target))) return(answer)
  stopifnot(length(target) == 1L)
  total <- as.numeric(sf::st_area(support))
  tolerance <- pmax(1, total * 1e-10)
  valid <- which(!sf::st_is_empty(support) & is.finite(total) & total > 0)
  answer[valid] <- "outside"
  # The original relation() tests inside <= tolerance first; preserve that
  # behaviour even for the unusual case of a whole support smaller than 1 m².
  valid <- valid[total[valid] > tolerance[valid]]
  if (!length(valid)) return(answer)
  hits <- valid[sf::st_intersects(target, support[valid])[[1L]]]
  if (!length(hits)) return(answer)
  full <- hits[sf::st_covers(target, support[hits])[[1L]]]
  answer[full] <- "full"
  boundary <- setdiff(hits, full)
  if (!length(boundary)) return(answer)
  # One indexed overlay for boundary cases, instead of one overlay per row.
  cut <- suppressWarnings(sf::st_intersection(support[boundary], target))
  inside <- numeric(length(boundary))
  if (length(cut)) inside[attr(cut, "idx")[, 1L]] <- as.numeric(sf::st_area(cut))
  answer[boundary] <- ifelse(inside <= tolerance[boundary], "outside",
    ifelse(total[boundary] - inside <= tolerance[boundary], "full", "partial"))
  answer
}


# 3. Aggregate once; reuse the small file on subsequent runs ----------------
prepare_landings <- function(input_file) {
  if (!rebuild_landings && exists("cl.all2", envir = .GlobalEnv, inherits = FALSE)) {
    message("  Reusing cl.all2 already in memory")
    source <- get("cl.all2", envir = .GlobalEnv)
  } else {
    message("  Loading the original RData once")
    loaded <- new.env(parent = emptyenv())
    load(input_file, envir = loaded)
    stopifnot(exists("cl.all2", envir = loaded, inherits = FALSE))
    source <- loaded$cl.all2
  }
  required <- c("CLyear", "CLarea", "CLstatisticalRectangle", "CLscientificWeight")
  stopifnot(all(required %in% names(source)), is.numeric(source$CLscientificWeight),
    any(c("CLspeciesCode", "CLspeciesFaoCode") %in% names(source)))
  source <- select(source, all_of(required), any_of(c("CLspeciesCode", "CLspeciesFaoCode", "CLcatchCategory")))
  if (exists("loaded", inherits = FALSE)) rm(loaded)
  if ("CLcatchCategory" %in% names(source)) {
    source <- filter(source, tolower(trimws(CLcatchCategory)) == "lan")
  } else message("  CLcatchCategory absent: input must already contain landings only")
  message("  Excluding ", sum(is.na(source$CLscientificWeight)), " rows with missing weight")
  source <- filter(source, !is.na(CLscientificWeight))
  stopifnot(nrow(source) > 0, all(is.finite(source$CLscientificWeight)), all(source$CLscientificWeight >= 0))
  for (nm in c("CLspeciesCode", "CLspeciesFaoCode")) if (!nm %in% names(source)) source[[nm]] <- NA_character_
  message("  Aggregating ", format(nrow(source), big.mark = ","), " landing rows; keeping ALL areas")
  result <- source |>
    mutate(CLyear = as.integer(as.character(CLyear)), CLarea = tolower(clean(CLarea)),
      CLstatisticalRectangle = toupper(clean(CLstatisticalRectangle)),
      supplied_aphia = clean(CLspeciesCode), CLspeciesFaoCode = toupper(clean(CLspeciesFaoCode))) |>
    summarise(CLscientificWeight = sum(CLscientificWeight),
      .by = c(CLyear, CLarea, CLstatisticalRectangle, supplied_aphia, CLspeciesFaoCode))
  stopifnot(!anyNA(result$CLyear), isTRUE(all.equal(sum(result$CLscientificWeight), sum(source$CLscientificWeight))))
  result
}
message("[3/10] Reading or preparing the compact landings table")
if (file.exists(compact_file) && !rebuild_landings) {
  dat <- readRDS(compact_file)
} else {
  dat <- prepare_landings(input_file)
  dir.create(dirname(compact_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(dat, compact_file, compress = FALSE)
}
# A recovered compact file must contain the same six aggregate columns.
dat <- select(dat, CLyear, CLarea, CLstatisticalRectangle, supplied_aphia, CLspeciesFaoCode, CLscientificWeight)
if (!is.null(analysis_years)) dat <- filter(dat, CLyear %in% analysis_years)
stopifnot(nrow(dat) > 0, !anyNA(dat$CLyear), !anyNA(dat$CLscientificWeight),
  all(is.finite(dat$CLscientificWeight)), all(dat$CLscientificWeight >= 0))
dat$row_id <- seq_len(nrow(dat))
input_kg <- sum(dat$CLscientificWeight)
analysed_years <- sort(unique(dat$CLyear))
message("  Using ", format(nrow(dat), big.mark = ","), " aggregate rows (",
  format(object.size(dat), units = "auto"), ")")

# 4. ICES stock, species and area relationships -----------------------------
message("[4/10] Reading stock and species references")
stock_list <- cached(paste0("stocks_", reference_year), icesSD::getSD(year = reference_year)) |>
  transmute(fishstock = as.character(StockKeyLabel), description = StockKeyDescription,
            SID_ecoregions = clean(EcoRegion)) |>
  distinct(fishstock, .keep_all = TRUE)
stopifnot(nrow(stock_list) > 0)
stock_species <- cached("stock_species", icesVocab::getCodeTypeRelation("ICES_StockCode", "SpecWoRMS")) |>
  transmute(fishstock = as.character(ICES_StockCode), aphia = as.character(SpecWoRMS)) |>
  semi_join(stock_list, by = "fishstock") |> distinct()
stock_areas <- cached("stock_areas", icesVocab::getCodeTypeRelation("ICES_StockCode", "ICES_Area")) |>
  transmute(fishstock = as.character(ICES_StockCode), area = tolower(clean(ICES_Area))) |>
  filter(!is.na(area)) |> semi_join(stock_list, by = "fishstock") |> distinct()
fao_lookup <- cached("asfis_worms", icesVocab::getCodeTypeRelation("SpecASFIS", "SpecWoRMS")) |>
  transmute(CLspeciesFaoCode = as.character(SpecASFIS), fallback = as.character(SpecWoRMS)) |>
  distinct() |> summarise(fallback = if (n_distinct(fallback) == 1L) first(fallback) else NA_character_,
    .by = CLspeciesFaoCode)
dat <- dat |> left_join(fao_lookup, by = "CLspeciesFaoCode", relationship = "many-to-one") |>
  mutate(taxon_conflict = !is.na(supplied_aphia) & !is.na(fallback) & supplied_aphia != fallback,
    aphia = if_else(taxon_conflict, NA_character_, coalesce(supplied_aphia, fallback)),
    aphia = if_else(!is.na(aphia) & grepl("^[0-9]+$", aphia), aphia, NA_character_))
# A genus/group code is not expanded into guessed species memberships.

# Scope by target taxa BEFORE any polygons, retaining ALL areas and competitors.
scope_stocks <- regional_stock_list(stock_list, scope_ecoregions)
selected_taxa <- unique(stock_species$aphia[stock_species$fishstock %in% scope_stocks$fishstock])
dat <- dat |> mutate(scope_status = case_when(
  is.null(scope_ecoregions) ~ "all_taxa_retained",
  is.na(aphia) ~ "unresolved_taxon_retained",
  aphia %in% selected_taxa ~ "selected_taxon_retained",
  TRUE ~ "other_taxon_excluded"))
scope_summary <- dat |> summarise(Rows = n(), Landings_kg = sum(CLscientificWeight),
  .by = c(CLyear, scope_status, CLspeciesFaoCode, aphia))
stopifnot(isTRUE(all.equal(sum(scope_summary$Landings_kg), input_kg)))
write.csv(scope_summary, file.path(output_dir, "landings_scope.csv"), row.names = FALSE)
dat <- filter(dat, scope_status != "other_taxon_excluded") |> select(-scope_status)
if (!nrow(dat)) stop("No landings remain in the selected taxon scope; see landings_scope.csv")
input_kg <- sum(dat$CLscientificWeight)
matching_scope <- list(ecoregions = scope_ecoregions, years = analysed_years,
  included_aphia = if (is.null(scope_ecoregions)) NULL else selected_taxa)
message("  SID scope: ", sum(scope_stocks$SID_selection == "SID ecoregion"),
  " tagged stocks, ", sum(scope_stocks$SID_selection == "documented addition"),
  " documented additions, ", sum(scope_stocks$SID_selection == "SID ecoregion missing"),
  " stocks with missing region tags")
message("  Retained ", nrow(dat), " aggregate rows across ALL areas; ",
  n_distinct(stock_species$fishstock[stock_species$aphia %in% dat$aphia]),
  " stocks to check, including same-taxon alternatives")

# 5. Official area, ecoregion and rectangle polygons ------------------------
message("[5/10] Reading ICES polygons")
message("  5a: ICES areas")
areas <- cached("areas", read_gis("ICES_Areas/MapServer/0"))
eco_name <- recode(Ecoregion, "Norwegian Sea" = "Norwegean Sea",
  "Oceanic Northeast Atlantic" = "Oceanic north-east Atlantic")
message("  5b: Target ecoregion")
ecoregion <- cached(paste0("ecoregion_", gsub("[^A-Za-z0-9]", "_", eco_name)),
  read_gis("ICES_Ecoregions/MapServer/0", paste0("Ecoregion='", gsub("'", "''", eco_name), "'")))
stopifnot(nrow(ecoregion) > 0)
fu25_rectangles <- c("15E0", "15E1", "16E0", "16E1", "17E1")
rectangle_ids <- sort(unique(na.omit(c(dat$CLstatisticalRectangle, fu25_rectangles))))
# Reuse polygons from existing batches even when scoping changes the batch IDs.
rectangle_cache_files <- list.files(cache_dir, "^rectangles_[0-9]+[.]rds$", full.names = TRUE)
rectangle_pool <- if (!refresh_reference && length(rectangle_cache_files))
  do.call(rbind, lapply(rectangle_cache_files, function(path) readRDS(path)$data)) else NULL
if (!is.null(rectangle_pool)) rectangle_pool <- rectangle_pool[!duplicated(rectangle_pool$ICESNAME), ]
rectangle_batches <- split(rectangle_ids, ceiling(seq_along(rectangle_ids) / 100))
message("  5c: ", length(rectangle_ids), " rectangles in ", length(rectangle_batches), " batches")
rectangles <- do.call(rbind, lapply(seq_along(rectangle_batches), function(i) {
  message("    Rectangle batch ", i, "/", length(rectangle_batches))
  ids <- rectangle_batches[[i]]
  path <- file.path(cache_dir, paste0("rectangles_", i, ".rds"))
  previous <- if (file.exists(path) && !refresh_reference) readRDS(path) else NULL
  if (!is.null(previous) && identical(previous$ids, ids)) return(previous$data)
  if (!is.null(rectangle_pool) && all(ids %in% rectangle_pool$ICESNAME))
    return(rectangle_pool[match(ids, rectangle_pool$ICESNAME), ])
  x <- read_gis("ICES_Statrec_complete/MapServer/0",
    paste0("ICESNAME IN ('", paste(gsub("'", "''", ids), collapse = "','"), "')"))
  saveRDS(list(ids = ids, data = x), path)
  x
}))
rm(rectangle_pool)
stopifnot(!anyDuplicated(rectangles$ICESNAME), all(fu25_rectangles %in% rectangles$ICESNAME))

# 6. Nephrops FU -> stock polygons ------------------------------------------
message("[6/10] Preparing Nephrops functional units")
single_fu <- c(5:17, 19, 22, 25, 30:34)
fu_key <- bind_rows(tibble(fu_no = paste0("FU", single_fu), fishstock = paste0("nep.fu.", single_fu)),
  tribble(~fu_no, ~fishstock,
    "4outFU", "nep.27.4outFU", "6aoutFU", "nep.27.6aoutFU", "7outFU", "nep.27.7outFU",
    "FU3", "nep.fu.3-4", "FU4", "nep.fu.3-4", "FU2021", "nep.fu.2021",
    "FU23", "nep.fu.2324", "FU24", "nep.fu.2324", "FU26", "nep.fu.2627", "FU27", "nep.fu.2627",
    "FU28", "nep.fu.2829", "FU29", "nep.fu.2829"))
nephrops_fu <- normalise(st_read(fu_file, quiet = TRUE)) |> select(fu_no)
stopifnot(!anyDuplicated(nephrops_fu$fu_no), !any(st_is_empty(nephrops_fu)),
  sum(nephrops_fu$fu_no == "FU25") == 1L)
# The supplied GIS FU25 omits the 2021 extension; use the May 2024 stock annex.
# https://ices-library.figshare.com/ndownloader/files/48791170
st_geometry(nephrops_fu)[nephrops_fu$fu_no == "FU25"] <-
  union_geometry(rectangles[match(fu25_rectangles, rectangles$ICESNAME), ])
nephrops_fu <- left_join(nephrops_fu, fu_key, by = "fu_no", relationship = "many-to-one")
nephrops_unmapped <- filter(nephrops_fu, is.na(fishstock)) |> st_drop_geometry() # FU1/FU18
nephrops <- nephrops_fu |> filter(!is.na(fishstock)) |>
  group_by(fishstock) |> summarise(do_union = TRUE, .groups = "drop")
nephrops$complete <- vapply(nephrops$fishstock, function(s)
  all(fu_key$fu_no[fu_key$fishstock == s] %in% nephrops_fu$fu_no), logical(1))

# 7. Build stock polygons and the reported geographic support ---------------
message("[7/10] Building geographic supports and stock boundaries")
used_stocks <- stock_species |> filter(aphia %in% dat$aphia) |> distinct(fishstock)
ordinary <- used_stocks$fishstock[!grepl("^(nep[.]|san[.]sa[.])", used_stocks$fishstock)]
has_sandeel <- any(startsWith(used_stocks$fishstock, "san.sa."))
area_codes <- unique(na.omit(c(dat$CLarea, stock_areas$area[stock_areas$fishstock %in% ordinary],
  if (has_sandeel) c("27.3", "27.4"))))
message("  Preparing ", length(area_codes), " area codes; sharing identical polygon unions")
atomic_parts <- setNames(lapply(area_codes, function(a) which(
  areas$Area_Full == a | startsWith(areas$Area_Full, paste0(a, ".")))), area_codes)
geometry_cache <- new.env(parent = emptyenv())
union_parts <- function(i) {
  i <- sort(unique(i))
  if (!length(i)) return(empty)
  key <- paste(i, collapse = ",")
  if (!exists(key, envir = geometry_cache, inherits = FALSE)) {
    message("    Combining ", length(i), " area polygons")
    assign(key, union_geometry(areas[i, ]), envir = geometry_cache)
  }
  get(key, envir = geometry_cache, inherits = FALSE)
}
area_polygons <- lapply(atomic_parts, union_parts)
san_envelope <- if (has_sandeel) union_geometry(c(area_polygons[["27.3"]], area_polygons[["27.4"]])) else empty
stock_polygons <- setNames(vector("list", nrow(used_stocks)), used_stocks$fishstock)
stock_complete <- setNames(rep(FALSE, nrow(used_stocks)), used_stocks$fishstock)
for (s in used_stocks$fishstock) {
  message("  Stock boundary: ", s)
  if (startsWith(s, "nep.")) {
    k <- match(s, nephrops$fishstock)
    stock_polygons[[s]] <- if (is.na(k)) empty else st_geometry(nephrops)[k]
    stock_complete[s] <- !is.na(k) && nephrops$complete[k]
    next
  }
  if (startsWith(s, "san.sa.")) { # broad candidate envelope, NOT a stock boundary
    stock_polygons[[s]] <- san_envelope
    stock_complete[s] <- TRUE
    next
  }
  codes <- stock_areas$area[stock_areas$fishstock == s]
  parts <- area_polygons[codes]
  stock_polygons[[s]] <- union_parts(unlist(atomic_parts[codes], use.names = FALSE))
  stock_complete[s] <- length(parts) > 0 && all(vapply(parts, function(g) !st_is_empty(g), logical(1)))
}
geography <- dat |> distinct(CLarea, CLstatisticalRectangle) |> mutate(geo_id = row_number())
support <- rep(list(empty), nrow(geography))
geography$geography_status <- "unknown"
area_groups <- split(seq_len(nrow(geography)), match(geography$CLarea, unique(geography$CLarea)))
# Process all rectangles in each area together; never dissolve all marine areas.
message("  Building ", nrow(geography), " geographic supports in ", length(area_groups), " area batches")
for (ids in area_groups) {
  area_name <- geography$CLarea[ids[1]]
  message("    Support batch: ", area_name, " (", length(ids), " locations)")
  a <- area_polygons[[area_name]]
  if (is.null(a)) a <- empty
  r <- match(geography$CLstatisticalRectangle[ids], rectangles$ICESNAME)
  if (!st_is_empty(a)) {
    support[ids[is.na(r)]] <- rep(list(a), sum(is.na(r)))
    geography$geography_status[ids] <- ifelse(is.na(r), "area_only", "area_rectangle_conflict")
  }
  has_r <- which(!is.na(r))
  if (!length(has_r)) next
  # For rectangle-only data, intersect with the area pieces before unioning.
  target <- if (!st_is_empty(a)) a else st_geometry(areas)
  cut <- suppressWarnings(st_intersection(st_geometry(rectangles)[r[has_r]], target))
  if (!length(cut)) next
  pieces <- split(seq_along(cut), attr(cut, "idx")[, 1L])
  for (key in names(pieces)) {
    g <- union_geometry(cut[pieces[[key]]])
    if (sum(as.numeric(st_area(g))) <= 0) next
    i <- ids[has_r[as.integer(key)]]
    support[[i]] <- g
    geography$geography_status[i] <- if (!st_is_empty(a)) "area_and_rectangle" else "rectangle_only"
  }
}
support_sfc <- do.call(c, support)

# 8. Match each species/geography to candidate stocks -----------------------
message("[8/10] Classifying geographic supports by stock")
# Compact, conservative review list. These stocks need additional boundaries or
# stock-origin rules; this version leaves them unresolved even in some areas
# where a more detailed rule could work. The list is not exhaustive.
# Source descriptions: https://stockdatabase.ices.dk/Manage/rollingissues.aspx
needs_review <- c("ane.27.9aS", "ane.27.9aW", "cap.27.1-2", "cap.27.2a514", "pra.27.4a", "pra.27.3a4a",
  "her.27.nirs", "her.27.irls", "her.27.6aN", "her.27.6aS7bc", "her.27.20-24", "her.27.3a47d",
  "her.27.1-24a514a", "her.27.5a", "cod.27.1-2", "cod.27.1-2.coastN", "cod.27.2.coastS",
  "cod.27.22-24", "cod.27.24-32", "cod.27.21", "cod.27.46a7d20", "reb.2127.dp", "reb.2127.sp",
  "reb.27.5a14", "reb.27.14b", "seh.27.1", "seh.27.125a14", "cyo.27.nea", "guq.27.nea")
dat <- left_join(dat, geography, by = c("CLarea", "CLstatisticalRectangle"), relationship = "many-to-one")
rules <- dat |> distinct(geo_id, aphia, CLarea) |> filter(!is.na(aphia)) |>
  inner_join(stock_species, by = "aphia", relationship = "many-to-many")
rules$evidence <- rep("unknown", nrow(rules))
for (s in unique(rules$fishstock)) {
  rows <- which(rules$fishstock == s)
  ids <- unique(rules$geo_id[rows])
  message("  Matching ", s, ": ", length(ids), " distinct locations")
  z <- relation_vec(support_sfc[ids], stock_polygons[[s]])
  rules$evidence[rows] <- z[match(rules$geo_id[rows], ids)]
}
rules$evidence <- vapply(seq_len(nrow(rules)), function(i) {
  s <- rules$fishstock[i]; a <- rules$CLarea[i]
  valid_a <- !is.na(a) && !is.null(area_polygons[[a]]) && !st_is_empty(area_polygons[[a]])
  if (valid_a && (a == "27" || startsWith(a, "27.")) && s %in% c("cod.21.1.isc", "cod.21.1.osc")) return("outside")
  if (valid_a && startsWith(a, "27.") && s == "cod.21.27.1.14" &&
      !(a == "27.14" || startsWith(a, "27.14."))) return("outside")
  z <- rules$evidence[i]
  if (z == "outside" && !stock_complete[s]) z <- "unknown"
  if (z != "outside" && (s %in% needs_review || startsWith(s, "san.sa."))) z <- "review"
  z
}, character(1))
candidates <- dat |> select(row_id, geo_id, aphia) |>
  inner_join(filter(rules, evidence != "outside") |> select(-CLarea),
    by = c("geo_id", "aphia"), relationship = "many-to-many")
matches <- candidates |> summarise(n_candidates = n(),
  candidate_stocks = paste(fishstock, collapse = ";"),
  candidate_evidence = paste(paste(fishstock, evidence, sep = ":"), collapse = ";"),
  fishstock = if (n() == 1L && first(evidence) == "full") first(fishstock) else NA_character_, .by = row_id)
dat <- left_join(dat, matches, by = "row_id", relationship = "one-to-one") |>
  mutate(n_candidates = coalesce(n_candidates, 0L), match_status = case_when(
    is.na(aphia) ~ "unresolved_taxon", !is.na(fishstock) ~ "unique",
    n_candidates == 0L ~ "unmatched", TRUE ~ "ambiguous_or_needs_review"))

# 9. Ecoregion classification and a short audit ------------------------------
message("[9/10] Classifying ecoregion and auditing the results")
# Boundary also includes coastline differences between reference polygons.
eco_relation <- relation_vec(support_sfc, union_geometry(ecoregion))
geography$CS <- recode(eco_relation, full = "inside", outside = "outside", partial = "boundary", unknown = "unknown")
dat <- left_join(dat, select(geography, geo_id, CS), by = "geo_id", relationship = "many-to-one")

audit_landings <- function(dat, candidates, stock_list, input_kg, years = sort(unique(dat$CLyear))) {
  stopifnot(!anyDuplicated(dat$row_id), !anyDuplicated(candidates[c("row_id", "fishstock")]),
    isTRUE(all.equal(sum(dat$CLscientificWeight), input_kg)))
  totals <- dat |> filter(match_status == "unique") |>
    summarise(Assigned_kg = sum(CLscientificWeight),
      Inside_kg = sum(CLscientificWeight[CS == "inside"]),
      Outside_kg = sum(CLscientificWeight[CS == "outside"]),
      Uncertain_kg = sum(CLscientificWeight[CS %in% c("boundary", "unknown")]),
      .by = c(CLyear, fishstock))
  exposure <- candidates |> select(row_id, fishstock) |>
    inner_join(dat |> filter(match_status != "unique") |> select(row_id, CLyear, CLscientificWeight),
      by = "row_id", relationship = "many-to-many") |>
    summarise(Unresolved_candidate_kg = sum(CLscientificWeight), .by = c(CLyear, fishstock))
  percent <- function(x, total) ifelse(total > 0, 100 * x / total, NA_real_)
  result <- crossing(CLyear = years, fishstock = stock_list$fishstock) |>
    left_join(totals, by = c("CLyear", "fishstock")) |>
    left_join(exposure, by = c("CLyear", "fishstock")) |>
    mutate(across(ends_with("_kg"), ~ coalesce(.x, 0)),
      ShareLandingsInEcoregion = percent(Inside_kg, Inside_kg + Outside_kg),
      AssignedShareLower_pct = percent(Inside_kg, Assigned_kg),
      AssignedShareUpper_pct = percent(Inside_kg + Uncertain_kg, Assigned_kg),
      SpatialCoverage_pct = percent(Inside_kg + Outside_kg, Assigned_kg),
      audit_status = case_when(Unresolved_candidate_kg > 0 ~ "unresolved_stock_assignment",
        Assigned_kg == 0 ~ "no_assigned_landings", Uncertain_kg > 0 ~ "uncertain_ecoregion",
        TRUE ~ "resolved_assigned_subset")) |>
    left_join(stock_list, by = "fishstock")
  stopifnot(isTRUE(all.equal(sum(result$Assigned_kg), sum(dat$CLscientificWeight[dat$match_status == "unique"]))))
  result
}
stock_summary <- audit_landings(dat, candidates, stock_list, input_kg, years = analysed_years) |>
  mutate(ReferenceYear = reference_year, Ecoregion = Ecoregion,
    Special_rule_required = fishstock %in% needs_review | startsWith(fishstock, "san.sa."),
    Missing_species_reference = !fishstock %in% stock_species$fishstock,
    Missing_area_reference = !fishstock %in% stock_areas$fishstock & !fishstock %in% nephrops$fishstock)
regional_stocks <- regional_stock_list(stock_list, Ecoregion)
stock_summary_ecoregion <- inner_join(stock_summary,
  select(regional_stocks, fishstock, SID_selection), by = "fishstock", relationship = "many-to-one")
matching_summary <- dat |> mutate(match_status = case_when(
  match_status == "unique" & fishstock %in% regional_stocks$fishstock ~ "unique_target_stock",
  match_status == "unique" ~ "unique_other_stock", TRUE ~ match_status)) |>
  summarise(Landings_kg = sum(CLscientificWeight), .by = c(CLyear, match_status))
unresolved_landings <- dat |> filter(match_status != "unique")

# 10. Save diagnostics and completed stock-matching checkpoint -------------
message("[10/10] Saving results and stock-matching checkpoint")

# This polygon-based stock summary is retained as a diagnostic only.
# The production ecoregion share is calculated later by the rectangle workflow.
write.csv(stock_summary_ecoregion,
          file.path(output_dir, "stock_summary_polygon_diagnostic.csv"),
          row.names = FALSE)
write.csv(matching_summary, file.path(output_dir, "matching_summary.csv"), row.names = FALSE)
write.csv(unresolved_landings, file.path(output_dir, "unresolved_landings.csv"), row.names = FALSE)

# Save exactly the stock-matching objects needed by the downstream rectangle
# workflow. The checkpoint contains the current ecoregion scope in matching_scope.
dir.create(dirname(checkpoint_file), recursive = TRUE, showWarnings = FALSE)
save(
  dat, candidates, stock_list, reference_year, matching_scope, input_kg,
  file = checkpoint_file,
  compress = FALSE
)

message(
  "Saved stock-matching checkpoint: ", checkpoint_file,
  "; years: ", paste(sort(unique(dat$CLyear)), collapse = ", "),
  "; reference year: ", reference_year
)

print(matching_summary)
print(stock_summary |> filter(CLyear == 2025, fishstock %in% c("mac.27.nea", "whb.27.1-91214", "nep.fu.25", "nep.fu.31")) |>
  select(fishstock, Assigned_kg, ShareLandingsInEcoregion, Unresolved_candidate_kg, audit_status))
# Inspect dat, candidates, rules and geography in R for individual records.
# Unresolved_candidate_kg is overlapping exposure: NEVER sum it across stocks.
# Shares and bounds apply only to uniquely assigned weight. Bounds do not cover
# records whose stock identity is unresolved or landings absent from RDBES.
# Optional diagnostic: icesVMS::get_csquare(stat_rec = "04D6") |> count(ecoregion)
# Those counts describe returned C-squares; they are not landings percentages.
# stock_summary remains the full internal audit in memory.
# stock_summary_polygon_diagnostic.csv is diagnostic only; do not use it as the production regional share.
# The downstream rectangle workflow reads checkpoint_file and calculates the production share.
# matching_scope records that this checkpoint was prepared for the selected ecoregion.
