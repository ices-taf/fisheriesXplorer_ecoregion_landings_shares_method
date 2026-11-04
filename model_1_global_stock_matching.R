# Global RDBES landings -> ICES stock matching
#
# Purpose:
#   Match RDBES landings to the fixed SID stock catalogue once, across all
#   FisheriesXplorer ecoregions. Regional shares are calculated later by
#   model_2_ecoregion_shares.R using the statistical-rectangle convention.
#
# Important:
#   - Uses one fixed SID stock catalogue for all RDBES years.
#   - Retains all areas for each selected taxon so stock denominators are global.
#   - Same-taxon competing stocks are retained during matching.
#   - The sardine 25E4/25E5 reassignment is an accepted production rule.
#   - No ecoregion share is calculated in this script.

library(dplyr)
library(tidyr)
library(sf)

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
if (!exists("reference_year", inherits = FALSE)) reference_year <- 2026L
if (!exists("analysis_years", inherits = FALSE)) analysis_years <- NULL
if (!exists("rebuild_landings", inherits = FALSE)) rebuild_landings <- FALSE
if (!exists("refresh_reference", inherits = FALSE)) refresh_reference <- FALSE

input_file <- "./data/extendedRDBES.RData"
compact_file <- "./data/rdbes_landings_compact.rds"
checkpoint_file <- "./model/rdbes_stock_matches.RData"
fu_file <- "./boot/initial/data/ICES_nephrops_FU/Nephrops Functional Units.shp"
cache_dir <- "./data/ices_cache_simple"
diagnostic_dir <- "./model/stock_matching_global"

dir.create("./model", recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diagnostic_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Helpers ----------------------------------------------------------------
clean <- function(x) na_if(trimws(as.character(x)), "")

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
  if (length(g) == 1L) return(g)
  normalise(st_union(g))
}

# Query explicit feature IDs so service limits cannot silently truncate layers.
read_gis <- function(layer, where = "1=1") {
  url <- paste0(
    "https://gis.ices.dk/gis/rest/services/ICES_reference_layers/",
    layer, "/query"
  )
  response <- httr::GET(
    url,
    query = list(where = where, returnIdsOnly = "true", f = "json"),
    httr::timeout(180)
  )
  httr::stop_for_status(response)
  info <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))
  if (!is.null(info$error) || is.null(info$objectIdFieldName)) {
    stop("GIS ID query failed: ", layer)
  }
  ids <- info$objectIds
  if (!length(ids)) return(st_sf(geometry = st_sfc(crs = 3857)))

  message("    ", layer, ": ", length(ids), " features to download")
  pages <- lapply(split(ids, ceiling(seq_along(ids) / 100)), function(wanted) {
    path <- tempfile(fileext = ".geojson")
    on.exit(unlink(path))
    r <- httr::GET(
      url,
      query = list(
        objectIds = paste(wanted, collapse = ","),
        outFields = "*", outSR = 3857, f = "geojson"
      ),
      httr::timeout(180),
      httr::write_disk(path)
    )
    httr::stop_for_status(r)
    x <- st_read(path, quiet = TRUE)
    stopifnot(
      nrow(x) == length(wanted),
      !is.na(st_crs(x)),
      setequal(
        as.character(x[[info$objectIdFieldName]]),
        as.character(wanted)
      )
    )
    normalise(x)
  })
  do.call(rbind, pages)
}

# Classify reported geographical supports against one stock boundary.
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
  valid <- valid[total[valid] > tolerance[valid]]
  if (!length(valid)) return(answer)

  hits <- valid[sf::st_intersects(target, support[valid])[[1L]]]
  if (!length(hits)) return(answer)
  full <- hits[sf::st_covers(target, support[hits])[[1L]]]
  answer[full] <- "full"

  boundary <- setdiff(hits, full)
  if (!length(boundary)) return(answer)
  cut <- suppressWarnings(sf::st_intersection(support[boundary], target))
  inside <- numeric(length(boundary))
  if (length(cut)) inside[attr(cut, "idx")[, 1L]] <- as.numeric(sf::st_area(cut))
  answer[boundary] <- ifelse(
    inside <= tolerance[boundary], "outside",
    ifelse(total[boundary] - inside <= tolerance[boundary], "full", "partial")
  )
  answer
}

split_sid_ecoregions <- function(x) {
  lapply(strsplit(coalesce(x, ""), ",", fixed = TRUE), trimws)
}

# 3. Aggregate landings once -------------------------------------------------
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
  stopifnot(
    all(required %in% names(source)),
    is.numeric(source$CLscientificWeight),
    any(c("CLspeciesCode", "CLspeciesFaoCode") %in% names(source))
  )

  source <- select(
    source,
    all_of(required),
    any_of(c("CLspeciesCode", "CLspeciesFaoCode", "CLcatchCategory"))
  )
  if (exists("loaded", inherits = FALSE)) rm(loaded)

  if ("CLcatchCategory" %in% names(source)) {
    source <- filter(source, tolower(trimws(CLcatchCategory)) == "lan")
  } else {
    message("  CLcatchCategory absent: input must already contain landings only")
  }

  source <- filter(source, !is.na(CLscientificWeight))
  stopifnot(
    nrow(source) > 0,
    all(is.finite(source$CLscientificWeight)),
    all(source$CLscientificWeight >= 0)
  )

  for (nm in c("CLspeciesCode", "CLspeciesFaoCode")) {
    if (!nm %in% names(source)) source[[nm]] <- NA_character_
  }

  message("  Aggregating ", format(nrow(source), big.mark = ","), " landing rows; keeping ALL areas")
  result <- source |>
    mutate(
      CLyear = as.integer(as.character(CLyear)),
      CLarea = tolower(clean(CLarea)),
      CLstatisticalRectangle = toupper(clean(CLstatisticalRectangle)),
      supplied_aphia = clean(CLspeciesCode),
      CLspeciesFaoCode = toupper(clean(CLspeciesFaoCode))
    ) |>
    summarise(
      CLscientificWeight = sum(CLscientificWeight),
      .by = c(
        CLyear, CLarea, CLstatisticalRectangle,
        supplied_aphia, CLspeciesFaoCode
      )
    )

  stopifnot(
    !anyNA(result$CLyear),
    isTRUE(all.equal(
      sum(result$CLscientificWeight),
      sum(source$CLscientificWeight)
    ))
  )
  result
}

message("[1/8] Reading or preparing compact RDBES landings")
compact_is_current <- file.exists(compact_file) &&
  file.exists(input_file) &&
  file.info(compact_file)$mtime >= file.info(input_file)$mtime

if (compact_is_current && !rebuild_landings) {
  dat <- readRDS(compact_file)
} else {
  if (file.exists(compact_file) && !rebuild_landings) {
    message("  extendedRDBES.RData is newer than the compact file; rebuilding compact landings")
  }
  dat <- prepare_landings(input_file)
  saveRDS(dat, compact_file, compress = FALSE)
}

dat <- select(
  dat,
  CLyear, CLarea, CLstatisticalRectangle,
  supplied_aphia, CLspeciesFaoCode, CLscientificWeight
)
if (!is.null(analysis_years)) dat <- filter(dat, CLyear %in% analysis_years)
stopifnot(
  nrow(dat) > 0,
  !anyNA(dat$CLyear),
  !anyNA(dat$CLscientificWeight),
  all(is.finite(dat$CLscientificWeight)),
  all(dat$CLscientificWeight >= 0)
)
dat$row_id <- seq_len(nrow(dat))
input_kg_before_scope <- sum(dat$CLscientificWeight)
analysed_years <- sort(unique(dat$CLyear))

# 4. Fixed SID catalogue and ICES Vocabulary relationships ------------------
message("[2/8] Reading fixed SID catalogue and ICES Vocabulary relationships")
stock_list <- cached(
  paste0("stocks_", reference_year),
  icesSD::getSD(year = reference_year)
) |>
  transmute(
    fishstock = as.character(StockKeyLabel),
    description = as.character(StockKeyDescription),
    SID_ecoregions = clean(EcoRegion)
  ) |>
  distinct(fishstock, .keep_all = TRUE)
stopifnot(nrow(stock_list) > 0L, !anyDuplicated(stock_list$fishstock))

stock_species <- cached(
  "stock_species",
  icesVocab::getCodeTypeRelation("ICES_StockCode", "SpecWoRMS")
) |>
  transmute(
    fishstock = as.character(ICES_StockCode),
    aphia = as.character(SpecWoRMS)
  ) |>
  semi_join(stock_list, by = "fishstock") |>
  distinct()

stock_areas <- cached(
  "stock_areas",
  icesVocab::getCodeTypeRelation("ICES_StockCode", "ICES_Area")
) |>
  transmute(
    fishstock = as.character(ICES_StockCode),
    area = tolower(clean(ICES_Area))
  ) |>
  filter(!is.na(area)) |>
  semi_join(stock_list, by = "fishstock") |>
  distinct()

fao_lookup <- cached(
  "asfis_worms",
  icesVocab::getCodeTypeRelation("SpecASFIS", "SpecWoRMS")
) |>
  transmute(
    CLspeciesFaoCode = as.character(SpecASFIS),
    fallback = as.character(SpecWoRMS)
  ) |>
  distinct() |>
  summarise(
    fallback = if (n_distinct(fallback) == 1L) first(fallback) else NA_character_,
    .by = CLspeciesFaoCode
  )

dat <- dat |>
  left_join(fao_lookup, by = "CLspeciesFaoCode", relationship = "many-to-one") |>
  mutate(
    taxon_conflict = !is.na(supplied_aphia) & !is.na(fallback) & supplied_aphia != fallback,
    aphia = if_else(taxon_conflict, NA_character_, coalesce(supplied_aphia, fallback)),
    aphia = if_else(!is.na(aphia) & grepl("^[0-9]+$", aphia), aphia, NA_character_)
  )

# Scope the expensive matching to taxa belonging to stocks in the 11 supported
# ecoregions, while keeping stocks with missing SID region tags visible and
# retaining same-taxon competitors later in the matching step.
sid_tags <- split_sid_ecoregions(stock_list$SID_ecoregions)
documented_scope_additions <- c("mon.27.78abd", "rjb.27.89a")
stock_list <- stock_list |>
  mutate(
    In_supported_ecoregions = vapply(
      sid_tags,
      function(z) any(z %in% supported_ecoregions),
      logical(1)
    ),
    SID_region_missing = is.na(SID_ecoregions) | SID_ecoregions == "",
    In_matching_scope = In_supported_ecoregions |
      SID_region_missing |
      fishstock %in% documented_scope_additions
  )

scope_stocks <- filter(stock_list, In_matching_scope)
selected_taxa <- unique(stock_species$aphia[stock_species$fishstock %in% scope_stocks$fishstock])

dat <- dat |>
  mutate(
    scope_status = case_when(
      is.na(aphia) ~ "unresolved_taxon_retained",
      aphia %in% selected_taxa ~ "selected_taxon_retained",
      TRUE ~ "other_taxon_excluded"
    )
  )

scope_summary <- dat |>
  summarise(
    Rows = n(),
    Landings_kg = sum(CLscientificWeight),
    .by = c(CLyear, scope_status, CLspeciesFaoCode, aphia)
  )
stopifnot(isTRUE(all.equal(sum(scope_summary$Landings_kg), input_kg_before_scope)))
write.csv(scope_summary, file.path(diagnostic_dir, "landings_scope.csv"), row.names = FALSE)
write.csv(scope_stocks, file.path(diagnostic_dir, "stock_scope_catalogue.csv"), row.names = FALSE)

dat <- filter(dat, scope_status != "other_taxon_excluded") |> select(-scope_status)
if (!nrow(dat)) stop("No RDBES landings remain in the supported ecoregion taxon scope.")
input_kg <- sum(dat$CLscientificWeight)
matching_scope <- list(
  ecoregions = supported_ecoregions,
  years = analysed_years,
  included_aphia = selected_taxa
)

message(
  "  Fixed SID ", reference_year, ": ", nrow(stock_list), " active stocks; ",
  nrow(scope_stocks), " stocks retained in matching scope"
)
message(
  "  Retained ", nrow(dat), " aggregate RDBES rows across ALL reported areas"
)

# 5. ICES areas and rectangles ----------------------------------------------
message("[3/8] Reading ICES area and statistical-rectangle polygons")
areas <- cached("areas", read_gis("ICES_Areas/MapServer/0"))

# FU25 correction uses the ICES stock annex boundary accepted in the earlier
# workflow. These rectangles are required even if absent from the current data.
fu25_rectangles <- c("15E0", "15E1", "16E0", "16E1", "17E1")
rectangle_ids <- sort(unique(na.omit(c(dat$CLstatisticalRectangle, fu25_rectangles))))

rectangle_cache_files <- list.files(
  cache_dir,
  "^rectangles_[0-9]+[.]rds$",
  full.names = TRUE
)
rectangle_pool <- if (!refresh_reference && length(rectangle_cache_files)) {
  do.call(rbind, lapply(rectangle_cache_files, function(path) readRDS(path)$data))
} else {
  NULL
}
if (!is.null(rectangle_pool)) {
  rectangle_pool <- rectangle_pool[!duplicated(rectangle_pool$ICESNAME), ]
}

rectangle_batches <- split(rectangle_ids, ceiling(seq_along(rectangle_ids) / 100))
message("  Preparing ", length(rectangle_ids), " rectangles in ", length(rectangle_batches), " batches")
rectangles <- do.call(rbind, lapply(seq_along(rectangle_batches), function(i) {
  ids <- rectangle_batches[[i]]
  path <- file.path(cache_dir, paste0("rectangles_", i, ".rds"))
  previous <- if (file.exists(path) && !refresh_reference) readRDS(path) else NULL
  if (!is.null(previous) && identical(previous$ids, ids)) return(previous$data)
  if (!is.null(rectangle_pool) && all(ids %in% rectangle_pool$ICESNAME)) {
    return(rectangle_pool[match(ids, rectangle_pool$ICESNAME), ])
  }
  x <- read_gis(
    "ICES_Statrec_complete/MapServer/0",
    paste0("ICESNAME IN ('", paste(gsub("'", "''", ids), collapse = "','"), "')")
  )
  saveRDS(list(ids = ids, data = x), path)
  x
}))
rm(rectangle_pool)
stopifnot(
  !anyDuplicated(rectangles$ICESNAME),
  all(fu25_rectangles %in% rectangles$ICESNAME)
)

# 6. Nephrops Functional Units ----------------------------------------------
message("[4/8] Preparing Nephrops Functional Unit stock boundaries")
single_fu <- c(5:17, 19, 22, 25, 30:34)
fu_key <- bind_rows(
  tibble(
    fu_no = paste0("FU", single_fu),
    fishstock = paste0("nep.fu.", single_fu)
  ),
  tribble(
    ~fu_no, ~fishstock,
    "4outFU", "nep.27.4outFU",
    "6aoutFU", "nep.27.6aoutFU",
    "7outFU", "nep.27.7outFU",
    "FU3", "nep.fu.3-4",
    "FU4", "nep.fu.3-4",
    "FU2021", "nep.fu.2021",
    "FU23", "nep.fu.2324",
    "FU24", "nep.fu.2324",
    "FU26", "nep.fu.2627",
    "FU27", "nep.fu.2627",
    "FU28", "nep.fu.2829",
    "FU29", "nep.fu.2829"
  )
)

nephrops_fu <- normalise(st_read(fu_file, quiet = TRUE)) |> select(fu_no)
stopifnot(
  !anyDuplicated(nephrops_fu$fu_no),
  !any(st_is_empty(nephrops_fu)),
  sum(nephrops_fu$fu_no == "FU25") == 1L
)

# Supplied FU25 GIS omits the 2021 extension; use the accepted five-rectangle boundary.
st_geometry(nephrops_fu)[nephrops_fu$fu_no == "FU25"] <-
  union_geometry(rectangles[match(fu25_rectangles, rectangles$ICESNAME), ])

nephrops_fu <- left_join(nephrops_fu, fu_key, by = "fu_no", relationship = "many-to-one")
nephrops_unmapped <- filter(nephrops_fu, is.na(fishstock)) |> st_drop_geometry()
write.csv(nephrops_unmapped, file.path(diagnostic_dir, "nephrops_unmapped_FU.csv"), row.names = FALSE)

nephrops <- nephrops_fu |>
  filter(!is.na(fishstock)) |>
  group_by(fishstock) |>
  summarise(do_union = TRUE, .groups = "drop")
nephrops$complete <- vapply(
  nephrops$fishstock,
  function(s) all(fu_key$fu_no[fu_key$fishstock == s] %in% nephrops_fu$fu_no),
  logical(1)
)

# 7. Build stock boundaries and reported geographical supports --------------
message("[5/8] Building stock boundaries and RDBES geographical supports")
used_stocks <- stock_species |>
  filter(aphia %in% dat$aphia) |>
  distinct(fishstock)

# Same-taxon alternatives from the whole fixed SID catalogue remain available,
# even if their own SID region tag is outside the 11 FisheriesXplorer regions.
ordinary <- used_stocks$fishstock[!grepl("^(nep[.]|san[.]sa[.])", used_stocks$fishstock)]
has_sandeel <- any(startsWith(used_stocks$fishstock, "san.sa."))

area_codes <- unique(na.omit(c(
  dat$CLarea,
  stock_areas$area[stock_areas$fishstock %in% ordinary],
  if (has_sandeel) c("27.3", "27.4")
)))
message("  Preparing ", length(area_codes), " ICES area codes")

atomic_parts <- setNames(
  lapply(area_codes, function(a) {
    which(areas$Area_Full == a | startsWith(areas$Area_Full, paste0(a, ".")))
  }),
  area_codes
)

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
san_envelope <- if (has_sandeel) {
  union_geometry(c(area_polygons[["27.3"]], area_polygons[["27.4"]]))
} else {
  empty
}

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
  if (startsWith(s, "san.sa.")) {
    # Broad candidate envelope only; not treated as a definitive stock boundary.
    stock_polygons[[s]] <- san_envelope
    stock_complete[s] <- TRUE
    next
  }
  codes <- stock_areas$area[stock_areas$fishstock == s]
  parts <- area_polygons[codes]
  stock_polygons[[s]] <- union_parts(unlist(atomic_parts[codes], use.names = FALSE))
  stock_complete[s] <- length(parts) > 0 && all(vapply(
    parts,
    function(g) !st_is_empty(g),
    logical(1)
  ))
}

geography <- dat |>
  distinct(CLarea, CLstatisticalRectangle) |>
  mutate(geo_id = row_number())
support <- rep(list(empty), nrow(geography))
geography$geography_status <- "unknown"
area_groups <- split(
  seq_len(nrow(geography)),
  match(geography$CLarea, unique(geography$CLarea))
)

message("  Building ", nrow(geography), " distinct reported geographical supports")
for (ids in area_groups) {
  area_name <- geography$CLarea[ids[1]]
  a <- area_polygons[[area_name]]
  if (is.null(a)) a <- empty
  r <- match(geography$CLstatisticalRectangle[ids], rectangles$ICESNAME)

  if (!st_is_empty(a)) {
    support[ids[is.na(r)]] <- rep(list(a), sum(is.na(r)))
    geography$geography_status[ids] <- ifelse(
      is.na(r), "area_only", "area_rectangle_conflict"
    )
  }

  has_r <- which(!is.na(r))
  if (!length(has_r)) next
  target <- if (!st_is_empty(a)) a else st_geometry(areas)
  cut <- suppressWarnings(st_intersection(st_geometry(rectangles)[r[has_r]], target))
  if (!length(cut)) next

  pieces <- split(seq_along(cut), attr(cut, "idx")[, 1L])
  for (key in names(pieces)) {
    g <- union_geometry(cut[pieces[[key]]])
    if (sum(as.numeric(st_area(g))) <= 0) next
    i <- ids[has_r[as.integer(key)]]
    support[[i]] <- g
    geography$geography_status[i] <- if (!st_is_empty(a)) {
      "area_and_rectangle"
    } else {
      "rectangle_only"
    }
  }
}
support_sfc <- do.call(c, support)

# 8. Match each species/location combination to candidate stocks ------------
message("[6/8] Matching RDBES species/locations to ICES stocks")
# These stock definitions need information beyond the current generic area
# relationships. They remain unresolved unless an explicit production rule exists.
needs_review <- c(
  "ane.27.9aS", "ane.27.9aW",
  "cap.27.1-2", "cap.27.2a514",
  "pra.27.4a", "pra.27.3a4a",
  "her.27.nirs", "her.27.irls", "her.27.6aN", "her.27.6aS7bc",
  "her.27.20-24", "her.27.3a47d", "her.27.1-24a514a", "her.27.5a",
  "cod.27.1-2", "cod.27.1-2.coastN", "cod.27.2.coastS",
  "cod.27.22-24", "cod.27.24-32", "cod.27.21", "cod.27.46a7d20",
  "reb.2127.dp", "reb.2127.sp", "reb.27.5a14", "reb.27.14b",
  "seh.27.1", "seh.27.125a14", "cyo.27.nea", "guq.27.nea"
)

dat <- left_join(
  dat,
  geography,
  by = c("CLarea", "CLstatisticalRectangle"),
  relationship = "many-to-one"
)

rules <- dat |>
  distinct(geo_id, aphia, CLarea) |>
  filter(!is.na(aphia)) |>
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
  s <- rules$fishstock[i]
  a <- rules$CLarea[i]
  valid_a <- !is.na(a) && !is.null(area_polygons[[a]]) && !st_is_empty(area_polygons[[a]])

  if (valid_a && (a == "27" || startsWith(a, "27.")) &&
      s %in% c("cod.21.1.isc", "cod.21.1.osc")) {
    return("outside")
  }
  if (valid_a && startsWith(a, "27.") && s == "cod.21.27.1.14" &&
      !(a == "27.14" || startsWith(a, "27.14."))) {
    return("outside")
  }

  z <- rules$evidence[i]
  if (z == "outside" && !stock_complete[s]) z <- "unknown"
  if (z != "outside" && (s %in% needs_review || startsWith(s, "san.sa."))) z <- "review"
  z
}, character(1))

candidates <- dat |>
  select(row_id, geo_id, aphia) |>
  inner_join(
    filter(rules, evidence != "outside") |> select(-CLarea),
    by = c("geo_id", "aphia"),
    relationship = "many-to-many"
  )

matches <- candidates |>
  summarise(
    n_candidates = n(),
    candidate_stocks = paste(fishstock, collapse = ";"),
    candidate_evidence = paste(paste(fishstock, evidence, sep = ":"), collapse = ";"),
    fishstock = if (n() == 1L && first(evidence) == "full") first(fishstock) else NA_character_,
    .by = row_id
  )

dat <- left_join(dat, matches, by = "row_id", relationship = "one-to-one") |>
  mutate(
    n_candidates = coalesce(n_candidates, 0L),
    match_status = case_when(
      is.na(aphia) ~ "unresolved_taxon",
      !is.na(fishstock) ~ "unique",
      n_candidates == 0L ~ "unmatched",
      TRUE ~ "ambiguous_or_needs_review"
    ),
    special_rule = "none"
  )

# 9. Accepted production stock-assignment rule: sardine 25E4 / 25E5 --------
message("[7/8] Applying accepted production stock-assignment rules")
if (!"pil.27.8abd" %in% stock_list$fishstock) {
  stop("Production sardine target stock pil.27.8abd is absent from the fixed SID catalogue.")
}

sardine_transfer <-
  dat$match_status == "unique" &
  dat$CLspeciesFaoCode == "PIL" &
  dat$fishstock == "pil.27.7" &
  dat$CLstatisticalRectangle %in% c("25E4", "25E5")

if (any(sardine_transfer)) {
  sardine_rule_audit <- dat[sardine_transfer, ] |>
    summarise(
      Landings_kg = sum(CLscientificWeight),
      .by = c(CLyear, CLarea, CLstatisticalRectangle, fishstock)
    ) |>
    mutate(
      Proposed_stock = "pil.27.8abd",
      special_rule = "sardine_25E4_25E5",
      Country_condition_verified = FALSE
    )
  write.csv(
    sardine_rule_audit,
    file.path(diagnostic_dir, "sardine_25E4_25E5_production_rule.csv"),
    row.names = FALSE
  )

  dat$fishstock[sardine_transfer] <- "pil.27.8abd"
  dat$special_rule[sardine_transfer] <- "sardine_25E4_25E5"
} else {
  warning("No PIL records met the sardine 25E4/25E5 production rule in this run.")
}

# 10. Diagnostics and checkpoint --------------------------------------------
message("[8/8] Saving global stock-matching diagnostics and checkpoint")
stopifnot(
  !anyDuplicated(dat$row_id),
  all(dat$fishstock[dat$match_status == "unique"] %in% stock_list$fishstock)
)

matching_summary <- dat |>
  summarise(
    Landings_kg = sum(CLscientificWeight),
    Rows = n(),
    .by = c(CLyear, match_status)
  )

assignment_summary <- dat |>
  filter(match_status == "unique", !is.na(fishstock)) |>
  summarise(
    Assigned_kg = sum(CLscientificWeight),
    .by = c(CLyear, fishstock, special_rule)
  )

unresolved_landings <- dat |>
  filter(match_status != "unique")

write.csv(matching_summary, file.path(diagnostic_dir, "matching_summary.csv"), row.names = FALSE)
write.csv(assignment_summary, file.path(diagnostic_dir, "stock_assignment_summary.csv"), row.names = FALSE)
write.csv(unresolved_landings, file.path(diagnostic_dir, "unresolved_landings.csv"), row.names = FALSE)

save(
  dat,
  candidates,
  stock_list,
  reference_year,
  matching_scope,
  input_kg,
  supported_ecoregions,
  file = checkpoint_file,
  compress = FALSE
)

writeLines(
  c(
    paste("Checkpoint:", normalizePath(checkpoint_file)),
    paste("Fixed SID reference year:", reference_year),
    paste("RDBES years:", paste(analysed_years, collapse = ", ")),
    paste("Supported ecoregions:", paste(supported_ecoregions, collapse = "; ")),
    "Stock matching is global across the supported ecoregion taxon scope.",
    "Same-taxon alternative stocks remain available as competitors during matching.",
    "Sardine PIL records in rectangles 25E4/25E5 assigned to pil.27.7 are reassigned to pil.27.8abd as a production rule.",
    "The sardine country condition cannot be verified because country is not retained in the compact matching table.",
    "No ecoregion share is calculated in this script."
  ),
  file.path(diagnostic_dir, "RUN_NOTES.txt")
)

print(matching_summary)
message("Saved global stock-matching checkpoint: ", checkpoint_file)
