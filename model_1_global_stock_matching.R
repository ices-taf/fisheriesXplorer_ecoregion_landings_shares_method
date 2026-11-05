# Global RDBES landings -> ICES stock matching (rectangle-first experiment)
#
# Purpose:
#   Match RDBES landings to the fixed SID stock catalogue once, across all
#   FisheriesXplorer ecoregions, using statistical rectangles and ICES area
#   codes rather than building polygon intersections for ordinary stocks.
#
# Experimental design relative to the validated baseline:
#   - The RDBES input, SID catalogue, ICES Vocabulary relationships, review
#     list, Nephrops source, and sardine production rule are unchanged.
#   - Ordinary-stock geography uses the official 2017 ICES statistical-
#     rectangle -> ICES-area mapping in boot/initial/data.
#   - Reported CLarea is used to disambiguate rectangles that cross area
#     boundaries and as a fallback when the rectangle is unavailable.
#   - Geometry is retained only for Nephrops Functional Units, which cannot be
#     represented by ordinary ICES area codes alone.
#   - No ecoregion share is calculated here; model_2 is unchanged.

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

if (!exists("reference_year", inherits = FALSE)) {
  reference_year <- 2026L
}

if (!exists("analysis_years", inherits = FALSE)) {
  analysis_years <- NULL
}

if (!exists("rebuild_landings", inherits = FALSE)) {
  rebuild_landings <- FALSE
}

if (!exists("refresh_reference", inherits = FALSE)) {
  refresh_reference <- FALSE
}

input_file <- "./data/extendedRDBES.RData"
compact_file <- "./data/rdbes_landings_compact.rds"
checkpoint_file <- "./model/rdbes_stock_matches.RData"

fu_file <- paste0(
  "./boot/initial/data/ICES_nephrops_FU/",
  "Nephrops Functional Units.shp"
)

statrec_map_dir <- "./boot/initial/data/ICES_StatRec_mapto_ICES_Areas"

cache_dir <- "./data/ices_cache_simple"

diagnostic_dir <- "./model/stock_matching_rectangle_first"

dir.create(
  "./model",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  cache_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  diagnostic_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# 2. Helpers ----------------------------------------------------------------

clean <- function(x) {
  na_if(trimws(as.character(x)), "")
}


cached <- function(name, expression) {

  path <- file.path(
    cache_dir,
    paste0(name, ".rds")
  )

  if (
    file.exists(path) &&
    !refresh_reference
  ) {

    message(
      "  Reading cache: ",
      name
    )

    return(
      readRDS(path)
    )
  }

  message(
    "  Preparing: ",
    name
  )

  started <- proc.time()[[3]]

  x <- force(expression)

  saveRDS(
    x,
    path
  )

  message(
    "  Finished ",
    name,
    " in ",
    round(
      proc.time()[[3]] - started,
      1
    ),
    " seconds"
  )

  x
}


normalise <- function(x) {

  st_make_valid(
    st_set_precision(
      st_transform(
        x,
        3857
      ),
      1000
    )
  )
}


empty <- st_sfc(
  st_geometrycollection(),
  crs = 3857
)


union_geometry <- function(x) {

  g <- st_geometry(x)

  if (
    !length(g) ||
    all(st_is_empty(g))
  ) {
    return(empty)
  }

  if (length(g) == 1L) {
    return(g)
  }

  normalise(
    st_union(g)
  )
}


# Read current ICES GIS layers.
# Used only for Nephrops in the hybrid workflow.
read_gis <- function(layer, where = "1=1") {

  url <- paste0(
    "https://gis.ices.dk/gis/rest/services/ICES_reference_layers/",
    layer,
    "/query"
  )

  response <- httr::GET(
    url,
    query = list(
      where = where,
      returnIdsOnly = "true",
      f = "json"
    ),
    httr::timeout(180)
  )

  httr::stop_for_status(response)

  info <- jsonlite::fromJSON(
    httr::content(
      response,
      "text",
      encoding = "UTF-8"
    )
  )

  if (
    !is.null(info$error) ||
    is.null(info$objectIdFieldName)
  ) {
    stop("GIS ID query failed: ", layer)
  }

  ids <- info$objectIds

  if (!length(ids)) {
    return(
      st_sf(
        geometry = st_sfc(crs = 3857)
      )
    )
  }

  pages <- lapply(
    split(
      ids,
      ceiling(seq_along(ids) / 100)
    ),
    function(wanted) {

      path <- tempfile(
        fileext = ".geojson"
      )

      on.exit(
        unlink(path)
      )

      r <- httr::GET(
        url,
        query = list(
          objectIds = paste(
            wanted,
            collapse = ","
          ),
          outFields = "*",
          outSR = 3857,
          f = "geojson"
        ),
        httr::timeout(180),
        httr::write_disk(path)
      )

      httr::stop_for_status(r)

      x <- st_read(
        path,
        quiet = TRUE
      )

      stopifnot(
        nrow(x) == length(wanted)
      )

      normalise(x)
    }
  )

  do.call(
    rbind,
    pages
  )
}

# Used only for rectangle-vs-Nephrops-FU geometry.
relation_vec <- function(
  support,
  target
) {

  support <- sf::st_geometry(support)
  target <- sf::st_geometry(target)

  answer <- rep(
    "unknown",
    length(support)
  )

  if (
    !length(support) ||
    !length(target) ||
    any(sf::st_is_empty(target))
  ) {
    return(answer)
  }

  stopifnot(
    length(target) == 1L
  )

  total <- as.numeric(
    sf::st_area(support)
  )

  tolerance <- pmax(
    1,
    total * 1e-10
  )

  valid <- which(
    !sf::st_is_empty(support) &
      is.finite(total) &
      total > 0
  )

  answer[valid] <- "outside"

  valid <- valid[
    total[valid] >
      tolerance[valid]
  ]

  if (!length(valid)) {
    return(answer)
  }

  hits <- valid[
    sf::st_intersects(
      target,
      support[valid]
    )[[1L]]
  ]

  if (!length(hits)) {
    return(answer)
  }

  full <- hits[
    sf::st_covers(
      target,
      support[hits]
    )[[1L]]
  ]

  answer[full] <- "full"

  boundary <- setdiff(
    hits,
    full
  )

  if (!length(boundary)) {
    return(answer)
  }

  cut <- suppressWarnings(
    sf::st_intersection(
      support[boundary],
      target
    )
  )

  inside <- numeric(
    length(boundary)
  )

  if (length(cut)) {

    inside[
      attr(cut, "idx")[, 1L]
    ] <- as.numeric(
      sf::st_area(cut)
    )
  }

  answer[boundary] <- ifelse(
    inside <= tolerance[boundary],
    "outside",
    ifelse(
      total[boundary] - inside <=
        tolerance[boundary],
      "full",
      "partial"
    )
  )

  answer
}


split_sid_ecoregions <- function(x) {

  lapply(
    strsplit(
      coalesce(x, ""),
      ",",
      fixed = TRUE
    ),
    trimws
  )
}


# The static statistical-rectangle map contains Area 27 subdivisions
# without the leading "27.", e.g. "4.b", "12.a.1".
normalise_static_area27 <- function(x) {

  x <- tolower(
    clean(x)
  )

  ifelse(
    is.na(x),
    NA_character_,
    ifelse(
      startsWith(x, "27."),
      x,
      paste0("27.", x)
    )
  )
}


# TRUE when 'parent' is the same ICES area as 'child',
# or a broader parent.
area_contains <- function(
  parent,
  child
) {

  if (
    is.na(parent) ||
    is.na(child)
  ) {
    return(FALSE)
  }

  identical(
    parent,
    child
  ) ||
    startsWith(
      child,
      paste0(
        parent,
        "."
      )
    )
}


# Given one area from the static rectangle map and one reported CLarea,
# return their most specific compatible representation.
# NA means conflict.
most_specific_overlap <- function(
  mapped_area,
  reported_area
) {

  if (
    is.na(mapped_area) &&
    is.na(reported_area)
  ) {
    return(NA_character_)
  }

  if (is.na(mapped_area)) {
    return(reported_area)
  }

  if (is.na(reported_area)) {
    return(mapped_area)
  }

  if (
    area_contains(
      mapped_area,
      reported_area
    )
  ) {
    return(reported_area)
  }

  if (
    area_contains(
      reported_area,
      mapped_area
    )
  ) {
    return(mapped_area)
  }

  NA_character_
}


# Classify a set of possible reported areas against a stock's
# ICES-area set.
#
# full:
#   every possible reported area is wholly covered by the stock areas.
#
# outside:
#   no possible reported area overlaps the stock areas hierarchically.
#
# partial:
#   some overlap exists, but the available area information is not
#   sufficient to prove that the whole reported location is in stock.
classify_area_set <- function(
  location_areas,
  stock_codes
) {

  location_areas <- unique(
    na.omit(
      location_areas
    )
  )

  stock_codes <- unique(
    na.omit(
      stock_codes
    )
  )

  if (
    !length(location_areas) ||
    !length(stock_codes)
  ) {
    return("unknown")
  }

  status <- vapply(
    location_areas,
    function(loc) {

      if (
        any(
          vapply(
            stock_codes,
            area_contains,
            logical(1),
            child = loc
          )
        )
      ) {
        return("full")
      }

      if (
        any(
          vapply(
            stock_codes,
            function(s) {
              area_contains(
                loc,
                s
              )
            },
            logical(1)
          )
        )
      ) {
        return("partial")
      }

      "outside"
    },
    character(1)
  )

  if (
    all(
      status == "full"
    )
  ) {
    return("full")
  }

  if (
    all(
      status == "outside"
    )
  ) {
    return("outside")
  }

  "partial"
}


# 3. Aggregate landings once -------------------------------------------------

prepare_landings <- function(
  input_file
) {

  if (
    !rebuild_landings &&
    exists(
      "cl.all2",
      envir = .GlobalEnv,
      inherits = FALSE
    )
  ) {

    message(
      "  Reusing cl.all2 already in memory"
    )

    source <- get(
      "cl.all2",
      envir = .GlobalEnv
    )

  } else {

    message(
      "  Loading the original RData once"
    )

    loaded <- new.env(
      parent = emptyenv()
    )

    load(
      input_file,
      envir = loaded
    )

    stopifnot(
      exists(
        "cl.all2",
        envir = loaded,
        inherits = FALSE
      )
    )

    source <- loaded$cl.all2
  }

  required <- c(
    "CLyear",
    "CLarea",
    "CLstatisticalRectangle",
    "CLscientificWeight"
  )

  stopifnot(
    all(
      required %in%
        names(source)
    ),
    is.numeric(
      source$CLscientificWeight
    ),
    any(
      c(
        "CLspeciesCode",
        "CLspeciesFaoCode"
      ) %in%
        names(source)
    )
  )

  source <- select(
    source,
    all_of(required),
    any_of(
      c(
        "CLspeciesCode",
        "CLspeciesFaoCode",
        "CLcatchCategory"
      )
    )
  )

  if (
    exists(
      "loaded",
      inherits = FALSE
    )
  ) {
    rm(loaded)
  }

  if (
    "CLcatchCategory" %in%
      names(source)
  ) {

    source <- filter(
      source,
      tolower(
        trimws(
          CLcatchCategory
        )
      ) == "lan"
    )

  } else {

    message(
      "  CLcatchCategory absent: input must already contain landings only"
    )
  }

  source <- filter(
    source,
    !is.na(
      CLscientificWeight
    )
  )

  stopifnot(
    nrow(source) > 0,
    all(
      is.finite(
        source$CLscientificWeight
      )
    ),
    all(
      source$CLscientificWeight >= 0
    )
  )

  for (
    nm in c(
      "CLspeciesCode",
      "CLspeciesFaoCode"
    )
  ) {

    if (
      !nm %in%
        names(source)
    ) {

      source[[nm]] <- NA_character_
    }
  }

  message(
    "  Aggregating ",
    format(
      nrow(source),
      big.mark = ","
    ),
    " landing rows; keeping ALL areas"
  )

  result <- source |>
    mutate(
      CLyear = as.integer(
        as.character(
          CLyear
        )
      ),
      CLarea = tolower(
        clean(
          CLarea
        )
      ),
      CLstatisticalRectangle = toupper(
        clean(
          CLstatisticalRectangle
        )
      ),
      supplied_aphia = clean(
        CLspeciesCode
      ),
      CLspeciesFaoCode = toupper(
        clean(
          CLspeciesFaoCode
        )
      )
    ) |>
    summarise(
      CLscientificWeight = sum(
        CLscientificWeight
      ),
      .by = c(
        CLyear,
        CLarea,
        CLstatisticalRectangle,
        supplied_aphia,
        CLspeciesFaoCode
      )
    )

  stopifnot(
    !anyNA(
      result$CLyear
    ),
    isTRUE(
      all.equal(
        sum(
          result$CLscientificWeight
        ),
        sum(
          source$CLscientificWeight
        )
      )
    )
  )

  result
}


message(
  "[1/8] Reading or preparing compact RDBES landings"
)

compact_is_current <-
  file.exists(
    compact_file
  ) &&
  file.exists(
    input_file
  ) &&
  file.info(
    compact_file
  )$mtime >=
  file.info(
    input_file
  )$mtime

if (
  compact_is_current &&
  !rebuild_landings
) {

  dat <- readRDS(
    compact_file
  )

} else {

  if (
    file.exists(
      compact_file
    ) &&
    !rebuild_landings
  ) {

    message(
      "  extendedRDBES.RData is newer than the compact file; rebuilding compact landings"
    )
  }

  dat <- prepare_landings(
    input_file
  )

  saveRDS(
    dat,
    compact_file,
    compress = FALSE
  )
}


dat <- select(
  dat,
  CLyear,
  CLarea,
  CLstatisticalRectangle,
  supplied_aphia,
  CLspeciesFaoCode,
  CLscientificWeight
)

if (
  !is.null(
    analysis_years
  )
) {

  dat <- filter(
    dat,
    CLyear %in%
      analysis_years
  )
}

stopifnot(
  nrow(dat) > 0,
  !anyNA(
    dat$CLyear
  ),
  !anyNA(
    dat$CLscientificWeight
  ),
  all(
    is.finite(
      dat$CLscientificWeight
    )
  ),
  all(
    dat$CLscientificWeight >= 0
  )
)

dat$row_id <- seq_len(
  nrow(dat)
)

input_kg_before_scope <- sum(
  dat$CLscientificWeight
)

analysed_years <- sort(
  unique(
    dat$CLyear
  )
)


# 4. Fixed SID catalogue and ICES Vocabulary relationships ------------------

message(
  "[2/8] Reading fixed SID catalogue and ICES Vocabulary relationships"
)

stock_list <- cached(
  paste0(
    "stocks_",
    reference_year
  ),
  icesSD::getSD(
    year = reference_year
  )
) |>
  transmute(
    fishstock = as.character(
      StockKeyLabel
    ),
    description = as.character(
      StockKeyDescription
    ),
    SID_ecoregions = clean(
      EcoRegion
    )
  ) |>
  distinct(
    fishstock,
    .keep_all = TRUE
  )

stopifnot(
  nrow(stock_list) > 0L,
  !anyDuplicated(
    stock_list$fishstock
  )
)


stock_species <- cached(
  "stock_species",
  icesVocab::getCodeTypeRelation(
    "ICES_StockCode",
    "SpecWoRMS"
  )
) |>
  transmute(
    fishstock = as.character(
      ICES_StockCode
    ),
    aphia = as.character(
      SpecWoRMS
    )
  ) |>
  semi_join(
    stock_list,
    by = "fishstock"
  ) |>
  distinct()


stock_areas <- cached(
  "stock_areas",
  icesVocab::getCodeTypeRelation(
    "ICES_StockCode",
    "ICES_Area"
  )
) |>
  transmute(
    fishstock = as.character(
      ICES_StockCode
    ),
    # Keep ICES Vocabulary area codes exactly as supplied.
    # Examples include broad areas such as "27" and "37",
    # as well as subdivisions such as "27.8.a".
    area = tolower(
      clean(
        ICES_Area
      )
    )
  ) |>
  filter(
    !is.na(area)
  ) |>
  semi_join(
    stock_list,
    by = "fishstock"
  ) |>
  distinct()


fao_lookup <- cached(
  "asfis_worms",
  icesVocab::getCodeTypeRelation(
    "SpecASFIS",
    "SpecWoRMS"
  )
) |>
  transmute(
    CLspeciesFaoCode = as.character(
      SpecASFIS
    ),
    fallback = as.character(
      SpecWoRMS
    )
  ) |>
  distinct() |>
  summarise(
    fallback = if (
      n_distinct(
        fallback
      ) == 1L
    ) {
      first(
        fallback
      )
    } else {
      NA_character_
    },
    .by = CLspeciesFaoCode
  )


dat <- dat |>
  left_join(
    fao_lookup,
    by = "CLspeciesFaoCode",
    relationship = "many-to-one"
  ) |>
  mutate(
    taxon_conflict =
      !is.na(
        supplied_aphia
      ) &
      !is.na(
        fallback
      ) &
      supplied_aphia !=
      fallback,

    aphia = if_else(
      taxon_conflict,
      NA_character_,
      coalesce(
        supplied_aphia,
        fallback
      )
    ),

    aphia = if_else(
      !is.na(
        aphia
      ) &
      grepl(
        "^[0-9]+$",
        aphia
      ),
      aphia,
      NA_character_
    )
  )


# Scope matching to taxa belonging to stocks in the 11 supported
# ecoregions, while keeping stocks with missing SID region tags
# and same-taxon competitors.

sid_tags <- split_sid_ecoregions(
  stock_list$SID_ecoregions
)

documented_scope_additions <- c(
  "mon.27.78abd",
  "rjb.27.89a"
)


stock_list <- stock_list |>
  mutate(
    In_supported_ecoregions = vapply(
      sid_tags,
      function(z) {
        any(
          z %in%
            supported_ecoregions
        )
      },
      logical(1)
    ),

    SID_region_missing =
      is.na(
        SID_ecoregions
      ) |
      SID_ecoregions == "",

    In_matching_scope =
      In_supported_ecoregions |
      SID_region_missing |
      fishstock %in%
      documented_scope_additions
  )


scope_stocks <- filter(
  stock_list,
  In_matching_scope
)


selected_taxa <- unique(
  stock_species$aphia[
    stock_species$fishstock %in%
      scope_stocks$fishstock
  ]
)


dat <- dat |>
  mutate(
    scope_status = case_when(
      is.na(
        aphia
      ) ~
        "unresolved_taxon_retained",

      aphia %in%
        selected_taxa ~
        "selected_taxon_retained",

      TRUE ~
        "other_taxon_excluded"
    )
  )


scope_summary <- dat |>
  summarise(
    Rows = n(),
    Landings_kg = sum(
      CLscientificWeight
    ),
    .by = c(
      CLyear,
      scope_status,
      CLspeciesFaoCode,
      aphia
    )
  )


stopifnot(
  isTRUE(
    all.equal(
      sum(
        scope_summary$Landings_kg
      ),
      input_kg_before_scope
    )
  )
)


write.csv(
  scope_summary,
  file.path(
    diagnostic_dir,
    "landings_scope.csv"
  ),
  row.names = FALSE
)


write.csv(
  scope_stocks,
  file.path(
    diagnostic_dir,
    "stock_scope_catalogue.csv"
  ),
  row.names = FALSE
)


dat <- filter(
  dat,
  scope_status !=
    "other_taxon_excluded"
) |>
  select(
    -scope_status
  )


if (!nrow(dat)) {
  stop(
    "No RDBES landings remain in the supported ecoregion taxon scope."
  )
}


input_kg <- sum(
  dat$CLscientificWeight
)


matching_scope <- list(
  ecoregions = supported_ecoregions,
  years = analysed_years,
  included_aphia = selected_taxa,
  geography_method = "rectangle_first_static_2017"
)


message(
  "  Fixed SID ",
  reference_year,
  ": ",
  nrow(stock_list),
  " active stocks; ",
  nrow(scope_stocks),
  " stocks retained in matching scope"
)


message(
  "  Retained ",
  nrow(dat),
  " aggregate RDBES rows across ALL reported areas"
)


# 5. Static statistical-rectangle -> ICES-area lookup -----------------------

message(
  "[3/8] Reading official static statistical-rectangle -> ICES-area mapping"
)


if (
  !dir.exists(
    statrec_map_dir
  )
) {

  stop(
    "Static rectangle mapping folder not found: ",
    statrec_map_dir
  )
}


statrec_shp <- list.files(
  statrec_map_dir,
  pattern = "[.]shp$",
  full.names = TRUE,
  ignore.case = TRUE
)


if (
  length(
    statrec_shp
  ) != 1L
) {

  stop(
    "Expected exactly one .shp in ",
    statrec_map_dir,
    "; found ",
    length(
      statrec_shp
    )
  )
}


statrec_map <- st_read(
  statrec_shp,
  quiet = TRUE
)


required_map_cols <- c(
  "ICESNAME",
  "AreasList"
)


if (
  !all(
    required_map_cols %in%
      names(statrec_map)
  )
) {

  stop(
    "Static rectangle map must contain: ",
    paste(
      required_map_cols,
      collapse = ", "
    )
  )
}


stopifnot(
  !anyDuplicated(
    statrec_map$ICESNAME
  )
)


rect_area_lookup <- statrec_map |>
  st_drop_geometry() |>
  transmute(
    CLstatisticalRectangle =
      toupper(
        clean(
          ICESNAME
        )
      ),
    area = as.character(
      AreasList
    )
  ) |>
  filter(
    !is.na(
      CLstatisticalRectangle
    ),
    !is.na(
      area
    ),
    nzchar(
      trimws(
        area
      )
    )
  ) |>
  separate_rows(
    area,
    sep = ","
  ) |>
  mutate(
    area = normalise_static_area27(
      area
    )
  ) |>
  filter(
    !is.na(
      area
    )
  ) |>
  distinct(
    CLstatisticalRectangle,
    area
  )


rect_area_lists <- split(
  rect_area_lookup$area,
  rect_area_lookup$CLstatisticalRectangle
)


boundary_rectangles <- rect_area_lookup |>
  summarise(
    n_areas = n_distinct(
      area
    ),
    .by = CLstatisticalRectangle
  ) |>
  filter(
    n_areas > 1L
  )


write.csv(
  boundary_rectangles,
  file.path(
    diagnostic_dir,
    "static_boundary_rectangles.csv"
  ),
  row.names = FALSE
)


# Resolve each distinct reported RDBES location to the most specific
# compatible ICES area(s).
#
# Most rectangles have one mapped area.
# For boundary rectangles, CLarea narrows the possible areas.
#
# If rectangle and CLarea conflict, no area is assumed and the stock
# evidence remains unknown.

geography <- dat |>
  distinct(
    CLarea,
    CLstatisticalRectangle
  ) |>
  mutate(
    geo_id = row_number(),
    location_areas = vector(
      "list",
      n()
    ),
    geography_status = NA_character_
  )


for (
  i in seq_len(
    nrow(
      geography
    )
  )
) {

  rec <- geography$CLstatisticalRectangle[i]
  reported <- geography$CLarea[i]

  mapped <- if (
    !is.na(rec) &&
    rec %in%
      names(
        rect_area_lists
      )
  ) {

    unique(
      rect_area_lists[[rec]]
    )

  } else {

    character()
  }


  if (length(mapped)) {

    if (
      !is.na(
        reported
      )
    ) {

      narrowed <- unique(
        na.omit(
          vapply(
            mapped,
            most_specific_overlap,
            character(1),
            reported_area = reported
          )
        )
      )


      if (
        length(
          narrowed
        )
      ) {

        geography$location_areas[[i]] <- narrowed

        geography$geography_status[i] <-
          if (
            length(
              mapped
            ) > 1L
          ) {

            "boundary_rectangle_resolved_by_CLarea"

          } else {

            "rectangle_and_CLarea"
          }

      } else {

        geography$location_areas[[i]] <- character()

        geography$geography_status[i] <-
          "rectangle_CLarea_conflict"
      }

    } else {

      geography$location_areas[[i]] <- mapped

      geography$geography_status[i] <-
        if (
          length(
            mapped
          ) > 1L
        ) {

          "boundary_rectangle_without_CLarea"

        } else {

          "rectangle_only"
        }
    }

  } else if (
    !is.na(
      reported
    )
  ) {

    geography$location_areas[[i]] <- reported

    geography$geography_status[i] <-
      if (
        is.na(rec)
      ) {

        "CLarea_only"

      } else {

        "rectangle_not_in_static_map_CLarea_fallback"
      }

  } else {

    geography$location_areas[[i]] <- character()

    geography$geography_status[i] <-
      "no_usable_geography"
  }
}


geography_audit <- geography |>
  mutate(
    location_area_count = lengths(
      location_areas
    )
  ) |>
  summarise(
    Locations = n(),
    .by = c(
      geography_status,
      location_area_count
    )
  ) |>
  arrange(
    geography_status,
    location_area_count
  )


write.csv(
  geography_audit,
  file.path(
    diagnostic_dir,
    "location_resolution_summary.csv"
  ),
  row.names = FALSE
)


used_rectangles <- unique(
  na.omit(
    dat$CLstatisticalRectangle
  )
)


missing_static_rectangles <- setdiff(
  used_rectangles,
  names(
    rect_area_lists
  )
)


write.csv(
  tibble(
    CLstatisticalRectangle =
      sort(
        missing_static_rectangles
      )
  ),
  file.path(
    diagnostic_dir,
    "RDBES_rectangles_missing_from_static_map.csv"
  ),
  row.names = FALSE
)


# 6. Nephrops Functional Units ---------------------------------------------
#
# Ordinary stocks use the simplified rectangle/ICES-area method.
#
# Nephrops deliberately retains the validated baseline spatial treatment:
#
# reported ICES area ∩ statistical rectangle
#                     ↓
#             geographic support
#                     ↓
#          Nephrops FU polygon
#
# This is required because Nephrops Functional Units do not correspond
# directly to ordinary ICES area codes.

message(
  "[4/8] Preparing Nephrops Functional Units using baseline spatial method"
)


# -------------------------------------------------------------------------
# 6a. Identify Nephrops stocks and RDBES locations
# -------------------------------------------------------------------------

single_fu <- c(
  5:17,
  19,
  22,
  25,
  30:34
)


fu_key <- bind_rows(

  tibble(
    fu_no = paste0(
      "FU",
      single_fu
    ),
    fishstock = paste0(
      "nep.fu.",
      single_fu
    )
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


nep_aphia <- unique(
  stock_species$aphia[
    startsWith(
      stock_species$fishstock,
      "nep."
    )
  ]
)


nep_geography <- dat |>
  filter(
    aphia %in% nep_aphia
  ) |>
  distinct(
    CLarea,
    CLstatisticalRectangle
  ) |>
  mutate(
    nep_geo_id = row_number()
  )


message(
  "  Nephrops: ",
  nrow(nep_geography),
  " distinct reported area/rectangle locations"
)


# -------------------------------------------------------------------------
# 6b. Load current ICES area polygons
# -------------------------------------------------------------------------

areas <- cached(
  "areas",
  read_gis(
    "ICES_Areas/MapServer/0"
  )
)


nep_area_codes <- unique(
  na.omit(
    nep_geography$CLarea
  )
)


nep_atomic_parts <- setNames(
  lapply(
    nep_area_codes,
    function(a) {

      which(
        areas$Area_Full == a |
          startsWith(
            areas$Area_Full,
            paste0(a, ".")
          )
      )
    }
  ),
  nep_area_codes
)


nep_area_polygons <- lapply(
  nep_atomic_parts,
  function(i) {

    if (!length(i)) {
      return(empty)
    }

    union_geometry(
      areas[i, ]
    )
  }
)


# -------------------------------------------------------------------------
# 6c. Load only statistical rectangles needed by Nephrops
# -------------------------------------------------------------------------

fu25_rectangles <- c(
  "15E0",
  "15E1",
  "16E0",
  "16E1",
  "17E1"
)


nep_rectangle_ids <- sort(
  unique(
    na.omit(
      c(
        nep_geography$CLstatisticalRectangle,
        fu25_rectangles
      )
    )
  )
)


# Reuse rectangles already downloaded by the baseline workflow whenever
# possible.

rectangle_cache_files <- list.files(
  cache_dir,
  "^rectangles_[0-9]+[.]rds$",
  full.names = TRUE
)


rectangle_pool <- if (
  !refresh_reference &&
  length(rectangle_cache_files)
) {

  do.call(
    rbind,
    lapply(
      rectangle_cache_files,
      function(path) {
        readRDS(path)$data
      }
    )
  )

} else {

  NULL
}


if (!is.null(rectangle_pool)) {

  rectangle_pool <-
    rectangle_pool[
      !duplicated(
        rectangle_pool$ICESNAME
      ),
    ]
}


if (
  !is.null(rectangle_pool) &&
  all(
    nep_rectangle_ids %in%
      rectangle_pool$ICESNAME
  )
) {

  message(
    "  Reusing cached ICES rectangle polygons for Nephrops"
  )

  nep_rectangles <- rectangle_pool[
    match(
      nep_rectangle_ids,
      rectangle_pool$ICESNAME
    ),
  ]

} else {

  message(
    "  Downloading ",
    length(nep_rectangle_ids),
    " ICES rectangle polygons needed for Nephrops"
  )

  rectangle_batches <- split(
    nep_rectangle_ids,
    ceiling(
      seq_along(nep_rectangle_ids) / 100
    )
  )


  nep_rectangles <- do.call(
    rbind,
    lapply(
      rectangle_batches,
      function(ids) {

        read_gis(
          "ICES_Statrec_complete/MapServer/0",
          paste0(
            "ICESNAME IN ('",
            paste(
              gsub(
                "'",
                "''",
                ids
              ),
              collapse = "','"
            ),
            "')"
          )
        )
      }
    )
  )
}


rm(rectangle_pool)


stopifnot(
  !anyDuplicated(
    nep_rectangles$ICESNAME
  ),
  all(
    fu25_rectangles %in%
      nep_rectangles$ICESNAME
  )
)


# -------------------------------------------------------------------------
# 6d. Prepare Nephrops FU polygons
# -------------------------------------------------------------------------

nephrops_fu <- normalise(
  st_read(
    fu_file,
    quiet = TRUE
  )
) |>
  select(
    fu_no
  )


stopifnot(
  !anyDuplicated(
    nephrops_fu$fu_no
  ),
  !any(
    st_is_empty(
      nephrops_fu
    )
  ),
  sum(
    nephrops_fu$fu_no ==
      "FU25"
  ) == 1L
)


# Same accepted FU25 correction used in the validated baseline.

st_geometry(
  nephrops_fu
)[
  nephrops_fu$fu_no ==
    "FU25"
] <- union_geometry(
  nep_rectangles[
    match(
      fu25_rectangles,
      nep_rectangles$ICESNAME
    ),
  ]
)


nephrops_fu <- left_join(
  nephrops_fu,
  fu_key,
  by = "fu_no",
  relationship = "many-to-one"
)


nephrops_unmapped <- filter(
  nephrops_fu,
  is.na(fishstock)
) |>
  st_drop_geometry()


write.csv(
  nephrops_unmapped,
  file.path(
    diagnostic_dir,
    "nephrops_unmapped_FU.csv"
  ),
  row.names = FALSE
)


nephrops <- nephrops_fu |>
  filter(
    !is.na(
      fishstock
    )
  ) |>
  group_by(
    fishstock
  ) |>
  summarise(
    do_union = TRUE,
    .groups = "drop"
  )


nephrops$complete <- vapply(
  nephrops$fishstock,
  function(s) {

    all(
      fu_key$fu_no[
        fu_key$fishstock ==
          s
      ] %in%
        nephrops_fu$fu_no
    )
  },
  logical(1)
)


# -------------------------------------------------------------------------
# 6e. Build the same area ∩ rectangle support used by baseline
# -------------------------------------------------------------------------

nep_support <- rep(
  list(empty),
  nrow(nep_geography)
)


nep_geography$geography_status <-
  "unknown"


nep_area_groups <- split(
  seq_len(
    nrow(nep_geography)
  ),
  match(
    nep_geography$CLarea,
    unique(
      nep_geography$CLarea
    )
  )
)


for (ids in nep_area_groups) {

  area_name <-
    nep_geography$CLarea[
      ids[1]
    ]


  a <- nep_area_polygons[[area_name]]


  if (
    is.null(a)
  ) {
    a <- empty
  }


  r <- match(
    nep_geography$CLstatisticalRectangle[
      ids
    ],
    nep_rectangles$ICESNAME
  )


  # Area exists but rectangle is missing:
  # baseline uses the reported ICES area alone.

  if (
    !st_is_empty(a)
  ) {

    missing_rectangle <- ids[
      is.na(r)
    ]

    if (
      length(
        missing_rectangle
      )
    ) {

      nep_support[
        missing_rectangle
      ] <- rep(
        list(a),
        length(
          missing_rectangle
        )
      )
    }

    nep_geography$geography_status[
      ids
    ] <- ifelse(
      is.na(r),
      "area_only",
      "area_rectangle_conflict"
    )
  }


  has_r <- which(
    !is.na(r)
  )


  if (!length(has_r)) {
    next
  }


  # If an ICES area is available, intersect the rectangle with it.
  # Otherwise use the rectangle against all ICES-area polygons.

  target <- if (
    !st_is_empty(a)
  ) {

    a

  } else {

    st_geometry(
      areas
    )
  }


  cut <- suppressWarnings(
    st_intersection(
      st_geometry(
        nep_rectangles
      )[
        r[has_r]
      ],
      target
    )
  )


  if (!length(cut)) {
    next
  }


  pieces <- split(
    seq_along(cut),
    attr(
      cut,
      "idx"
    )[, 1L]
  )


  for (
    key in names(pieces)
  ) {

    g <- union_geometry(
      cut[
        pieces[[key]]
      ]
    )


    if (
      sum(
        as.numeric(
          st_area(g)
        )
      ) <= 0
    ) {
      next
    }


    i <- ids[
      has_r[
        as.integer(key)
      ]
    ]


    nep_support[[i]] <- g


    nep_geography$geography_status[i] <-
      if (
        !st_is_empty(a)
      ) {

        "area_and_rectangle"

      } else {

        "rectangle_only"
      }
  }
}


nep_support_sfc <- do.call(
  c,
  nep_support
)


stopifnot(
  length(
    nep_support_sfc
  ) ==
    nrow(
      nep_geography
    )
)


# -------------------------------------------------------------------------
# 6f. Compare each precise support with each Nephrops FU
# -------------------------------------------------------------------------

nep_evidence <- bind_rows(
  lapply(
    nephrops$fishstock,
    function(s) {

      k <- match(
        s,
        nephrops$fishstock
      )


      z <- relation_vec(
        nep_support_sfc,
        st_geometry(
          nephrops
        )[k]
      )


      # Preserve conservative baseline treatment for incomplete FU maps.

      if (
        !isTRUE(
          nephrops$complete[k]
        )
      ) {

        z[
          z == "outside"
        ] <- "unknown"
      }


      tibble(
        CLarea =
          nep_geography$CLarea,

        CLstatisticalRectangle =
          nep_geography$CLstatisticalRectangle,

        fishstock = s,

        evidence = z
      )
    }
  )
)


write.csv(
  filter(
    nep_evidence,
    evidence !=
      "outside"
  ),
  file.path(
    diagnostic_dir,
    "nephrops_location_candidate_evidence.csv"
  ),
  row.names = FALSE
)


write.csv(
  nep_geography |>
    st_drop_geometry(),
  file.path(
    diagnostic_dir,
    "nephrops_location_resolution.csv"
  ),
  row.names = FALSE
)


# 7. Match species + resolved rectangle/area location to candidate stocks ----

message(
  "[5/8] Matching RDBES species/locations to ICES stocks using area codes"
)


needs_review <- c(

  "ane.27.9aS",
  "ane.27.9aW",

  "cap.27.1-2",
  "cap.27.2a514",

  "pra.27.4a",
  "pra.27.3a4a",

  "her.27.nirs",
  "her.27.irls",
  "her.27.6aN",
  "her.27.6aS7bc",
  "her.27.20-24",
  "her.27.3a47d",
  "her.27.1-24a514a",
  "her.27.5a",

  "cod.27.1-2",
  "cod.27.1-2.coastN",
  "cod.27.2.coastS",
  "cod.27.22-24",
  "cod.27.24-32",
  "cod.27.21",
  "cod.27.46a7d20",

  "reb.2127.dp",
  "reb.2127.sp",
  "reb.27.5a14",
  "reb.27.14b",

  "seh.27.1",
  "seh.27.125a14",

  "cyo.27.nea",
  "guq.27.nea"
)


# Preserve geography list-column only long enough to calculate evidence.

dat <- left_join(
  dat,
  select(
    geography,
    CLarea,
    CLstatisticalRectangle,
    geo_id,
    geography_status
  ),
  by = c(
    "CLarea",
    "CLstatisticalRectangle"
  ),
  relationship = "many-to-one"
)


stock_area_lists <- split(
  stock_areas$area,
  stock_areas$fishstock
)


location_area_lists <- setNames(
  geography$location_areas,
  geography$geo_id
)


rectangle_by_geo <- setNames(
  geography$CLstatisticalRectangle,
  geography$geo_id
)


reported_area_by_geo <- setNames(
  geography$CLarea,
  geography$geo_id
)


rules <- dat |>
  distinct(
    geo_id,
    aphia
  ) |>
  filter(
    !is.na(
      aphia
    )
  ) |>
  inner_join(
    stock_species,
    by = "aphia",
    relationship = "many-to-many"
  )


rules$evidence <- vapply(
  seq_len(
    nrow(
      rules
    )
  ),
  function(i) {

    s <- rules$fishstock[i]

    gid <- as.character(
      rules$geo_id[i]
    )

    loc_areas <-
      location_area_lists[[gid]]

    reported_area <-
      reported_area_by_geo[[gid]]

    rec <-
      rectangle_by_geo[[gid]]


    # Nephrops is the only stock family still using geometry.
    # Evidence is based on the whole statistical rectangle relative
    # to the FU polygon.

    if (
  startsWith(
    s,
    "nep."
  )
) {

  same_area <-
    (
      is.na(nep_evidence$CLarea) &
        is.na(reported_area)
    ) |
    (
      !is.na(nep_evidence$CLarea) &
        !is.na(reported_area) &
        nep_evidence$CLarea ==
        reported_area
    )


  same_rectangle <-
    (
      is.na(
        nep_evidence$CLstatisticalRectangle
      ) &
        is.na(rec)
    ) |
    (
      !is.na(
        nep_evidence$CLstatisticalRectangle
      ) &
        !is.na(rec) &
        nep_evidence$CLstatisticalRectangle ==
        rec
    )


  k <- which(
    nep_evidence$fishstock ==
      s &
      same_area &
      same_rectangle
  )


  if (!length(k)) {
    return("unknown")
  }


  return(
    nep_evidence$evidence[k[1]]
  )
}


    # Sandeel assessment areas are not ordinary ICES_Area stock
    # boundaries in this workflow.
    #
    # Preserve the baseline behaviour: keep locations in 27.3/27.4
    # as review candidates rather than claiming a unique assignment.

    if (
      startsWith(
        s,
        "san.sa."
      )
    ) {

      z <- classify_area_set(
        loc_areas,
        c(
          "27.3",
          "27.4"
        )
      )

      return(
        if (
          z == "outside"
        ) {
          "outside"
        } else {
          "review"
        }
      )
    }


    # Preserve the baseline Iceland/NEA cod exclusions based on
    # reported CLarea.

    if (
      !is.na(
        reported_area
      ) &&
      startsWith(
        reported_area,
        "27."
      ) &&
      s %in%
      c(
        "cod.21.1.isc",
        "cod.21.1.osc"
      )
    ) {

      return(
        "outside"
      )
    }


    if (
      !is.na(
        reported_area
      ) &&
      startsWith(
        reported_area,
        "27."
      ) &&
      s ==
      "cod.21.27.1.14" &&
      !(
        reported_area ==
          "27.14" ||
        startsWith(
          reported_area,
          "27.14."
        )
      )
    ) {

      return(
        "outside"
      )
    }


    z <- classify_area_set(
      loc_areas,
      stock_area_lists[[s]]
    )


    if (
      z != "outside" &&
      s %in%
      needs_review
    ) {

      z <- "review"
    }


    z
  },
  character(1)
)


candidates <- dat |>
  select(
    row_id,
    geo_id,
    aphia
  ) |>
  inner_join(
    filter(
      rules,
      evidence !=
        "outside"
    ),
    by = c(
      "geo_id",
      "aphia"
    ),
    relationship = "many-to-many"
  ) |>
  select(
    row_id,
    geo_id,
    aphia,
    fishstock,
    evidence
  )


matches <- candidates |>
  summarise(
    n_candidates = n(),

    candidate_stocks = paste(
      fishstock,
      collapse = ";"
    ),

    candidate_evidence = paste(
      paste(
        fishstock,
        evidence,
        sep = ":"
      ),
      collapse = ";"
    ),

    fishstock = if (
      n() == 1L &&
      first(
        evidence
      ) == "full"
    ) {
      first(
        fishstock
      )
    } else {
      NA_character_
    },

    .by = row_id
  )


dat <- left_join(
  dat,
  matches,
  by = "row_id",
  relationship = "one-to-one"
) |>
  mutate(
    n_candidates = coalesce(
      n_candidates,
      0L
    ),

    match_status = case_when(

      is.na(
        aphia
      ) ~
        "unresolved_taxon",

      !is.na(
        fishstock
      ) ~
        "unique",

      n_candidates ==
        0L ~
        "unmatched",

      TRUE ~
        "ambiguous_or_needs_review"
    ),

    special_rule = "none"
  )


# 8. Accepted production stock-assignment rule:
#    sardine 25E4 / 25E5 ----------------------------------------------------

message(
  "[6/8] Applying accepted production stock-assignment rules"
)


if (
  !"pil.27.8abd" %in%
    stock_list$fishstock
) {

  stop(
    "Production sardine target stock pil.27.8abd is absent from the fixed SID catalogue."
  )
}


sardine_transfer <-

  dat$match_status ==
    "unique" &

  dat$CLspeciesFaoCode ==
    "PIL" &

  dat$fishstock ==
    "pil.27.7" &

  dat$CLstatisticalRectangle %in%
    c(
      "25E4",
      "25E5"
    )


if (
  any(
    sardine_transfer
  )
) {

  sardine_rule_audit <-
    dat[
      sardine_transfer,
    ] |>
    summarise(
      Landings_kg = sum(
        CLscientificWeight
      ),
      .by = c(
        CLyear,
        CLarea,
        CLstatisticalRectangle,
        fishstock
      )
    ) |>
    mutate(
      Proposed_stock =
        "pil.27.8abd",

      special_rule =
        "sardine_25E4_25E5",

      Country_condition_verified =
        FALSE
    )


  write.csv(
    sardine_rule_audit,
    file.path(
      diagnostic_dir,
      "sardine_25E4_25E5_production_rule.csv"
    ),
    row.names = FALSE
  )


  dat$fishstock[
    sardine_transfer
  ] <- "pil.27.8abd"


  dat$special_rule[
    sardine_transfer
  ] <- "sardine_25E4_25E5"

} else {

  warning(
    "No PIL records met the sardine 25E4/25E5 production rule in this run."
  )
}


message(
  "[7/8] Preparing matching diagnostics"
)


stopifnot(
  !anyDuplicated(
    dat$row_id
  ),

  all(
    dat$fishstock[
      dat$match_status ==
        "unique"
    ] %in%
      stock_list$fishstock
  )
)


matching_summary <- dat |>
  summarise(
    Landings_kg = sum(
      CLscientificWeight
    ),
    Rows = n(),
    .by = c(
      CLyear,
      match_status
    )
  )


assignment_summary <- dat |>
  filter(
    match_status ==
      "unique",
    !is.na(
      fishstock
    )
  ) |>
  summarise(
    Assigned_kg = sum(
      CLscientificWeight
    ),
    .by = c(
      CLyear,
      fishstock,
      special_rule
    )
  )


unresolved_landings <- dat |>
  filter(
    match_status !=
      "unique"
  )


write.csv(
  matching_summary,
  file.path(
    diagnostic_dir,
    "matching_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  assignment_summary,
  file.path(
    diagnostic_dir,
    "stock_assignment_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  unresolved_landings,
  file.path(
    diagnostic_dir,
    "unresolved_landings.csv"
  ),
  row.names = FALSE
)


# Preserve the checkpoint contract expected by model_2.

message(
  "[8/8] Saving rectangle-first stock-matching checkpoint"
)


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
    paste(
      "Checkpoint:",
      normalizePath(
        checkpoint_file
      )
    ),

    paste(
      "Fixed SID reference year:",
      reference_year
    ),

    paste(
      "RDBES years:",
      paste(
        analysed_years,
        collapse = ", "
      )
    ),

    paste(
      "Supported ecoregions:",
      paste(
        supported_ecoregions,
        collapse = "; "
      )
    ),

    paste(
      "EXPERIMENTAL matching method:",
      "rectangle-first using the official 2017 static ICES",
      "statistical-rectangle -> ICES-area map."
    ),

    paste(
      "Ordinary stocks are matched by species plus hierarchical",
      "ICES-area compatibility; ordinary stock polygons are not built."
    ),

    paste(
      "Reported CLarea disambiguates boundary rectangles and is used",
      "when a rectangle is missing from the static map."
    ),

    paste(
      "If rectangle and CLarea conflict, the location is left",
      "geographically unresolved rather than guessed."
    ),

    paste(
      "Nephrops retains geometry because Functional Units are not",
      "ordinary ICES areas; rectangle geometry is compared directly",
      "with FU polygons."
    ),

    "The accepted five-rectangle FU25 correction is retained.",

    paste(
      "Same-taxon alternative stocks remain available as competitors",
      "during matching."
    ),

    paste(
      "The conservative needs_review list and sandeel review treatment",
      "are retained."
    ),

    paste(
      "Sardine PIL records in rectangles 25E4/25E5 assigned to",
      "pil.27.7 are reassigned to pil.27.8abd as a production rule."
    ),

    paste(
      "The sardine country condition cannot be verified because",
      "country is not retained in the compact matching table."
    ),

    paste(
      "No ecoregion share is calculated in this script;",
      "model_2 is unchanged."
    )
  ),

  file.path(
    diagnostic_dir,
    "RUN_NOTES.txt"
  )
)


print(
  matching_summary
)


message(
  "Saved rectangle-first stock-matching checkpoint: ",
  checkpoint_file
)