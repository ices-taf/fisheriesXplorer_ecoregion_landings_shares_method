#' Load ASFIS species reference table
#'
#' Reads the local ASFIS ASFIS\_sp CSV file and returns a data frame with
#' English common name, scientific name, and FAO three-letter species code.
#'
#' @details
#' This function expects the file `./data-raw/ASFIS_sp_2025.csv` to be present
#' relative to the project root (or package root, if used inside a package).
#' Only the columns `English_name`, `Scientific_Name`, and `Alpha3_Code`
#' are kept.
#'
#' @return
#' A data frame (or tibble, if used with dplyr) with three character columns:
#' \describe{
#'   \item{English_name}{English common name of the species.}
#'   \item{Scientific_Name}{Scientific (Latin) name of the species.}
#'   \item{Alpha3_Code}{FAO three-letter species code.}
#' }
#'
#' @examples
#' \dontrun{
#' species <- load_asfis_species()
#' head(species)
#' }
#'
load_asfis_species <- function() {

  # 1. Define the FAO URL
  url <- "http://www.fao.org/fishery/static/ASFIS/ASFIS_sp_2026.1.zip"

  # 2. Create a temporary directory
  temp_dir <- tempfile()
  dir.create(temp_dir)

  # Automatically remove temporary files when the function finishes
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  # 3. Download the ZIP file
  zip_file <- file.path(temp_dir, "ASFIS_sp_2026.1.zip")

  utils::download.file(
    url = url,
    destfile = zip_file,
    mode = "wb",
    quiet = TRUE
  )

  # 4. Identify the CSV file inside the ZIP
  zip_contents <- utils::unzip(zip_file, list = TRUE)$Name

  csv_file <- zip_contents[
    grepl("ASFIS_sp_2026\\.1\\.csv$", zip_contents)
  ]

  if (length(csv_file) != 1) {
    stop("Could not identify the ASFIS CSV file in the ZIP.")
  }

  # 5. Extract the CSV file
  utils::unzip(
    zipfile = zip_file,
    files = csv_file,
    exdir = temp_dir
  )

  # 6. Read the CSV (original workflow)
  species <- read.csv(
    file.path(temp_dir, csv_file),
    na.strings = "",
    stringsAsFactors = FALSE
  )

  # 7. Keep only the required columns (original workflow)
  species <- dplyr::select(
    species,
    English_name,
    Scientific_Name,
    Alpha3_Code
  )

  return(species)
}


# 8. Execute function and save TAF boot output
species <- load_asfis_species()

write.csv(
  species,
  file = "./boot/ASFIS_species.csv",
  row.names = FALSE,
  na = ""
)
