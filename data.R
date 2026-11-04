# Prepare data, write CSV data tables

# Before:
# After:

library(icesTAF)
library(dplyr)
library(stringr)
library(tidyr)
library(icesVocab)

mkdir("data")



# Import data---------------------------
RDBESwd <- "./boot/initial/data/"

# # sspCorresp <- data.table::fread(paste0("./bootstrap/data/", "RDBES_especes_fao.csv")) %>% 
# sspCorresp <- data.table::fread(paste0(RDBESwd, "RDBES_especes_fao.csv")) %>% 
#   select(ESP_COD_FAO, CLspeciesCode = ICES_APHIA) 

# ### !! some duplicated CLspeciesCode
# sspCorresp %>% group_by(CLspeciesCode) %>% summarise(n= n()) %>% filter(n>1)
# ### fix, take the first one not to duplicate cl latter
# sspCorresp <- sspCorresp %>% distinct(CLspeciesCode, .keep_all = TRUE)

# Species reference -------------------------------------------------------

# Current official ASFIS table prepared in boot
asfis_species <- data.table::fread(
  file.path("./boot/ASFIS_species.csv")
) %>%
  transmute(
    Alpha3_Code = toupper(trimws(as.character(Alpha3_Code))),
    Scientific_Name,
    English_name
  ) %>%
  filter(!is.na(Alpha3_Code), Alpha3_Code != "") %>%
  distinct(Alpha3_Code, .keep_all = TRUE)


# ICES Vocabulary crosswalk: ASFIS code <-> WoRMS/Aphia ID
asfis_worms <- icesVocab::getCodeTypeRelation(
  "SpecASFIS",
  "SpecWoRMS"
) %>%
  transmute(
    Alpha3_Code = toupper(trimws(as.character(SpecASFIS))),
    CLspeciesCode = as.character(SpecWoRMS)
  ) %>%
  filter(
    !is.na(Alpha3_Code),
    !is.na(CLspeciesCode)
  ) %>%
  # Only retain codes that exist in the current ASFIS table
  semi_join(asfis_species, by = "Alpha3_Code") %>%
  distinct()


# Identify Aphia IDs that map to exactly one current ASFIS code
species_crosswalk <- asfis_worms %>%
  group_by(CLspeciesCode) %>%
  summarise(
    n_asfis_codes = n_distinct(Alpha3_Code),
    Alpha3_Code = if_else(
      n_asfis_codes == 1L,
      first(Alpha3_Code),
      NA_character_
    ),
    .groups = "drop"
  )


# Keep ambiguous mappings visible for auditing
ambiguous_species_crosswalk <- species_crosswalk %>%
  filter(n_asfis_codes > 1)

write.csv(
  ambiguous_species_crosswalk,
  "./data/species_crosswalk_ambiguous.csv",
  row.names = FALSE
)

### import CL and CE data
cl.all <- data.table::fread(file = paste0(RDBESwd, "/RDBES Landings CL 2021 2025/CommercialLanding.csv"), quote = "")

unique(cl.all$CLyear)

# ce.all <- read.csv(file = paste0(RDBESwd, "/RDBES Effort CE 2021 2025/CommercialEffort.csv"), quote = "")
# dim(cl.all); dim(ce.all)


### add CLspeciesFaoCode when missing in cl

# cl_old <- left_join(cl.all, sspCorresp)

# cl_old <- cl_old %>%
#   mutate(CLspeciesFaoCode = case_when(
#     CLspeciesFaoCode=="" ~ ESP_COD_FAO,
#     TRUE ~ CLspeciesFaoCode
#   ))
# Fill missing FAO species codes -----------------------------------------



###################### new part
cl.all <- cl.all %>%
  mutate(
    CLspeciesCode = as.character(CLspeciesCode),
    CLspeciesFaoCode = na_if(
      toupper(trimws(as.character(CLspeciesFaoCode))),
      ""
    )
  )

n_rows_before_species_join <- nrow(cl.all)

cl.all <- cl.all %>%
  left_join(
    species_crosswalk %>%
      select(
        CLspeciesCode,
        ASFIS_from_Aphia = Alpha3_Code
      ),
    by = "CLspeciesCode",
    relationship = "many-to-one"
  ) %>%
  mutate(
    # Never replace an existing RDBES FAO code.
    # Use the Aphia-derived code only when the mapping was unique.
    CLspeciesFaoCode = coalesce(
      CLspeciesFaoCode,
      ASFIS_from_Aphia
    )
  )

# The reference join must never duplicate RDBES landing records
stopifnot(nrow(cl.all) == n_rows_before_species_join)


# create a database splitting lines where several vessels have been input in one strata. The hypothesis here that we divide the values (landins and effort) by the number of vessels

## cl--------------
cl.all <- cl.all %>%
  mutate(CLlandingsValue = as.numeric(CLlandingsValue))


specials <- str_extract_all(cl.all$CLencryptedVesselIds, "[^A-Za-z0-9]")
# Pour obtenir tous les caractères spéciaux uniques :
unique(unlist(specials))



#### differnet special characters per countries...
## FI EI LT DE GB-ENG GB-WLS GB-NIR EE SE FR PL NL DK PT are using ";" as a separator for CLencryptedVesselIds
## ES uses "/" as a separator for CLencryptedVesselIds
## LV uses "-" as a separator for CLencryptedVesselIds


cl.all2 <- cl.all  %>%
  mutate(numVessel = str_count(CLencryptedVesselIds, ";|/|-")) %>%
  mutate(numVessel = case_when(
    numVessel > 0 ~ numVessel + 1,
    .default = 1
  )) %>%
  mutate(
    CLscientificWeight_perVessel = CLscientificWeight / numVessel, 
    CLlandingsValue_perVessel = CLlandingsValue / numVessel
  ) %>%
  mutate(CLencryptedVesselIds = strsplit(CLencryptedVesselIds, ";|/|-")) %>%
  unnest(CLencryptedVesselIds) %>%
  select(-c("CLscientificWeight", "CLlandingsValue", "numVessel")) %>%
  rename(CLscientificWeight = CLscientificWeight_perVessel, 
         CLlandingsValue = CLlandingsValue_perVessel)

# checks
dim(cl.all); dim(cl.all2)
cl.all  %>% summarise((sum(CLlandingsValue, na.rm=T))) / cl.all2 %>% summarise((sum(CLlandingsValue, na.rm=T)))
cl.all  %>% summarise((sum(CLscientificWeight, na.rm=T))) / cl.all2 %>% summarise((sum(CLscientificWeight, na.rm=T)))


save(cl.all2, file = "./data/extendedRDBES.RData")
