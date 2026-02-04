# Eldridge Wisely
# This copy modified for function use by Ella Crotty

library(usethis)
library(rgbif)
library(rfishbase)
library(dplyr)
library(readr)  
library(worrms)
library(taxize)
library(robis)
library(tidyverse)
library(taxonomizr)
library(readr)
library(here)


#https://docs.ropensci.org/rgbif/articles/gbif_credentials.html
#usethis::edit_r_environ()

# ------------------- Setup Instructions -------------------
#> Make an R project and a folder called "databases" in the 
#> working directory, this will make the 
#> here() commands work
#> make a folder called custom_db
#> for wkt geometry draw it on this webapp: https://wktmap.com, and copy the polygon
#> into your location parameters
#> usethis::edit_r_environ()
#> input GBIF login information
#> devtools::install_github("james-thorson/FishLife")
#> devtools::install_github("cfree14/freeR")
# ------------------- Setup Instructions -------------------

# ------------------- Parameter setting -------------------
darwin_fish_db <- # works
  "https://datazone.darwinfoundation.org/media/pdf/checklist/2016Aug24_Tirado-Sanchez_et_al_Galapagos_Pisces_Checklist.csv"
darwin_crust_db <- # works
  "https://datazone.darwinfoundation.org/media/pdf/checklist/2016Sep30_Tirado-Sanchez_et_al_Galapagos_Marine_crustaceans_Checklist.csv"
obis_taxa <- # works
  c("Agnatha", "Chondrichthyes", "Osteichthyes") 
worms_taxa <- # used because Osteichthyes doesn't have a worms id and Actinopterygii does, TODO fix
  c("Agnatha", "Chondrichthyes", "Actinopterygii") 
taxa2 <-
  "Crustacea"
local_poly <- "POLYGON ((-93.339844 -3.162456, -93.339844 2.547988, -87.1875 2.547988, -87.1875 -3.162456, -93.339844 -3.162456))"
regional_poly <- "POLYGON ((-117.421875 31.952162, -91.933594 -6.315299, -81.386719 -6.315299, -76.113281 7.710992, -82.089844 8.581021, -87.011719 13.581921, -104.238281 20.303418, -112.5 32.249974, -117.421875 31.952162))"
# ------------------- Parameter setting -------------------

if(is.na(worms_taxa[1])) {
  worms_taxa <- obis_taxa
} else {
  print("Using worms_taxa")
}

#OBIS----
####Galapagos checklist####
obis_local <- checklist(obis_taxa, 
                        geometry = local_poly)

#> # Source - https://stackoverflow.com/a/58737047
# Posted by R me matey, modified by community. See post 'Timeline' for change history
# Retrieved 2026-02-03, License - CC BY-SA 4.0

#library(dplyr)
#y <- ""
#data.frame(x = 1:5) %>% 
#  filter(case_when(y=="" ~ x > 3, #When y == "", x > 3
#                   T ~ x<3) #Otherwise, x < 3
#  ) %>% 
#  tail(1)


#### Pull out only the marine species list for the local database
obis_local_species <- obis_local %>% # this is not working - it removes all terrestrial but 
  # still gets some that are freshwater
  filter(case_when(
    !is.na(marine) ~ is_marine == marine,
    !is.na(terrestrial) ~ is_terrestrial == terrestrial, # this is outranking everything below it somehow
    !is.na(freshwater) ~ is_freshwater == freshwater,
    !is.na(brackish) ~ is_brackish == brackish,
    .default = taxonRank == "Species"
  )) %>% 
  select(scientificName) %>%
  rename(Species = scientificName) %>% 
  mutate(Source = "OBIS")

#### Get the local marine taxa list for checking the results against
obis_local_taxa <- obis_local %>%
  filter(case_when(
    !is.na(marine) ~ is_marine == marine,
    !is.na(terrestrial) ~ is_terrestrial == terrestrial,
    !is.na(freshwater) ~ is_freshwater == freshwater,
    !is.na(brackish) ~ is_brackish == brackish,
    .default = taxonRank == "Species"
  )) %>%
  select(scientificName) %>%
  rename(Taxa=scientificName)

####TEP checklist####
obis_tep <- checklist(obis_taxa,
                      geometry = regional_poly)
# Much larger!

#### Get the regional marine taxa list for checking the results against
obis_tep_taxa <- obis_tep %>%
  filter(case_when(
    !is.na(marine) ~ is_marine == marine,
    !is.na(terrestrial) ~ is_terrestrial == terrestrial,
    !is.na(freshwater) ~ is_freshwater == freshwater,
    !is.na(brackish) ~ is_brackish == brackish,
    .default = taxonRank == "Species"
  )) %>%
  select(scientificName) %>%
  rename(Taxa = scientificName)
