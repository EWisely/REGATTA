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

# read.csv(outputs)
# ------------------- Parameter setting -------------------

if(is.na(worms_taxa[1])) {
  worms_taxa <- obis_taxa
} else {
  print("Using worms_taxa")
}

#Combine OBIS, GBIF and Darwin lists----

Species_local_comb <- unique(rbind(obis_local_species,
                                   #GBIF_species, 
                                   Darwin_species))
Species_local_comb <- Species_local_comb %>%
  filter(Species!="NA")

#order alphabetically
Species_local_comb$Species <- 
  Species_local_comb$Species[order(Species_local_comb$Species)]

##Write  list to a txt file to use as CRABS input----
write_delim(Species_local_comb, 
            paste(here(),
                  sep = "",
                  "/custom_db/comprehensive_sp_list.txt"), 
            delim = '\t', 
            col_names = FALSE)
nrow(Species_local_comb)

print("Insert CRABS instructions here")
