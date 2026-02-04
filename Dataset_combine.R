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
regional_poly <- 
  "POLYGON ((-117.421875 31.952162, -91.933594 -6.315299, -81.386719 -6.315299, -76.113281 7.710992, -82.089844 8.581021, -87.011719 13.581921, -104.238281 20.303418, -112.5 32.249974, -117.421875 31.952162))"

comb_inputnames <- c("Local_species", "OBIS_Species", "GBIF_Species") # allows for running the input functions multiple times
# if desired
comb_outputname <- "Comprehensive_species_list"
# ------------------- Parameter setting -------------------

Species_local_comb <- data.frame()

for (i in comb_inputnames) {
  db <- read.csv(paste(here("datasets", i), 
                       sep = "", ".csv"), 
                 fileEncoding = "latin1", 
                 check.names = FALSE)
  Species_local_comb <- unique(rbind(Species_local_comb, db))
}

Species_local_comb <- Species_local_comb %>%
  filter(Species!="NA")

##Write  list to a txt file to use as CRABS input----
write_delim(Species_local_comb, 
            paste(here("custom_db"),
                  sep = "",
                  "/",
                  comb_outputname,
                  ".txt"), 
            delim = '\t')

print("Insert CRABS instructions here")
