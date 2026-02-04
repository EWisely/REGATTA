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

#Darwin foundation ----

#### Galapagos Marine Fish List 2016 ####
x <- 1 # for naming database csv files
Darwin_species <- data.frame() # for putting species into 
Darwin_taxa <- data.frame()

for (i in darwin_dbs) {
  print(i)
  curl::curl_download(url = i, 
                      destfile = paste(here("databases"), 
                                       sep = "",
                                       "/", 
                                       "darwin_checklist_", 
                                       x, 
                                       ".csv"))
  darwin_db <- read.csv(paste(here("databases"), 
                              sep = "",
                              "/", 
                              "darwin_checklist_", 
                              x, 
                              ".csv"), 
                        fileEncoding = "latin1", 
                        check.names = FALSE)
  
  # Combine Genus and species columns when species is not just sp. and filter everything else
  darwin_sp <- darwin_db %>%
    filter(`Specific Epithtet`!= "sp.")%>%
    unite(Genus_species, c("Genus", "Specific Epithtet"), sep =" ")%>%
    select(Genus_species)%>%
    rename(Species=Genus_species) %>% 
    mutate(Source = "Darwin")
  
  # Get names at all taxonomic levels
  Darwin_taxa1 <- darwin_sp %>% 
    select(Species) %>% 
    rename(taxa = Species)
  Darwin_taxa2 <- darwin_db %>% 
    select(Genus) %>%
    rename(taxa = Genus)
  Darwin_taxa3 <- darwin_db %>% 
    select(Family) %>%
    rename(taxa=Family)
  Darwin_taxa4 <- darwin_db %>% 
    select(Order) %>%
    rename(taxa=Order)
  Darwin_taxa5<-darwin_db %>% 
    select(Class) %>%
    rename(taxa=Class)
  Darwin_taxaComb <- unique(rbind(Darwin_taxa1,
                                  Darwin_taxa2,
                                  Darwin_taxa3,
                                  Darwin_taxa4,
                                  Darwin_taxa5)) %>% 
    mutate(Source = "Darwin")
  
  Darwin_species <- rbind(Darwin_species, darwin_sp)
  Darwin_taxa <- rbind(Darwin_taxa, Darwin_taxaComb)
  
  x <- x+1
}