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

#> this takes inputnames <- [a vector of csv files containing Genus and Species as separate columns] and  
#> outputname as a string which your output file will be named, and it returns a file named outputname.csv
#> and containing the species latin names as one column and the source as "Local_csv" as another column.
#> This output can then be input into Dataset_combine.


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
obis_taxa <- # works
  c("Agnatha", "Chondrichthyes", "Osteichthyes") 
worms_taxa <- # used because Osteichthyes doesn't have a worms id and Actinopterygii does, TODO fix
  c("Agnatha", "Chondrichthyes", "Actinopterygii") 
taxa2 <-
  "Crustacea"
regional_poly <- 
  "POLYGON ((-117.421875 31.952162, -91.933594 -6.315299, -81.386719 -6.315299, -76.113281 7.710992, -82.089844 8.581021, -87.011719 13.581921, -104.238281 20.303418, -112.5 32.249974, -117.421875 31.952162))"
loc_csvs <- 
  c("2016Aug24_Tirado-Sanchez_et_al_Galapagos_Pisces_Checklist.csv")
loc_outputname <-
  "Local_species"
# ------------------- Parameter setting -------------------

if(is.na(worms_taxa[1])) {
  worms_taxa <- obis_taxa
} else {
  print("Using worms_taxa")
}


#### Galapagos Marine Fish List 2016 ####
x <- 1 # for naming database csv files
Local_species <- data.frame() # for putting species into 
Local_taxa <- data.frame()

for (i in loc_csvs) {
  print(paste("Processing", i))
  db <- data.frame(read.csv(here("datasets", i),
                 fileEncoding = "latin1", 
                 check.names = FALSE))  
  
  # Combine Genus and species columns when species is not just sp. and filter everything else
  sp <- db %>%
    filter(`Species`!= "sp.") %>%
    unite(Genus_species, c("Genus", "Species"), sep =" ")%>%
    select(Genus_species)%>%
    rename(Species=Genus_species) %>% 
    mutate(Source = "Local_csv")
  
  # Get names at all taxonomic levels
  taxa1 <- sp %>% 
    select(Species) %>% 
    rename(taxa = Species)
  taxa2 <- db %>% 
    select(Genus) %>%
    rename(taxa = Genus)
  taxa3 <- db %>% 
    select(Family) %>%
    rename(taxa=Family)
  taxa4 <- db %>% 
    select(Order) %>%
    rename(taxa=Order)
  taxa5 <-db %>% 
    select(Class) %>%
    rename(taxa=Class)
  local_taxaComb <- unique(rbind(taxa1,
                                  taxa2,
                                  taxa3,
                                  taxa4,
                                  taxa5)) %>% 
    mutate(Source = "Local_csv")
  
  Local_species <- rbind(Local_species, sp)
  Local_taxa <- rbind(Local_taxa, local_taxaComb)
  
  x <- x+1
}

write.csv(Local_species, paste(here("datasets"), sep = "", "/", loc_outputname, ".csv"), row.names = F)
print("Local species csv processing complete. Check your datasets folder for the output.")