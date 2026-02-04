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
#> 
#> for wkt geometry draw it on this webapp: https://wktmap.com, and copy the polygon
#> into your location parameters
#> 
#> usethis::edit_r_environ()
#> input GBIF login information
#> Edit your .Renviron to look like this:
#> GBIF_USER="jwaller"
#> GBIF_PWD="safe_fake_password_123"
#> GBIF_EMAIL="jwaller@gbif.org"
#> 
#> devtools::install_github("james-thorson/FishLife") # not available for this version of R
#> devtools::install_github("cfree14/freeR")
#> crul::set_opts(http_version = 2)
#> https://github.com/ropensci/crul/issues/174#issuecomment-1599112561
#> options(timeout=1500)
#> https://stackoverflow.com/questions/35282928/how-do-i-set-a-timeout-for-utilsdownload-file-in-r 
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
gbif_outputname <- 
  "GBIF_Species"
# ------------------- Parameter setting -------------------

if(is.na(worms_taxa[1])) {
  worms_taxa <- obis_taxa
} else {
  print("Using worms_taxa")
}

#GBIF----

## Download GBIF dataset of occurences by geometry ####

# Add just backbone key for a certain taxon, like fish or crustaceans

#### Find GBIF keys for taxa of interest ####
#"Agnatha", "Chondrichthyes", "Osteichthyes"
#"Osteichtheyes" has a wormsID but no gbif id, so using Actinopterygii instead, which is why worms_taxa is an option

#Find list of Classes within 
#Worms (World Register of Marine Species)

classes <- c()

for (i in worms_taxa) {
  wid <- worrms::wm_name2id(name = i)
  classlist <- taxize::worms_downstream(id = wid, downto = "class")
  classes <- rbind(classes, classlist)
}

class_names<-classes$name

backbone_keys <- class_names %>% 
  name_backbone_checklist() %>% # match to backbone 
  filter(!matchType == "NONE") %>% # get matched names
  pull(usageKey) 

#### Download Galapagos Fish Checklist from GBIF ####
download <- occ_download(
   pred_within("POLYGON ((-93.339844 -3.162456, -93.339844 2.547988, -87.1875 2.547988, -87.1875 -3.162456, -93.339844 -3.162456))"),
   pred_in("taxonKey", backbone_keys), # important to use pred_in
   pred("hasCoordinate", TRUE),
   pred("hasGeospatialIssue", FALSE),
   format = "SPECIES_LIST"
 )

download

# Get download ID
download_output <- capture.output(download)
print(download_output)
download_id <- str_sub(download_output[15], start = 17, end = 41)
print(download_id) 

# Source - https://stackoverflow.com/a/55851721
# Posted by sckott
# Retrieved 2026-02-03, License - CC BY-SA 4.0
# Theoretically pauses the code here until the download is done but it is not working
still_running <- TRUE
status_ping <- 9
while (still_running) {
  print(occ_download_meta(key = download_id))
  meta <- occ_download_meta(key = download_id)
  print(meta)
  print(meta$status)
  status <- meta$status
  print(status)
  still_running <- !(status %in% c("SUCCEEDED", "KILLED"))
  print(still_running)
  Sys.sleep(status_ping) # sleep between pings
}

##### Import Galapagos Fish and get the species list ####

GBIF_list <- occ_download_get(download_id) %>%
  occ_download_import()

system("say O C C download complete")

GBIF_species <- GBIF_list %>% 
  filter(taxonRank == "SPECIES") %>%
  select(species) %>%
  rename(Species = species) %>% 
  mutate(Source = "GBIF")

#get list including higher taxonomic levels to compare the results file with
GBIF_taxa1 <- GBIF_list %>% 
  filter(taxonRank == "SPECIES") %>%
  select(species) %>%
  rename(taxa = species)
GBIF_taxa2 <- GBIF_list %>% 
  filter(taxonRank == "GENUS") %>%
  select(genus) %>%
  rename(taxa = genus)
GBIF_taxa3 <- GBIF_list %>% 
  filter(taxonRank == "FAMILY") %>%
  select(family) %>%
  rename(taxa = family)
GBIF_taxa4 <- GBIF_list%>% 
  filter(taxonRank == "ORDER") %>%
  select(order) %>%
  rename(taxa = order)
GBIF_taxa5 <- GBIF_list%>% 
  filter(taxonRank == "CLASS") %>%
  select(class) %>%
  rename(taxa = class)

GBIF_taxa <- unique(rbind(GBIF_taxa1,
                        GBIF_taxa2,
                        GBIF_taxa3,
                        GBIF_taxa4,
                        GBIF_taxa5)) %>% 
  mutate(Source = "GBIF")

write.csv(GBIF_species, paste(here("datasets"), sep = "", "/", gbif_outputname, ".csv"), row.names = F)