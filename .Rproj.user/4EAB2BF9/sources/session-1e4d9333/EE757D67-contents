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

#GBIF----

## Download GBIF dataset of occurences within the Galapagos region by geometry ####

#All Galapagos Species (not just fish and crustaceans)

# Add just backbone key for fish or crustaceans!


#### Find GBIF keys for crustaceans ####
#Find list of Classes within Crustacea subphylum
#Worms (World Register of Marine Species)

worrms_id <- worrms::wm_name2id(name = taxa2)
#1066
worrms::wm_external(id = worrms_id, type = "ncbi")
#6657

taxize::worms_downstream(id = worrms_id, downto = "class")
#id          name  rank
#1    1070   Maxillopoda class
#2    1069  Branchiopoda class

crustacean_classes <- 
  taxize::worms_downstream(id = worrms_id, downto = "class")
crustacean_class_names <- 
  crustacean_classes$name

Crustacean_backbone_keys <- crustacean_class_names %>% 
  name_backbone_checklist() %>% # match to backbone 
  filter(!matchType == "NONE") %>% # get matched names
  pull(usageKey) 

#### Download Galapagos Crustacean list ####

# occ_download(
#   pred_within("POLYGON ((-93.339844 -3.162456, -93.339844 2.547988, -87.1875 2.547988, -87.1875 -3.162456, -93.339844 -3.162456))"),
#   pred_in("taxonKey", Crustacean_backbone_keys), # important to use pred_in
#   pred("hasCoordinate", TRUE),
#   pred("hasGeospatialIssue", FALSE),
#   format = "SPECIES_LIST"
# )

##### Download TEP Crustacean list #####


# occ_download(
#   pred_within("POLYGON ((-117.421875 31.952162, -91.933594 -6.315299, -81.386719 -6.315299, -76.113281 7.710992, -82.089844 8.581021, -87.011719 13.581921, -104.238281 20.303418, -112.5 32.249974, -117.421875 31.952162))"),
#   pred_in("taxonKey", Crustacean_backbone_keys), # important to use pred_in
#   pred("hasCoordinate", TRUE),
#   pred("hasGeospatialIssue", FALSE),
#   format = "SPECIES_LIST"
# )

#### Find GBIF keys for fish ####
#"Agnatha", "Chondrichthyes", "Osteichthyes"
#"Osteichtheyes" has a wormsID but no gbif id, so using Actinopterygii instead.

#Find list of Classes within 
#Worms (World Register of Marine Species)

fish_classes <- c()

for (i in worms_taxa) {
  wid <- worrms::wm_name2id(name = i)
  classlist <- taxize::worms_downstream(id = wid, downto = "class")
  fish_classes <- rbind(fish_classes, classlist)
}

fish_class_names<-fish_classes$name

Fish_backbone_keys <- fish_class_names %>% 
  name_backbone_checklist() %>% # match to backbone 
  filter(!matchType == "NONE") %>% # get matched names
  pull(usageKey) 

#### Download Galapagos Fish Checklist from GBIF ####
# occ_download(
#   pred_within("POLYGON ((-93.339844 -3.162456, -93.339844 2.547988, -87.1875 2.547988, -87.1875 -3.162456, -93.339844 -3.162456))"),
#   pred_in("taxonKey", Fish_backbone_keys), # important to use pred_in
#   pred("hasCoordinate", TRUE),
#   pred("hasGeospatialIssue", FALSE),
#   format = "SPECIES_LIST"
# )


##### Download TEP Fish list #####


# occ_download(
#   pred_within("POLYGON ((-117.421875 31.952162, -91.933594 -6.315299, -81.386719 -6.315299, -76.113281 7.710992, -82.089844 8.581021, -87.011719 13.581921, -104.238281 20.303418, -112.5 32.249974, -117.421875 31.952162))"),
#   pred_in("taxonKey", Fish_backbone_keys), # important to use pred_in
#   pred("hasCoo33rdinate", TRUE),
#   pred("hasGeospatialIssue", FALSE),
#   format = "SPECIES_LIST"
# )

# <<gbif download>>
# Your download is being processed by GBIF:
#   https://www.gbif.org/occurrence/download/0172907-240321170329656
# After it finishes, use
# d <- occ_download_get('0172907-240321170329656') %>%
#   occ_download_import()
# to retrieve your download.

##### Import Galapagos Crustaceans and get the species list ####

GBIF_Gal_crustaceans_list <- occ_download_get('0172650-240321170329656') %>%
  occ_download_import()

GBIF_Gal_crustaceans_species <- GBIF_Gal_crustaceans_list%>% 
  filter(taxonRank=="SPECIES")%>%
  select(species)%>%
  rename(Species=species)


#get list including higher taxonomic levels to compare the results file with
GBIF_Gal_crustaceans_taxa1<-GBIF_Gal_crustaceans_list%>% 
  filter(taxonRank=="SPECIES")%>%
  select(species)%>%
  rename(taxa=species)
GBIF_Gal_crustaceans_taxa2<-GBIF_Gal_crustaceans_list%>% 
  filter(taxonRank=="GENUS")%>%
  select(genus)%>%
  rename(taxa=genus)
GBIF_Gal_crustaceans_taxa3<-GBIF_Gal_crustaceans_list%>% 
  filter(taxonRank=="FAMILY")%>%
  select(family)%>%
  rename(taxa=family)
GBIF_Gal_crustaceans_taxa4<-GBIF_Gal_crustaceans_list%>% 
  filter(taxonRank=="ORDER")%>%
  select(order)%>%
  rename(taxa=order)
GBIF_Gal_crustaceans_taxa5<-GBIF_Gal_crustaceans_list%>% 
  filter(taxonRank=="CLASS")%>%
  select(class)%>%
  rename(taxa=class)
GBIF_Gal_crustaceans_taxa<-unique(rbind(GBIF_Gal_crustaceans_taxa1,GBIF_Gal_crustaceans_taxa2,GBIF_Gal_crustaceans_taxa3,GBIF_Gal_crustaceans_taxa4, GBIF_Gal_crustaceans_taxa5))


##### Import TEP Crustaceans and get the taxa list ####

GBIF_TEP_crustaceans_list <- occ_download_get('0172684-240321170329656') %>%
  occ_download_import()

GBIF_TEP_crustaceans_taxa1<-GBIF_TEP_crustaceans_list%>% 
  filter(taxonRank=="SPECIES")%>%
  select(species)%>%
  rename(taxa=species)
GBIF_TEP_crustaceans_taxa2<-GBIF_TEP_crustaceans_list%>% 
  filter(taxonRank=="GENUS")%>%
  select(genus)%>%
  rename(taxa=genus)
GBIF_TEP_crustaceans_taxa3<-GBIF_TEP_crustaceans_list%>% 
  filter(taxonRank=="FAMILY")%>%
  select(family)%>%
  rename(taxa=family)
GBIF_TEP_crustaceans_taxa4<-GBIF_TEP_crustaceans_list%>% 
  filter(taxonRank=="ORDER")%>%
  select(order)%>%
  rename(taxa=order)
GBIF_TEP_crustaceans_taxa5<-GBIF_TEP_crustaceans_list%>% 
  filter(taxonRank=="CLASS")%>%
  select(class)%>%
  rename(taxa=class)
GBIF_TEP_crustaceans_taxa<-unique(rbind(GBIF_TEP_crustaceans_taxa1,GBIF_TEP_crustaceans_taxa2,GBIF_TEP_crustaceans_taxa3,GBIF_TEP_crustaceans_taxa4, GBIF_TEP_crustaceans_taxa5))


##### Import Galapagos Fish and get the species list ####

GBIF_Gal_fish_list <- occ_download_get('0172899-240321170329656') %>%
  occ_download_import()

GBIF_Gal_fish_species<-GBIF_Gal_fish_list%>% 
  filter(taxonRank=="SPECIES")%>%
  select(species)%>%
  rename(Species=species)

#get list including higher taxonomic levels to compare the results file with
GBIF_Gal_fish_taxa1<-GBIF_Gal_fish_list%>% 
  filter(taxonRank=="SPECIES")%>%
  select(species)%>%
  rename(taxa=species)
GBIF_Gal_fish_taxa2<-GBIF_Gal_fish_list%>% 
  filter(taxonRank=="GENUS")%>%
  select(genus)%>%
  rename(taxa=genus)
GBIF_Gal_fish_taxa3<-GBIF_Gal_fish_list%>% 
  filter(taxonRank=="FAMILY")%>%
  select(family)%>%
  rename(taxa=family)
GBIF_Gal_fish_taxa4<-GBIF_Gal_fish_list%>% 
  filter(taxonRank=="ORDER")%>%
  select(order)%>%
  rename(taxa=order)
GBIF_Gal_fish_taxa5<-GBIF_Gal_fish_list%>% 
  filter(taxonRank=="CLASS")%>%
  select(class)%>%
  rename(taxa=class)
GBIF_Gal_fish_taxa<-unique(rbind(GBIF_Gal_fish_taxa1,GBIF_Gal_fish_taxa2,GBIF_Gal_fish_taxa3,GBIF_Gal_fish_taxa4, GBIF_Gal_fish_taxa5))


##### Import TEP Fish and get the taxa list ####
GBIF_TEP_fish_list<-occ_download_get('0172907-240321170329656') %>%
  occ_download_import()

GBIF_TEP_fish_taxa1<-GBIF_TEP_fish_list%>% 
  filter(taxonRank=="SPECIES")%>%
  select(species)%>%
  rename(taxa=species)
GBIF_TEP_fish_taxa2<-GBIF_TEP_fish_list%>% 
  filter(taxonRank=="GENUS")%>%
  select(genus)%>%
  rename(taxa=genus)
GBIF_TEP_fish_taxa3<-GBIF_TEP_fish_list%>% 
  filter(taxonRank=="FAMILY")%>%
  select(family)%>%
  rename(taxa=family)
GBIF_TEP_fish_taxa4<-GBIF_TEP_fish_list%>% 
  filter(taxonRank=="ORDER")%>%
  select(order)%>%
  rename(taxa=order)
GBIF_TEP_fish_taxa5<-GBIF_TEP_fish_list%>% 
  filter(taxonRank=="CLASS")%>%
  select(class)%>%
  rename(taxa=class)
GBIF_TEP_fish_taxa<-unique(rbind(GBIF_TEP_fish_taxa1,GBIF_TEP_fish_taxa2,GBIF_TEP_fish_taxa3,GBIF_TEP_fish_taxa4, GBIF_TEP_fish_taxa5))
