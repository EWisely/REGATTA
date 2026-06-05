# Eldridge Wisely
# This copy modified for function use by Ella Crotty

#' Download a GBIF species list inside a WKT polygon
#'
#' Resolves the given high-level taxa via WoRMS, finds GBIF backbone
#' keys for the resulting classes, and submits a GBIF occurrence
#' download restricted to the polygon. Writes the result to
#' `datasets/<gbif_outputname>.csv` and prints the file path. GBIF
#' login credentials must already be set in `~/.Renviron`; see the
#' Setup section of the README.
#'
#' @param obis_taxa A character vector of taxon names at **class level or
#'   broader** (GBIF won't accept order/family/genus/species here).
#' @param worms_taxa A character vector of substitute taxon names to use
#'   instead of `obis_taxa` when looking up WoRMS IDs. NA reuses
#'   `obis_taxa`. Useful when WoRMS labels differ from GBIF backbone
#'   labels (e.g. WoRMS `Osteichthyes` vs. GBIF `Actinopterygii`).
#' @param regional_poly A WKT POLYGON string of the form
#'   `"POLYGON ((long lat, long lat, ...))"`. Draw a region on
#'   [wktmap.com](https://wktmap.com) and copy the generated polygon.
#' @param gbif_outputname Basename (no `.csv` extension) for the output
#'   file under `datasets/`. Default `"GBIF_Species"`.
#'
#' @details
#' Takes ~15 minutes to run end to end. The function prints status
#' ("PREPARING" / "RUNNING") in the console while polling the GBIF
#' download API; this is normal. GBIF coverage is uneven — see the
#' Common troubleshooting section of the README for known issues
#' (Osteichthyes not recognized, occasional timeouts, etc.).
#'
#' @return Invisibly NULL. Writes
#'   `datasets/<gbif_outputname>.csv` as a side effect.
#'
#' @examples
#' \dontrun{
#' GBIF_download(
#'   obis_taxa       = c("Osteichthyes", "Multicrustacea"),
#'   worms_taxa      = c("Actinopterygii", "Multicrustacea"),
#'   regional_poly   = "POLYGON ((-124.85 51.75, -129.20 51.15, -128.14 41.87, -122.34 42.00, -121.86 44.62, -122.70 46.59, -121.55 47.87, -124.85 51.75))",
#'   gbif_outputname = "GBIF_OCNMS_fish_crust"
#' )
#' }
#'
#' @export
GBIF_download <- function(obis_taxa,
                          worms_taxa = NA,
                          regional_poly,
                          gbif_outputname = "GBIF_Species"
    ) {
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
  
  if(is.na(worms_taxa[1])) {
    worms_taxa <- obis_taxa
  } else {
    print("Using worms_taxa")
  }
  
  # throw an error if polygon is incorrect
  if(grepl("POLYGON \\(\\([-. |[[:digit:]]|,]*\\)\\)", # finds POLYGON (( [-, space, period, 
           # comma or any number] continue for any length, ))
           regional_poly)) { # make sure regional_poly is formatted correctly
    corrects <- T
  } else {
    corrects <- F
  }
  ## Download GBIF dataset of occurences by geometry ####
  # Add just backbone key for a certain taxon, like fish or crustaceans
  #### Find GBIF keys for taxa of interest ####
  #"Osteichtheyes" has a wormsID but no gbif id, so using Actinopterygii instead, which is why worms_taxa is an option
  
  #Find list of Classes within 
  #Worms (World Register of Marine Species)
  
  classes <- c()
  
  for (i in worms_taxa) {
    print(paste("Getting WoRMS ID for", i))
    
    tryCatch( { # Test the code first to see if it works
      worrms::wm_name2id(name = i)
    },
    error = function(e) {
      print(e)
      message(paste("WoRMS ID Not Found, Check That Your Listed Taxon is in WoRMS:", i))
    }
    )
    
    wid <- worrms::wm_name2id(name = i)
    print(paste("Getting WoRMS Downstream for", i))
    
    tryCatch( { # Test the code first to see if it works
      taxize::worms_downstream(id = wid, downto = "class")
    },
    error = function(e) {
      print(e)
      message(paste("WoRMS Downstream Not Found, Check That Your Listed Taxon is in WoRMS:", i))
    }
    )
    
    classlist <- taxize::worms_downstream(id = wid, downto = "class")
    print(paste(i, sep = " ", "class list got:"))
    print(head(classlist))
    classes <- rbind(classes, classlist)
  }
  
  class_names<-classes$name
  
  print("Full class list:")
  print(class_names)
  classnameno <- length(class_names)
  print("Getting backbone keys")
  
  backbone_keys <- class_names %>% 
    name_backbone_checklist() %>% # match to backbone 
    filter(!matchType == "NONE") %>% # get matched names
    pull(usageKey) 
  
  print("GBIF backbone keys:")
  print(backbone_keys)
  backbonekeyno <- length(backbone_keys)
  
  print("Downloading Database from GBIF")
  #### Download Checklist from GBIF ####
  download <- occ_download(
    pred_within(regional_poly),
    pred_in("taxonKey", backbone_keys), # important to use pred_in
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE),
    format = "SPECIES_LIST"
  )
  
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
  print(occ_download_meta(key = download_id)) # print info about download
  
  while (still_running) { # starts true
    meta <- occ_download_meta(key = download_id) # get metadata about download
    print(meta$status) # print just the status, the rest of it won't have changed since previous print 
    status <- meta$status 
    still_running <- !(status %in% c("SUCCEEDED", "KILLED")) # SUCCEEDED or KILLED are end states so they stop the loop
    Sys.sleep(status_ping) # sleep between pings
  }
  
  ##### Import and get the species list ####
  
  GBIF_list <- occ_download_get(download_id) %>%
    occ_download_import()
  
  print("Download Complete")
  
  # Reprint all previous output because the download status stuff clogged the console
  print("Full class list:")
  print(class_names)
  print("Getting backbone keys")
  print("GBIF backbone keys:")
  print(backbone_keys)
  
  
  GBIF_species <- GBIF_list %>% 
    filter(taxonRank == "SPECIES") %>%
    select(species) %>%
    rename(Species = species) %>% 
    mutate(Source = "GBIF")
  
  print("Species list/output first few rows:")
  print(head(GBIF_species))
  print("Getting higher taxonomic levels")
  
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
  
  print("First few lines of higher taxonomic levels:")
  print(head(GBIF_taxa))
  
  write.csv(GBIF_species, paste(here("datasets"), sep = "", "/", gbif_outputname, ".csv"), row.names = F)
  
  print("GBIF Download & Export Complete. Check your datasets folder for the output.")
  print(paste("Output at:", sep = " ", paste(here("datasets"), sep = "", "/", gbif_outputname, ".csv")))
  print(paste(backbonekeyno, sep = " ", "backbone keys found out of", classnameno, "classes in list"))
}

