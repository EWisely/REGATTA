# Eldridge Wisely
# This copy modified for function use by Ella Crotty

#' Download an OBIS species list inside a WKT polygon
#'
#' Pulls an OBIS checklist for the given taxa within a regional WKT
#' polygon, with optional habitat filters. Writes the result to
#' `datasets/<obis_outputname>.csv` and prints the file path.
#'
#' @param obis_taxa A character vector of taxon names at any level (OBIS
#'   recognizes class, order, family, genus, species, etc.).
#' @param worms_taxa A character vector of substitute taxon names to use
#'   instead of `obis_taxa` when looking up WoRMS IDs. NA reuses
#'   `obis_taxa`.
#' @param regional_poly A WKT POLYGON string of the form
#'   `"POLYGON ((long lat, long lat, ...))"`. Draw a region on
#'   [wktmap.com](https://wktmap.com) and copy the generated polygon.
#' @param obis_outputname Basename (no `.csv` extension) for the output
#'   file under `datasets/`. Default `"OBIS_Species"`. Choose a short
#'   distinctive name — you will pass it into
#'   [build_regional_checklist()] later.
#' @param marine,freshwater,terrestrial,brackish Habitat filters. TRUE
#'   returns only species marked in that habitat; FALSE excludes them;
#'   NA (default) means no filter.
#'
#' @details
#' The four habitat filters **stack as intersections, not unions**. So
#' `marine = TRUE, freshwater = TRUE` returns only species marked both
#' marine and freshwater (i.e. anadromous). To get marine OR
#' freshwater, set `terrestrial = FALSE` and leave marine/freshwater
#' as NA.
#'
#' @return Invisibly NULL. Writes
#'   `datasets/<obis_outputname>.csv` as a side effect.
#'
#' @examples
#' \dontrun{
#' # Pacific NW (OCNMS) salmonids and copepods, intentionally allowing
#' # salmon to slip through the habitat filter. The same pattern works
#' # for any taxa over any WKT polygon — freshwater insects in
#' # Germany, terrestrial mammals over a Sonoran Desert polygon, etc.
#' # obis_outputname is yours to pick.
#' OBIS_download(
#'   obis_taxa       = c("Salmonidae", "Copepoda"),
#'   regional_poly   = "POLYGON ((-124.85 51.75, -129.20 51.15, -128.14 41.87, -122.34 42.00, -121.86 44.62, -122.70 46.59, -121.55 47.87, -124.85 51.75))",
#'   obis_outputname = "OBIS_my_region_my_group",
#'   marine          = TRUE,
#'   freshwater      = NA,    # needed to keep salmon, which spans habitats
#'   terrestrial     = FALSE,
#'   brackish        = NA
#' )
#' }
#'
#' @export
OBIS_download <- function(obis_taxa,
                          worms_taxa = NA,
                          regional_poly,
                          obis_outputname = "OBIS_Species",
                          marine = NA,
                          freshwater = NA,
                          terrestrial = NA,
                          brackish = NA) {
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
  
  print("WARNING: Check the documentation for setup 
        steps before running this function!")
  
  # throw an error if polygon is incorrect
  if(grepl("POLYGON \\(\\([-. |[[:digit:]]|,]*\\)\\)", # finds POLYGON (( [-, space, period, 
           # comma or any number] continue for any length, ))
           regional_poly)) { # make sure regional_poly is formatted correctly
    corrects <- T
  } else {
    corrects <- F
  }
  if(corrects == F) stop("regional_poly must be in POLYGON format: POLYGON ((longitude latitude, longitude latitude etc.))")
  
  if(is.na(worms_taxa[1])) {
    worms_taxa <- obis_taxa
  } else {
    print("Using worms_taxa")
  }
  
  #OBIS----
  obis_sp_in <- checklist(obis_taxa, 
                          geometry = regional_poly)
  
  print(paste("Searching for taxa:", obis_taxa, sep = " "))
  
  #### Pull out only the marine species list for the regional database
  # This did successfully return all marine or brackish species and no freshwater or terrestrial
  # species when used with marine = NA, brackish = NA, freshwater = F, and terrestrial = F
  # Ran tests: Brackish only success, marine only success, marine T brackish agnostic success
  if (!is.na(marine)) {
    obis_sp_in <- obis_sp_in %>% 
      filter(is_marine == marine)
  } else {
    print("marine agnostic")
  }
  if (!is.na(brackish)) {
    obis_sp_in <- obis_sp_in %>% 
      filter(is_brackish == brackish)
  } else {
    print("brackish agnostic")
  }
  if (!is.na(freshwater)) {
    obis_sp_in <- obis_sp_in %>% 
      filter(is_freshwater == freshwater)
  } else {
    print("freshwater agnostic")
  }
  if(!is.na(terrestrial)) {
    obis_sp_in <- obis_sp_in %>% 
      filter(is_terrestrial == terrestrial)
  } else {
    print("terrestrial agnostic")
  }
  
  
  obis_sp <- obis_sp_in %>% 
    filter(taxonRank == "Species"
    ) %>% 
    select(scientificName
    ) %>%
    rename(Species = scientificName
    ) %>% 
    mutate(Source = "OBIS"
    )
  
  print("First few rows of output: ")
  print(head(obis_sp))
  write.csv(obis_sp, paste(here("datasets"), sep = "", "/", obis_outputname, ".csv"), row.names = F)
  
  print("OBIS download complete.  Check your datasets folder for the output.")
  print(paste("Find your output at:", paste(here("datasets"), sep = "", "/", obis_outputname, ".csv")), sep = " ")
}

