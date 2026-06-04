# Eldridge Wisely
# This copy modified for function use by Ella Crotty

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

