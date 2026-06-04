# Eldridge Wisely
# This copy modified for function use by Ella Crotty

#' Download a regional species checklist from OBIS
#' 
#' Downloads a list of all organisms in the OBIS (Ocean Biodiversity Information System) 
#' database in the specified taxa within the specified polygon. 
#' on a world map, and returns a list of the species it finds that is compatible with 
#' Build_regional_checklist(). This is not useful if you are not seeking ocean organisms.
#'
#' @param obis_taxa character vector. A list of taxon names, any level.
#' @param worms_taxa character vector. A list of substitute taxon names 
#' to use instead of obis_taxa, any level.
#' @param regional_poly string. Polygon describing the region of 
#' the globe to take species observations from in the OBIS database. 
#' Must be in the format POLYGON ((long lat, long lat))
#' @param obis_outputname string. Name of output file (not including filetype extension .csv).
#' @param marine logical (T or F). If T, only returns marine species. If F, only 
#' returns non-marine species. If NA, returns marine and non-marine species.
#' @param freshwater logical (T or F). If T, only returns freshwater species. If F, only 
#' returns non-freshwater species. If NA, returns freshwater and non-freshwater species.
#' @param terrestrial logical (T or F). If T, only returns terrestrial species. If F, 
#' only returns non-terrestrial species. If NA, returns terrestrial and non-terrestrial species.
#' @param brackish logical (T or F). If T, only returns brackish species. If F, only returns 
#' non-brackish species. If NA, returns brackish and non-brackish species.
#'
#' @returns Prints the filepath to the output CSV.Saves a CSV file of local species fetched from OBIS that is compatible with Build_regional_checklist.
#' 
#' @export 
#'
#' @examples
#' OBIS_download(obis_taxa = c("Salmonidae", "Copepoda"), 
#' regional_poly = "POLYGON ((-117.421875 31.952162, -91.933594 -6.315299, 
#' -81.386719 -6.315299, -76.113281 7.710992, -82.089844 8.581021, 
#' -87.011719 13.581921, -104.238281 20.303418, -112.5 32.249974, 
#' -117.421875 31.952162))", # region around the state of Washington 
#' obis_outputname = "OBIS_Species_FunctionTest", # choose this name based on 
#' # what you want your output to be. We recommend choosing something 
#' # shorter than this, because you'll have to input it into Build_regional_checklist(). 
#' # If you leave it blank, it will give the default name "OBIS_Species.csv" and 
#' # Build_regional_checklist() will also default to this name so the pipeline 
#' # will work.  This is not recommended if you are running this function 
#' multiple times for any reason. 
#' marine = T, 
#' freshwater = NA, # necessary to get salmon because it's also freshwater 
#' terrestrial = F, 
#' brackish = NA # necessary to get salmon because it's also brackish 
#' )

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

