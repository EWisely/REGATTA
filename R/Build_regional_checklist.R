# Eldridge Wisely
# This copy modified for function use by Ella Crotty

#' Combine species checklists from multiple sources
#' 
#' Takes the names of output of Local_csv_download, OBIS_download, and/or GBIF_download
#' and combines them into a text file that lists all of the species in your local area. 
#' All three initial datasets are optional in this function.
#' 
#' Each input name is optional, so doing both GBIF and OBIS is not required.
#'
#' @param comb_inputnames character vector. A list of names (without file extensions) 
#' of the CSV files output by OBIS_download(), GBIF_download(), 
#' and/or Local_csv_download()
#' @param comb_outputname string. A name (without file extension) of 
#' the TXT file listing local species that will be output by the function.
#'
#' @returns The function saves a TXT file where each row lists a Latin binomial species name and the source (OBIS, GBIF, or local_csv). It also prints the filepath to the TXT file.
#' @export
#'
#' @examples
#' Build_regional_checklist(comb_inputnames = c("GBIF_Species_FunctionTest", "Local_Species_FunctionTest", 
#' "OBIS_Species_FunctionTest"), # defaults to default outputs from GBIF, OBIS, 
#' # and local csv downloads
#' comb_outputname = "Comprehensive_species_list_FunctionTest" # whatever you want your full species
#' # list to be named. You will call this during Validate_local_assignment().
#' )

Build_regional_checklist <- function(comb_inputnames = c("GBIF_Species, Local_Species, OBIS_Species"), # defaults to default outputs
                            comb_outputname = "Comprehensive_species_list") {
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

  Species_local_comb <- data.frame()

  for (i in comb_inputnames) {
    db <- read.csv(paste(here("datasets", i),
                         sep = "", ".csv"),
                   fileEncoding = "latin1",
                   check.names = FALSE)
    Species_local_comb <- unique(rbind(Species_local_comb, db)) # idk what it does if there's an OBIS and GBIF record
  }

  Species_local_comb <- Species_local_comb %>%
    filter(Species!="NA")

  print(paste(nrow(Species_local_comb), sep = " ", "unique species found"))

  ##Write  list to a txt file to use as CRABS input----
  write_delim(Species_local_comb,
              paste(here("custom_db"),
                    sep = "",
                    "/",
                    comb_outputname,
                    ".txt"),
              delim = '\t')

  print("Your combined dataset should now be a .txt file in the 'custom_db' folder. See documentation for next steps.")
  print(paste("Output file at: ", sep = " ", paste(here("custom_db"), sep = "", "/", comb_outputname, ".txt")))
}
