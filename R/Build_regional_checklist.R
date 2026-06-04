# Eldridge Wisely
# This copy modified for function use by Ella Crotty

dataset_combine <- function(comb_inputnames = c("GBIF_Species, Local_Species, OBIS_Species"), # defaults to default outputs
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
