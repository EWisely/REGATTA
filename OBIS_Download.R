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
  
  print("WARNING: Check the documentation or run setup_explainer() for setup 
        steps before running this function!")
  
  if(is.na(worms_taxa[1])) {
    worms_taxa <- obis_taxa
  } else {
    print("Using worms_taxa")
  }
  
  #OBIS----
  obis_sp_in <- checklist(obis_taxa, 
                          geometry = regional_poly)
  
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
  
  write.csv(obis_sp, paste(here("datasets"), sep = "", "/", obis_outputname, ".csv"), row.names = F)
  
  print("OBIS download complete.  Check your datasets folder for the output.")
}

# Test Run

#OBIS_download(obis_taxa = c("Agnatha", "Chondrichthyes", "Osteichthyes"),
#              worms_taxa = c("Agnatha", "Chondrichthyes", "Actinopterygii"),
#              regional_poly = "POLYGON ((-117.421875 31.952162, -91.933594 -6.315299, -81.386719 -6.315299, -76.113281 7.710992, -82.089844 8.581021, -87.011719 13.581921, -104.238281 20.303418, -112.5 32.249974, -117.421875 31.952162))",
#              obis_outputname ="OBIS_Species_FunctionTest1",
#              marine = NA,
#              freshwater = F,
#              terrestrial = F,
#              brackish = NA
#)

# Test Run 2
OBIS_download(obis_taxa = c("Salmonidae", "Copepoda"),
              regional_poly = "POLYGON ((-117.421875 31.952162, -91.933594 -6.315299, -81.386719 -6.315299, -76.113281 7.710992, -82.089844 8.581021, -87.011719 13.581921, -104.238281 20.303418, -112.5 32.249974, -117.421875 31.952162))",
              obis_outputname ="OBIS_OCNMS_Species_FunctionTest",
              marine = T,
              freshwater = NA, # necessary to get salmon probably (yup it's freshwater, brackish, and marine)
              terrestrial = F,
              brackish = NA
              )
# returned a reasonable list