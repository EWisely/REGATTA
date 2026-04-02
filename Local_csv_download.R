# Eldridge Wisely
# This copy modified for function use by Ella Crotty

Local_csv_download <- function(loc_csvs,
                               loc_outputname = "Local_Species") {
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
  
  print("WARNING: Make sure your csv of local species is in the datasets folder of your working directory before running this function!")
  
  #### Galapagos Marine Fish List 2016 ####
  x <- 1 # for naming database csv files
  Local_species <- data.frame() # for putting species into 
  Local_taxa <- data.frame()
  
  for (i in loc_csvs) {
    print(paste("Processing", i))
  
    db <- data.frame(read.csv(here("datasets", i),
                              fileEncoding = "latin1", 
                              check.names = FALSE))  
    
    if("Genus" %in% colnames(db)) { # make sure Genus and Species column names are correct
      correctg <- T
    } else {
      correctg <- F
    }
    if(correctg == F) stop('"Genus" column must exist and be named exactly that')
    
    if("Species" %in% colnames(db)) { # make sure Genus and Species column names are correct
      corrects <- T
    } else {
      corrects <- F
    }
    if(corrects == F) stop('"Species" column must exist and be named exactly that')
    
    # Combine Genus and species columns when species is not just sp. and filter everything else
    sp <- db %>%
      filter(`Species`!= "sp.") %>%
      unite(Genus_species, c("Genus", "Species"), sep =" ")%>%
      select(Genus_species)%>%
      rename(Species=Genus_species) %>% 
      mutate(Source = "Local_csv")
    
    # Get names at both taxonomic levels
    taxa1 <- sp %>% 
      select(Species) %>% 
      rename(taxa = Species)
    taxa2 <- db %>% 
      select(Genus) %>%
      rename(taxa = Genus)
    local_taxaComb <- unique(rbind(taxa1,
                                   taxa2)) %>% 
      mutate(Source = "Local_csv")
    
    Local_species <- rbind(Local_species, sp)
    Local_taxa <- rbind(Local_taxa, local_taxaComb)
    
    x <- x+1
    
    print(paste(i, "processed"))
  }
  
  print("First few lines of output:")
  print(head(Local_species))
  write.csv(Local_species, paste(here("datasets"), sep = "", "/", loc_outputname, ".csv"), row.names = F)
  print("Local species csv processing complete. Check your datasets folder for the output.")
  print(paste("Output in: ", sep = " ", paste(here("datasets"), sep = "", "/", loc_outputname, ".csv")))
}

# test (error)
#Local_csv_download(obis_taxa = c("Agnatha", "Chondrichthyes", "Osteichthyes"),
#                   loc_csvs = c("2016Aug24_Tirado-Sanchez_et_al_Galapagos_Pisces_Checklist.csv", "badcheck.csv"),
#                   loc_outputname = "Local_Species_FunctionTest1")

# test (success)
#Local_csv_download(loc_csvs = c("2016Aug24_Tirado-Sanchez_et_al_Galapagos_Pisces_Checklist.csv"),
#                   loc_outputname = "Local_Species_FunctionTest1")

# test (error)
#Local_csv_download(obis_taxa = c("Salmonidae", "Copepoda"),
#                   loc_csvs = c("OCNMS_local_sp_badtest.csv", "OCNMS_local_invert_badtest.csv"),
#                  loc_outputname = "OCNMS_Local_Species_FunctionTest1")

# test (success)
Local_csv_download(loc_csvs = c("OCNMS_local_sp_test.csv", "OCNMS_local_invert_test.csv"),
                   loc_outputname = "OCNMS_Local_Species_FunctionTest")