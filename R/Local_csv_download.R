# Eldridge Wisely
# This copy modified for function use by Ella Crotty

#' Read user-supplied local checklist CSVs
#'
#' Reads one or more CSV files of locally-known species, drops rows
#' where `Species` is `"sp."`, unites `Genus` and `Species` into a
#' single `Genus_species` column, and writes the deduplicated result
#' to `datasets/<loc_outputname>.csv` along with a `Source` column
#' marking each row as `"Local_csv"`.
#'
#' @param loc_csvs Character vector of CSV file names (**including** the
#'   `.csv` extension) inside the `datasets/` folder. Each CSV must have
#'   columns named exactly `Genus` and `Species` (capitalization matters).
#' @param loc_outputname Basename (no `.csv` extension) for the output
#'   file under `datasets/`. Default `"Local_Species"`. Choose a
#'   distinctive name — you will pass it into
#'   [build_regional_checklist()] later.
#'
#' @details
#' Input CSVs live in `datasets/` (alongside the outputs of
#' [GBIF_download()] and [OBIS_download()]). Rows with `Species ==
#' "sp."` are dropped because they cannot be matched against
#' species-level NCBI taxonomy downstream.
#'
#' @return Invisibly NULL. Writes
#'   `datasets/<loc_outputname>.csv` as a side effect.
#'
#' @examples
#' \dontrun{
#' # Whatever local checklists you have in datasets/ that carry Genus
#' # and Species columns — fish, mammals, plants, microbes, whatever.
#' # The loc_outputname is yours to pick:
#' Local_csv_download(
#'   loc_csvs       = c("regional_checklist_A.csv", "regional_checklist_B.csv"),
#'   loc_outputname = "Local_<my_region>_<my_group>"
#' )
#' }
#'
#' @export
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

