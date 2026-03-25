#Eldridge Wisely
# Modified by Ella Crotty
#10-21-24

#after doing a global database assignment (with obitools3) and a local database assignment (with vsearch, using a sintax-formatted database with both the lca function and the userout option), this script will reconcile the two results.

#Compare percent ID columns of local and global assignments, and pick the lowest taxonomic rank shared between any winning global assignments and the local species checklist.

validate_local_assignment <- function(Primer = "Please Define Primer", #troubleshooting trying to get the damn thing to take arbitrary strings
                            local_database = "Please Define local_database",
                            Local_advantage = T,  #Local advantage controls the behavior of the program when global_pctid = local_pctid.  If set to TRUE, then the local name becomes the preferred name when the percent ID values are equal.  If not, the global name remains the preferred name, and the script will check for the presence of the global species, genus, and family in the local checklist made by combining information from GBIF, OBIS, and the Darwin center, and the assignment will be downgraded to the lowest shared taxonomic rank between the global assignment and the local checklist.
                            userout = "Please Define userout",
                            lcafile = "Please Define lcafile",
                            obitoolsfile = "Please Define obitoolsfile",
                            asvsfile = "Please Define asvsfile"
                            # sql database from before
                            ) {

  library(prettyunits)
  library(taxonomizr)
  library(here)
  library(metabaR)
  library(tidyverse)
  
  # prepping empty cells in stats report
  
  gl_pctloc <- NA
  gl_post_loc_pref <- NA
  gl_post_gl_pref <- NA
  gl_change_asvs <- NA
  gl_gl_pctgl <- NA
  gl_loc_pref <- NA
  gl_gl_pref <- NA
  Loc_pctloc <- NA
  Loc_post_loc_pref <- NA
  Loc_post_gl_pref <- NA
  Loc_change_asvs <- NA
  Loc_gl_pctgl <- NA
  Loc_loc_pref <- NA
  Loc_gl_pref <- NA
  
  ##################### Begin script ################

  #This takes a long time and plenty of space, but is necessary for this script
  prepareDatabase('accessionTaxa.sql')
  print("SQL database prepared")
  #accessed 4-26-24

  vsearch_userout_file <- here("script_03_inputs", userout)
  
  vsearch_lca_file <- here("script_03_inputs", lcafile)
  
  obitools_results_file <- here("script_03_inputs", obitoolsfile)
  
  ASVs_file <- here("script_03_inputs", asvsfile)
  print("Input files imported")

  #Load global obitools database results and local vsearch database results ----
  
  
  ##Global EMBL database results in obitools3 output format----
  obi_result95<-readr::read_delim(obitools_results_file)
  colnames(obi_result95)
  obi_result95<-obi_result95%>%
    select(c(ID,TAXID,SCIENTIFIC_NAME,BEST_IDENTITY,COUNT))
  print("OBI results read")
  ##Load local database vsearch results after a lowest common ancestor (LCA) analysis (using the in-silico PCR'ed mito file added to lsu and 16S without insilico PCR)-----

  lca_vsearch<-readr::read_delim(vsearch_lca_file, col_names = c("ID","sintax"), delim = '\t')
  print("LCA inputs read")
  
  ###parse the sintax taxonomy column into 7 new columns.----
  
  lca_vsearch<-lca_vsearch%>%
    separate_wider_delim(sintax, delim = ',', names = c("Domain","Phylum","Class","Order","Family","Genus","Species"),too_few = "align_start")
  
  
  lca_vsearch$Domain<-str_remove_all(lca_vsearch$Domain,pattern = "d:")
  lca_vsearch$Phylum<-str_remove_all(lca_vsearch$Phylum, pattern = "p:")
  lca_vsearch$Class<-str_remove_all(lca_vsearch$Class, pattern = "c:")
  lca_vsearch$Order<-str_remove_all(lca_vsearch$Order, pattern = "o:")
  lca_vsearch$Family<-str_remove_all(lca_vsearch$Family, pattern = "f:")
  lca_vsearch$Genus<-str_remove_all(lca_vsearch$Genus, pattern = "g:")
  lca_vsearch$Species<-str_remove_all(lca_vsearch$Species, pattern = "s:")
  lca_vsearch$Species<-str_replace_all(lca_vsearch$Species, pattern = "_", replacement = " ")
  print("LCA vsearch taxonomy levelscleaned")
  #read in Local database vsearch results in "userout" format to get the percent identities from the local assignments
  
  pct_vsearch<-readr::read_delim(vsearch_userout_file, col_names =FALSE, delim = '\t')
  
  pct_vsearch<-pct_vsearch%>%
    separate(col = X2, sep = ";",into = c("ACC","sintax"))%>%
    dplyr::rename(ID=X1, vsearch_pctid=X3, vsearch_alnlen=X4, vsearch_mism=X5, vsearch_opens=X6)
  
  print("Pct vsearch calculated")
  #Import the obitools ASVs file so that the sequences of each can be used in making the MOTUS table for MetabaR later.
  
  
  query_ASVs<-read_delim(ASVs_file,col_names = FALSE)
  
  
  query_ASVs<-query_ASVs%>%
    separate(X1, into = c("ID","Count"), sep=" ")%>%
    dplyr::rename(Sequence=X2)%>%
    select(c("ID","Sequence"))
  
  print("ASV results read")
  
  #combine the LCA taxonomy with the percent_id column from the userout file (because in vsearch I specified only the top matches so the multiples all share the same pctid with each other, so we can pick any of them)
  vsearch_results<- left_join(lca_vsearch, pct_vsearch, by ="ID", multiple ="any")
  
  vsearch_results_seq<-full_join(vsearch_results, query_ASVs, by="ID")
  
  vsearch_results_all_matches<-full_join(lca_vsearch,pct_vsearch, by ="ID")
  
  #Inspect the local database results---- 
  ## Need to add code here to save these results for export or presentation----
  #filter the local database results to check how many sequences were assigned to each taxonomic level.
  
  print("Vsearch results merged, beginning local only")
  
  #*****Begin Local_only
  #*
  print("Calculating local stats")
  
  View(vsearch_results_all_matches)
  
  #total ASVs 
  Loc_tot=nrow(vsearch_results_all_matches)
  Loc_tot # pre-regatta local
  
  #Number of total ASVs assigned to a taxon
  Loc_asg<-nrow(vsearch_results_all_matches[is.na(vsearch_results_all_matches$sintax) ==FALSE, ])
  Loc_asg
  # sintax has any taxonomy that is available
  
  #Percent of all ASVs assigned to a taxon
  Loc_pctasg <- (Loc_asg/nrow(vsearch_results_all_matches))*100

  k_only<-vsearch_results_all_matches%>%
    dplyr::filter(is.na(Phylum))%>%
    dplyr::filter(Domain =="Eukaryota")
  Loc_kin <- nrow(k_only)
  
  phylum_only<-vsearch_results_all_matches%>%
    dplyr::filter(is.na(Class))%>%
    dplyr::filter(Domain =="Eukaryota")
  Loc_phy <- nrow(phylum_only)
  #41 observations LCA'ed to just Phylum
  
  cla_only<-vsearch_results_all_matches%>%
    dplyr::filter(is.na(Order))%>%
    dplyr::filter(Domain =="Eukaryota")
  Loc_cla <- nrow(cla_only)
  
  ord_only<-vsearch_results_all_matches%>%
    dplyr::filter(is.na(Family))%>%
    dplyr::filter(Domain =="Eukaryota")
  Loc_ord <- nrow(ord_only)
  
  fam_only<-vsearch_results_all_matches%>%
    dplyr::filter(is.na(Genus))%>%
    dplyr::filter(Domain =="Eukaryota")
  Loc_fam <- nrow(fam_only)
  
  gen_only<-vsearch_results_all_matches%>%
    dplyr::filter(is.na(Species))%>%
    dplyr::filter(Domain =="Eukaryota")
  Loc_gen <- nrow(gen_only)
  
  sp_only<-vsearch_results_all_matches%>%
    dplyr::filter(!is.na(Species))%>%
    dplyr::filter(Domain =="Eukaryota")
  Loc_sp <- nrow(sp_only)
  
  unique_sp<-vsearch_results_all_matches%>% # currently returning empty
    select(Species)%>%
    unique()
  nrow(unique_sp)
  Loc_nsp <- nrow(unique_sp)
  #46 including NA
  
  unique_genuses<-vsearch_results_all_matches%>%
    select(Genus)%>%
    unique()
  nrow(unique_genuses)
  #44 including NA
  Loc_ngen <- nrow(unique_genuses)
  
  unique_families<-vsearch_results_all_matches%>%
    select(Family)%>%
    unique()
  nrow(unique_families)
  #44 including NA
  Loc_nfam <- nrow(unique_families)
  
  unique_orders<-vsearch_results_all_matches%>%
    select(Order)%>%
    unique()
  nrow(unique_orders)
  #44 including NA
  Loc_nord <- nrow(unique_orders)
  
  unique_classes<-vsearch_results_all_matches%>%
    select(Class)%>%
    unique()
  nrow(unique_classes)
  #44 including NA
  Loc_ncla <- nrow(unique_classes)
  
  unique_phyla<-vsearch_results_all_matches%>%
    select(Phylum)%>%
    unique()
  nrow(unique_phyla)
  #44 including NA
  Loc_nphy <- nrow(unique_phyla)
  
  print("Local stats calculated")
  #Combine results from local database assignments with vsearch, and global database assignments with obitools3----
  
  #combine results of obitools and vsearch_global with lca
  
  #join by ID column
  
  lca_obi_combined<- full_join(vsearch_results_seq, obi_result95, by ="ID")
  
  lca_obi_combined<-lca_obi_combined%>%
    mutate(lca_name = coalesce(Species,Genus,Family,Order,Class,Phylum,Domain))
  
  print("Local and global assignment databases combined")
  #Compare global and local database assignments and pick the best assignment for each ASV----
  
  #if BEST_IDENTITY*100 (percent identity of global assignment) is greater than vsearch_pctid (percent identity of local assignment) then mutate fish_combined$preferred_name is SCIENTIFIC_NAME.
  
  #if vsearch_pctid (percent identity of local assignment) is greater than or equal to BEST_IDENTITY*100, then mutate fish_combined$preferred_name is lca_name (if Local_advantage ==TRUE)
  
  if (Local_advantage ==TRUE){
  best_ID_combined<-lca_obi_combined%>%
    mutate(local_pctid=
             if_else(is.na(vsearch_pctid),
                     0,
                     vsearch_pctid))%>%
    mutate(global_pctid=
             if_else(BEST_IDENTITY=="NA",
                     0,
                     signif(BEST_IDENTITY,3)*100))%>%
    mutate(preferred_pctid =
             if_else(local_pctid >= global_pctid, 
                     local_pctid, 
                     global_pctid))%>%
    mutate(preferred_name=
             if_else(local_pctid >= global_pctid, 
                     lca_name, 
                     SCIENTIFIC_NAME))%>%
    mutate(database=
             if_else(local_pctid >= global_pctid,
                     "local",
                     "global"))
  }else{
    best_ID_combined<-lca_obi_combined%>%
      mutate(local_pctid=
               if_else(is.na(vsearch_pctid),
                       0,
                       vsearch_pctid))%>%
      mutate(global_pctid=
               if_else(BEST_IDENTITY=="NA",
                       0,
                       signif(BEST_IDENTITY,3)*100))%>%
      mutate(preferred_pctid =
               if_else(local_pctid > global_pctid, 
                       local_pctid, 
                       global_pctid))%>%
      mutate(preferred_name=
               if_else(local_pctid > global_pctid, 
                       lca_name, 
                       SCIENTIFIC_NAME))%>%
      mutate(database=
               if_else(local_pctid > global_pctid,
                       "local",
                       "global"))
    
    
  }
  
  View(best_ID_combined)
  
  print("Local and global assignments chosen")
  #summarize changes to taxonomic classifications----
  ## Again, this needs to be made into an exportable table 
  
  
  #count how many got assigned to local vs. global and number of unassigned.
  # After REGATTA section
  #total ASVs 
  View(obi_result95)
  gl_tot=nrow(obi_result95)
  gl_tot # pre-regatta global
  
  ### Some more global stats ###
  #Number of total ASVs assigned to a taxon global
  gl_asg<-nrow(obi_result95[is.na(obi_result95$SCIENTIFIC_NAME) ==FALSE, ])
  # sintax has any taxonomy that is available
  
  #Percent of all ASVs assigned to a taxon global
  gl_pctasg <- (gl_asg/nrow(obi_result95))*100
  ### Some more global stats ###
  
  #to make sure we didn't lose any
  post_tot <- nrow(best_ID_combined)
  post_tot # after regatta
  #1438
  
  #Updated number of total ASVs assigned to a taxon after regatta, previously known as "assigned"
  post_asg<-nrow(best_ID_combined[is.na(best_ID_combined$preferred_name) ==FALSE, ])
  post_asg
  #546
  
  #Percent of all ASVs assigned to a taxon after combining vsearch and obitools
  post_pctasg <- (post_asg/nrow(best_ID_combined))*100
  #37.72233
  
  #compared to just obitools
  obi_assigned<-nrow(best_ID_combined[is.na(best_ID_combined$SCIENTIFIC_NAME) ==FALSE, ])
  obi_assigned
  #542
  
  #increase in assigned ASVs from global only
  post_change_asvs <- post_asg - obi_assigned
  # spreadsheet column: change in number of ASVs assigned after local reconciliation
  #4
  
  
  lca_vsearch_assigned<-nrow(best_ID_combined[is.na(best_ID_combined$lca_name) ==FALSE, ])
  lca_vsearch_assigned # assigned local only
  #188
  
  #After comparing the taxonomic assignments of global and local to get preferred names----
  
  #count the number of times the global database was used for the preferred assignment
  post_gl_pref<-nrow(best_ID_combined[best_ID_combined$database == 'global'& is.na(best_ID_combined$preferred_name) ==FALSE, ])
  
  post_gl_pref
  #403
  
  # percent local assignments preferred
  post_pctgl <- (post_gl_pref/nrow(best_ID_combined[is.na(best_ID_combined$SCIENTIFIC_NAME)==FALSE, ]))*100
  # percentage of all assignments (rows with any scientific name) that are local
  
  #put these in a new dataframe
  global_preferred_assignments<-best_ID_combined%>%
    filter(database=="global")
  
  #print global_preferred_assignments to a csv file to maybe use as supplement.
  write.csv(global_preferred_assignments, here(paste0(Primer,"_output/global_preferred_assignments_before_final_LCA.csv")))
  
  #percentage of all ASVs assigned to the global database
  (post_gl_pref/nrow(best_ID_combined))*100
  #28.02503%
  #percentage of taxonomically identified ASVs assigned to the global database
  (post_gl_pref/post_asg)*100
  #73.80952%
  
  #count the number of times the local database was used for the preferred assignment
  locally_assigned<-nrow(best_ID_combined[best_ID_combined$database == 'local'& is.na(best_ID_combined$preferred_name) ==FALSE, ])
  locally_assigned
  #23
  #percentage of all ASVs
  (locally_assigned/nrow(best_ID_combined))*100
  #1.599444%
  
  #percentage of identified sequences
  (locally_assigned/post_asg)*100
  #4.212454% 
  
  
  #number of occurrences where local assignment changed the existing global assignment
  post_loc_pref<-nrow(best_ID_combined[best_ID_combined$database == 'local'& is.na(best_ID_combined$SCIENTIFIC_NAME) ==FALSE, ])
  post_loc_pref # count of local assignment preferred
  #19 taxa updated from existing global assignments
  
  # percent local assignments preferred
  post_pctloc <- (post_loc_pref/nrow(best_ID_combined[is.na(best_ID_combined$SCIENTIFIC_NAME)==FALSE, ]))*100
  # percentage of all assignments (rows with any scientific name) that are local
  
  print("Calculating post-REGATTA taxon stats")
  
  k_only<-best_ID_combined%>%
    dplyr::filter(is.na(Phylum))%>%
    dplyr::filter(Domain =="Eukaryota")
  post_kin <- nrow(k_only)
  
  phylum_only<-best_ID_combined%>%
    dplyr::filter(is.na(Class))%>%
    dplyr::filter(Domain =="Eukaryota")
  post_phy <- nrow(phylum_only)
  #41 observations LCA'ed to just Phylum
  
  cla_only<-best_ID_combined%>%
    dplyr::filter(is.na(Order))%>%
    dplyr::filter(Domain =="Eukaryota")
  post_cla <- nrow(cla_only)
  
  ord_only<-best_ID_combined%>%
    dplyr::filter(is.na(Family))%>%
    dplyr::filter(Domain =="Eukaryota")
  post_ord <- nrow(ord_only)
  
  fam_only<-best_ID_combined%>%
    dplyr::filter(is.na(Genus))%>%
    dplyr::filter(Domain =="Eukaryota")
  post_fam <- nrow(fam_only)
  
  gen_only<-best_ID_combined%>%
    dplyr::filter(is.na(Species))%>%
    dplyr::filter(Domain =="Eukaryota")
  post_gen <- nrow(gen_only)
  
  sp_only<-best_ID_combined%>%
    dplyr::filter(!is.na(Species))%>%
    dplyr::filter(Domain =="Eukaryota")
  post_sp <- nrow(sp_only)
  
  unique_sp<-best_ID_combined%>% # currently returning empty
    select(Species)%>%
    unique()
  nrow(unique_sp)
  post_nsp <- nrow(unique_sp)
  #46 including NA
  
  unique_genuses<-best_ID_combined%>%
    select(Genus)%>%
    unique()
  nrow(unique_genuses)
  #44 including NA
  post_ngen <- nrow(unique_genuses)
  
  unique_families<-best_ID_combined%>%
    select(Family)%>%
    unique()
  nrow(unique_families)
  #44 including NA
  post_nfam <- nrow(unique_families)
  
  unique_orders<-best_ID_combined%>%
    select(Order)%>%
    unique()
  nrow(unique_orders)
  #44 including NA
  post_nord <- nrow(unique_orders)
  
  unique_classes<-best_ID_combined%>%
    select(Class)%>%
    unique()
  nrow(unique_classes)
  #44 including NA
  post_ncla <- nrow(unique_classes)
  
  unique_phyla<-best_ID_combined%>%
    select(Phylum)%>%
    unique()
  nrow(unique_phyla)
  #44 including NA
  post_nphy <- nrow(unique_phyla)
  
  #number of unique taxa after choosing the best assignment between global and local:
  
  summary_best_ID_combined<-best_ID_combined%>%
    summarise(n_globally_IDed_taxa= n_distinct(SCIENTIFIC_NAME),
              n_locally_IDed_taxa= n_distinct(lca_name),
              n_combined_IDed_taxa= n_distinct(preferred_name)
  
  )

  print("Post-merge statistics calculated")
  print("Finding LCAs between local and global IDs")
  
  #Resuming the analysis----
  
  ##Finding the lowest common ancestor between 1) locally vs globally ID'ed taxa  and 2) globally ID'ed taxa and the local species list from step 1.----
  
  ###LCA of global and local ID'ed taxa-----
  
  #need the taxonomizr database to be already prepared
  meta_best_ID_combined<-best_ID_combined%>%
    separate(sintax, into=c("sintax-to-genus", "vsearch_species"),sep = "s:", remove = FALSE)
  
  v_species<-gsub("_"," ",meta_best_ID_combined$vsearch_species)
  meta_best_ID_combined$vsearch_species<-v_species
  
  taxaId<-getId(v_species,'accessionTaxa.sql')
  print(taxaId)
  
  print("Taxa IDs got")
  
  meta_best_ID_combined$v_taxID<-taxaId
  print("vsearch TaxIDs got")
  #get the taxonomy for the vsearch TaxIDs (even though we have the sintax format already, this is getting it ready to do LCA in this program)
  
  local_levels<-getTaxonomy(taxaId,'accessionTaxa.sql')
  print(local_levels)
  print("Local TaxIDs got") # current error before this: Error: no such table: nodes

  global_taxaId<-meta_best_ID_combined$TAXID
  global_levels<-getTaxonomy(global_taxaId,'accessionTaxa.sql')
  global_levels <- as.data.frame(global_levels)
  print("Global TaxIDs got")
  
  global_taxa<-meta_best_ID_combined$SCIENTIFIC_NAME
  
  print("Calculating global stats")
  
  #*****Begin Global_only
  
  View(global_levels)
  
  k_only<-global_levels%>%
    dplyr::filter(is.na(phylum))%>%
    dplyr::filter(domain =="Eukaryota")
  gl_kin <- nrow(k_only)
  
  phylum_only<-global_levels%>%
    dplyr::filter(is.na(class))%>%
    dplyr::filter(domain =="Eukaryota")
  gl_phy <- nrow(phylum_only)
  #41 observations LCA'ed to just Phylum
  
  cla_only<-global_levels%>%
    dplyr::filter(is.na(order))%>%
    dplyr::filter(domain =="Eukaryota")
  gl_cla <- nrow(cla_only)
  
  ord_only<-global_levels%>%
    dplyr::filter(is.na(family))%>%
    dplyr::filter(domain =="Eukaryota")
  gl_ord <- nrow(ord_only)
  
  fam_only<-global_levels%>%
    dplyr::filter(is.na(genus))%>%
    dplyr::filter(domain =="Eukaryota")
  gl_fam <- nrow(fam_only)
  
  gen_only<-global_levels%>%
    dplyr::filter(is.na(species))%>%
    dplyr::filter(domain =="Eukaryota")
  gl_gen <- nrow(gen_only)
  
  sp_only<-global_levels%>%
    dplyr::filter(!is.na(species))%>%
    dplyr::filter(domain =="Eukaryota")
  gl_sp <- nrow(sp_only)
  
  unique_sp<-global_levels%>% # currently returning empty
    select(species)%>%
    unique()
  nrow(unique_sp)
  gl_nsp <- nrow(unique_sp)
  #46 including NA
  
  unique_genuses<-global_levels%>%
    select(genus)%>%
    unique()
  nrow(unique_genuses)
  #44 including NA
  gl_ngen <- nrow(unique_genuses)
  
  unique_families<-global_levels%>%
    select(family)%>%
    unique()
  nrow(unique_families)
  #44 including NA
  gl_nfam <- nrow(unique_families)
  
  unique_orders<-global_levels%>%
    select(order)%>%
    unique()
  nrow(unique_orders)
  #44 including NA
  gl_nord <- nrow(unique_orders)
  
  unique_classes<-global_levels%>%
    select(class)%>%
    unique()
  nrow(unique_classes)
  #44 including NA
  gl_ncla <- nrow(unique_classes)
  
  unique_phyla<-global_levels%>%
    select(phylum)%>%
    unique()
  nrow(unique_phyla)
  #44 including NA
  gl_nphy <- nrow(unique_phyla)
  
  print("Global stats calculated")
  
  
  #condenseTaxa(taxa)
  ##   superkingdom phylum     class      order family genus species
  ## 1 "Eukaryota"  "Chordata" "Mammalia" NA    NA     NA    NA
  #This function can also be fed a large number of grouped hits, e.g. BLAST hits for high throughput sequencing reads after filtering for the best hits for each read, and output a condensed taxonomy for each grouping:
  
  
  IDs<-meta_best_ID_combined$ID
  
  local_taxas<-as.data.frame(local_levels, row.names = paste(IDs,"_local"))
  global_taxas<-as.data.frame(global_levels, row.names = paste(IDs, "_global"))
  taxas<-rbind(local_taxas, global_taxas)
  
  doubleIDs<-append(IDs,IDs)
  condensed_lca_global_vs_local<-condenseTaxa(taxas, groupings = doubleIDs)
  
  condensed_lca_global_vs_local<-rownames_to_column(as.data.frame(condensed_lca_global_vs_local), var = "ID")
  
  condensed_lca_global_vs_local<-condensed_lca_global_vs_local%>%
    mutate(global_v_local_lca_name = coalesce(species,genus,family,order,class,phylum)) # error: no superkingdom found, deleted
  # Given a set of vectors, coalesce() finds the first non-missing value at each position.
  
  #add global_v_local_lca_name column to best_ID_combined
  
  lca_global_vs_local_assignments<-condensed_lca_global_vs_local%>%
    select(ID,global_v_local_lca_name)
  
  print("best_ID_combined")
  print(best_ID_combined[26,])
  print("lca_global_vs_local_assignments")
  print(lca_global_vs_local_assignments[26,])
  
  best_ID_combined <- full_join(best_ID_combined, lca_global_vs_local_assignments, by ="ID")
  
  print("Joined best ID and LCA dataframes")
  #after visually looking at best_ID_combined for local database sequences that were assigned to 
  # very different taxa by global and local databases, and even within the local database assignments 
  # (phylum only dataframe), there were over 3,000 sequences that matched to Calanus sinicus with a 
  # few HQ619236 HQ619232 HQ619230 HQ619237 HQ619234 HQ619232 HQ619231 HQ619235 HQ619228 from the same study.  
  # I'll check these sequences in the database file to see what's going on.  I also need to remove nans 
  # from the database file. #3214 Calanus sinicus sequences, #41 ASVs reduced ID to eukaryota after global_v_local 
  # lca step.  Also remove d:Bacteria, and check MT872704, MT672041
  
  #re-ran this code with the new cleaned vsearch local database file results and nothing was phylum only!  
  # the best_ID_combined file looks like they're converging on genus level agreements now too!
  
  
  ##When global is the database instead of local, use the LCA between global and local as preferred name----
  if (Local_advantage ==TRUE){
  best_ID_combined<-best_ID_combined%>%
    mutate(preferred_name=
             if_else(global_pctid > local_pctid & is.na(lca_name)==FALSE,
                     global_v_local_lca_name, 
                     preferred_name))%>%
    mutate(database=
             if_else(global_pctid > local_pctid & is.na(lca_name)==FALSE,
                     "lca_global_v_local",
                     database))
  }else{
    best_ID_combined<-best_ID_combined%>%
      mutate(preferred_name=
               if_else(global_pctid >= local_pctid & is.na(lca_name)==FALSE,
                       global_v_local_lca_name, 
                       preferred_name))%>%
      mutate(database=
               if_else(global_pctid >= local_pctid & is.na(lca_name)==FALSE,
                       "lca_global_v_local",
                       database))
  }
  
  #put the global names and lca reassignments in a new dataframe
  global_preferred_assignments<-best_ID_combined%>%
    filter(database%in%c("lca_global_v_local","global"))
  
  #print global_preferred_assignments to a csv file to maybe use as supplement.
  write.csv(global_preferred_assignments, here(paste0(Primer,"_output/global_preferred_assignments_after_obi-vsearch_LCA.csv")))
  
  
  #updated database resulted in #440 taxa updated from existing global assignments (down from 480), 44.11474% locally assigned sequences (out of assigned sequences), 3.94% of all ASVs, 446 ASVs had greater or equal percent ID as the global database and were therefore preferred (down from 3738 (although 3214 of those were Calanus sinicus problematic sequences, likely with lots of Ns)), 55% of assigned ASVs are still from the global database, 4.9% of all ASVs.
  #6 ASVs were assigned locally that had no assignment with the global obitools database.
  
  
  #This can happen with a fairly incomplete local database.  Next steps:
  
  
  #Compare the global taxonomic assignments and the local species checklist----
  
  
  ##When global is still the database after comparing pctid, and global_v_local LCA name is NA, find the genus or family in the local checklist and reduce the ID to that level.----
  
  
  #read in the comprehensive local checklist and isolate the genus column
  
  print(here())
  
  local_checklist<-read.csv(here("custom_db", local_database), header = FALSE)
    
    #clean it up
  local_checklist<-local_checklist%>%
      dplyr::mutate(V1=gsub("Gen. ", "",V1))%>%
      dplyr::mutate(V1=gsub("indet. ", "",V1))%>%
      dplyr::mutate(V1=gsub("\"", "",V1))%>%
      dplyr::mutate(V1=gsub("sp. ", "",V1))%>%
      dplyr::mutate(V1=gsub("cf. ", "",V1))%>%
      dplyr::rename(Scientific_name=V1)
  
  print("Local checklist cleaned")
  #separate genus from species
  local_genuses<-local_checklist%>%
    separate(Scientific_name, into = c("listed_genus", "listed_species"), sep=" ",remove = FALSE)
  
  #print("local_genuses")
  #View(local_genuses)
  #print("local_checklist")
  #View(local_checklist)
  local_checklist<-unique(full_join(local_checklist, local_genuses, 
                                    by="Scientific_name", relationship = "many-to-many"))
  ##get taxonomy for everything in the comprehensive Galapagos species checklist from taxonomizr-----
  ## "Alepisaurus ferox" gets two rows for some reason but many-to-many + unique fixed that
  #View(local_checklist)
  
  #need the taxonomizr database to be already prepared
  checklist_taxIDs<-getId(local_checklist$Scientific_name,'accessionTaxa.sql')
  #print(checklist_taxIDs)
  
  local_checklist$taxID<-checklist_taxIDs
  
  #get the taxonomy for the checklist TaxIDs
  
  checklist_levels<-getTaxonomy(checklist_taxIDs,'accessionTaxa.sql')
  #print(checklist_levels)
  checklist_levels<-as.data.frame(checklist_levels)%>%
    dplyr::rename(Scientific_name=species)
    
  
  local_checklist<-full_join(local_checklist, checklist_levels, by="Scientific_name")
  
  # If the database=="global", check if the preferred_name is in the local_checklist, 
  # if it is, then leave it.  If not, get the global_taxas entry for that ID and check 
  # if the genus is in the local_checklist$genus column, go up the columns until a match.
  
  
  #put the global names and lca reassignments in a new dataframe
  global_preferred<-best_ID_combined%>%
    filter(database=="global")
  
  checklist_best_ID_combined<-global_preferred%>%
    select(ID,Sequence,TAXID,SCIENTIFIC_NAME,global_pctid,preferred_name,preferred_pctid)
  
  global_preferred_levels<-getTaxonomy(global_preferred$TAXID,'accessionTaxa.sql')
  global_preferred_IDs<-global_preferred$ID
  
  global_preferred_levels<-as.data.frame(global_preferred_levels)%>%
    mutate(ID=global_preferred_IDs)
  
  global_preferred_levels<-full_join(global_preferred_levels,checklist_best_ID_combined)
  
  global_preferred_levels$local_relative<-NA
  #stop at family and everything beyond that level becomes NA
  global_preferred_levels<-global_preferred_levels%>%
    mutate(local_relative=
             if_else(species%in%local_checklist$Scientific_name & is.na(species)==FALSE,
                     species, 
                     if_else(genus%in%local_checklist$listed_genus& is.na(genus)==FALSE,
                             genus,
                             if_else(genus%in%local_checklist$genus& is.na(genus)==FALSE,
                                     genus,
                                     if_else(family%in%local_checklist$family& is.na(family)==FALSE,
                                             family,
                                             NA)))))
  
  
  local_relative_df<-global_preferred_levels%>%
    select(ID,local_relative)
  
  
  #put the local_relative into the best_ID_combined dataframe
  best_ID_combined<-full_join(best_ID_combined,local_relative_df, by="ID")
  best_ID_combined<-best_ID_combined%>%
    mutate(preferred_name=
             if_else(database=="global",
                     local_relative,
                     preferred_name))
  
  
  off_target_global_preferred_list<-global_preferred_levels
  
  off_target_global_preferred_list<-off_target_global_preferred_list%>%
    mutate(preferred_name=local_relative)%>%
    mutate(database="global")
  
  #print global_preferred_assignments to a csv file to maybe use as supplement.
  write.csv(off_target_global_preferred_list, here(paste0(Primer,"_output/global_preferred_assignments_after_local_db_and_checklist_LCA.csv")))
  
  
  #re-summarize changes after LCA of global vs. local.  Now all should be local.-----
  
  #count how many got assigned to local vs. global and number of unassigned.
  
  # post-checklist and LCA 
  #total ASVs
  total_ASVs=nrow(obi_result95)
  print(paste("Total ASVS =", total_ASVs))
  #to make sure we didn't lose any
  final_ASVs=nrow(best_ID_combined)
  print(paste("Final ASVS =", final_ASVs))
  #1438
  
  #Updated number of total ASVs assigned to a taxon
  assigned<-nrow(best_ID_combined[is.na(best_ID_combined$preferred_name) ==FALSE, ])
  print(paste("Number of total ASVs assigned to a taxon =", assigned))
  #545
  
  #Percent of all ASVs assigned to a taxon after combining vsearch and obitools
  (assigned/nrow(best_ID_combined))*100
  print(paste("Percent of all ASVs assigned to a taxon after combining vsearch and obitools =", (assigned/nrow(best_ID_combined))*100))
  #37.89986
  
  #compared to just obitools
  obi_assigned<-nrow(best_ID_combined[is.na(best_ID_combined$SCIENTIFIC_NAME) ==FALSE, ])
  print(paste("OBItools assigned =", obi_assigned))
  #542
  
  #this is a decrease in assignment of 526 crustacean ASVs (because of off-target amplification of insects and bryozoans, and cnidarians)
  assigned - obi_assigned
  print(paste("Number of total ASVs assigned to a taxon minus OBItools assigned =", assigned - obi_assigned))
  
  #3
  
  lca_vsearch_assigned<-nrow(best_ID_combined[is.na(best_ID_combined$lca_name) ==FALSE, ])
  print(paste("LCA VSearch assigned =", lca_vsearch_assigned))
  #188
  
  #After comparing the taxonomic assignments of global and local to get preferred names----
  
  #count the number of times the global database was used for the preferred assignment
  lca_global_v_local_assigned<-nrow(best_ID_combined[best_ID_combined$database == 'lca_global_v_local'& is.na(best_ID_combined$preferred_name) ==FALSE, ])
  
  print(paste("WHAT IS THIS =", lca_global_v_local_assigned))
  #45
  
  globally_assigned<-nrow(best_ID_combined[is.na(best_ID_combined$preferred_name) ==FALSE& best_ID_combined$preferred_name==best_ID_combined$SCIENTIFIC_NAME, ])
  print(paste("Number of times the global database was used for the preferred assignment =", globally_assigned))
  #357
  
  
  #count the number of times the local database was used for the preferred assignment
  locally_assigned<-nrow(best_ID_combined[is.na(best_ID_combined$preferred_name) ==FALSE& best_ID_combined$preferred_name==best_ID_combined$lca_name, ])
  print(paste("Number of times the local database was used for the preferred assignment =", locally_assigned))
  
  #143
  
  #percentage of all ASVs
  (locally_assigned/nrow(best_ID_combined))*100
  print(paste("Percentage of all ASVs locally assigned =", (locally_assigned/nrow(best_ID_combined))*100))
  #9.944367%
  
  #percentage of identified sequences
  (locally_assigned/assigned)*100
  print(paste("Percentage of all identified sequences locally assigned =", (locally_assigned/assigned)*100))
  #26.23853%
  
  
  #number of occurences where local assignment changed the existing global assignment
  global_to_local_assigment<-nrow(best_ID_combined[best_ID_combined$database == 'local'& is.na(best_ID_combined$SCIENTIFIC_NAME) ==FALSE, ])
  print(paste("Number of occurences where local assignment changed the existing global assignment =", global_to_local_assigment))
  #139 taxa updated from existing global assignments
  
  
  #number of global taxa that got updated with local checklist
  global_to_local_relative<-nrow(best_ID_combined[best_ID_combined$preferred_name == best_ID_combined$local_relative & best_ID_combined$preferred_name!= best_ID_combined$SCIENTIFIC_NAME &is.na(best_ID_combined$preferred_name)==FALSE, ])
  print(paste("Number of global taxa that got updated with local checklist =", global_to_local_relative))
  
  
  #Make new motu (taxa) table for MetabaR ----
  print("Making new MOTU table for MetabaR")
  
  ## MOTUs characteristics table
  
  final_names<-best_ID_combined$preferred_name
  
  final_taxaId<-getId(final_names,'accessionTaxa.sql')
  #print(final_taxaId)
  
  final_taxa<-getTaxonomy(final_taxaId,'accessionTaxa.sql')
  #print(final_taxa)
  class(final_taxa)
  
  #make a new dataframe with the final taxa
  
  final_taxa<-as.data.frame(final_taxa)
  IDs<-best_ID_combined$ID
  
  rownames(final_taxa)<-IDs
  
  nrow(final_taxa)
  print(paste("Final number of taxa =", nrow(final_taxa)))
  #1438
  
  #select only the columns I want for metabaR and rename them (lowercase sequence) taking out "preferred" for brevity.
  final_taxa$pct_id<-best_ID_combined$preferred_pctid
  
  final_taxa$Scientific_name<-final_names
  final_taxa$taxID<-final_taxaId
  final_taxa$ASV<-IDs
  final_taxa$sequence<-best_ID_combined$Sequence
  final_taxa$database<-best_ID_combined$database
  final_taxa$COUNT<-best_ID_combined$COUNT
  
  
  final_taxa<-final_taxa%>%
    mutate(pct_id=
             if_else(is.na(Scientific_name)==TRUE,
                     NA,
                     pct_id))%>%
    mutate(database=
             if_else(is.na(Scientific_name)==TRUE,
                     NA,
                     database))
  
  print("final_taxa assembled")
  class(final_taxa)
  write.csv(final_taxa,here(paste0("06_local_vs_global_results/",Primer,"_Menu_ready_for_MetabaR.csv")))
  
  write.csv(best_ID_combined,here(paste0("06_local_vs_global_results/",Primer,"_best_ID_combined.csv")))
  
  print(paste0("Finished cross validating taxonomic assignments for ", 
               Primer, 
               " with Local advantage set to ", 
               Local_advantage))
  
  
  
  post_summary_best_ID_combined<-best_ID_combined%>%
    summarise(n_global_v_checklist=n_distinct(local_relative),
              n_combined_IDed_taxa= n_distinct(preferred_name))
              
    
  post_summary_best_ID_combined
  globalvlocal<-best_ID_combined%>%
    filter(database=="lca_global_v_local")
  summaryglobalvlocal<-globalvlocal%>%
    summarise(n_global_v_local= n_distinct(preferred_name))
  summaryglobalvlocal
  
  finallocal<-best_ID_combined%>%
    filter(preferred_name==lca_name)
  summaryfinallocal<-finallocal%>%
    summarise(n_finallocal= n_distinct(preferred_name))
  summaryfinallocal
  
  finalglobal<-best_ID_combined%>%
    filter(preferred_name==SCIENTIFIC_NAME)
  summaryfinalglobal<-finalglobal%>%
    summarise(n_finalglobal= n_distinct(preferred_name))
  summaryfinalglobal
  
  print("Generating stats report")
  
  row_names <- c("total ASVs",
                 "assigned ASVs",
                 "Percentage of ASVs assigned",
                 "count of local assignment preferred",
                 "percent local assignments",
                 "count of global assignment preferred",
                 "percent global assignments",
                 "change in number of ASVs assigned",
                 "ID'ed to kingdom only",
                 "ID'ed to phylum only",
                 "ID'ed to class only",
                 "ID'ed to order only",
                 "ID'ed to family only",
                 "ID'ed to genus only",
                 "ID'ed to species",
                 "Number of phyla",
                 "Number of classes",
                 "Number of orders",
                 "Number of families",
                 "Number of genera",
                 "Number of species")
  
  global_stats <- c(gl_tot, gl_asg, gl_pctasg, gl_loc_pref, gl_pctloc,
                    gl_gl_pref, gl_gl_pctgl, gl_change_asvs,
                    gl_kin, gl_phy, gl_cla, gl_ord, gl_fam, gl_gen, gl_sp,
                    gl_nphy, gl_ncla, gl_nord, gl_nfam, gl_ngen, gl_nsp)
  
  local_stats <- c(Loc_tot, Loc_asg, Loc_pctasg, Loc_loc_pref, Loc_pctloc,
                   Loc_gl_pref, Loc_gl_pctgl, Loc_change_asvs,
                   Loc_kin, Loc_phy, Loc_cla, Loc_ord, Loc_fam, Loc_gen, Loc_sp,
                   Loc_nphy, Loc_ncla, Loc_nord, Loc_nfam, Loc_ngen, Loc_nsp)
  
  post_stats <- c(post_tot, post_asg, post_pctasg, post_loc_pref, post_pctloc,
                  post_gl_pref, post_pctgl, post_change_asvs,
                  post_kin, post_phy, post_cla, post_ord, post_fam, post_gen, post_sp,
                  post_nphy, post_ncla, post_nord, post_nfam, post_ngen, post_nsp)
  
  stats_summary <- data.frame(row_names, global_stats, local_stats, post_stats)
  View(stats_summary)
  write.csv(stats_summary, 
            here(paste0(Primer,"_output/stats_summary.csv")))
  
  
  "Validate_local_assignment complete"
  return(final_taxa)
}

mf <- "MiFish"
lc_db <- "comprehensive_galapagos_fish_list.txt" # for troubleshooting

# Test
best_ID_combined_test <- validate_local_assignment(Primer = mf, 
                          local_database = lc_db, 
                          Local_advantage = T,
                          userout = "userout_MiFish5-1db_Galapagos_top5_comprehensive_galapagos_results.txt",
                          lcafile = "lca_MiFish5-1db_Galapagos_top5_comprehensive_galapagos_results.txt",
                          obitoolsfile = "MiFish_Menu_95_named.tab",
                          asvsfile = "upper_MiFish_Menu_95_named_cleared_tags.tab"
                          )