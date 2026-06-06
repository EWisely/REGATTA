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
#'   broader** (GBIF won't accept order/family/genus/species here). List the
#'   specific classes you want -- **avoid broad ambiguous names**. For example
#'   `"Vertebrata"` is ambiguous in WoRMS (the vertebrate subphylum *and* a
#'   red-algae genus) and errors here; pass the classes instead
#'   (e.g. `c("Actinopterygii", "Chondrichthyes", "Mammalia", "Aves",
#'   "Reptilia")`). Names are validated and disambiguated by [resolve_taxa()]
#'   (kingdom-aware), so you can type any common synonym -- e.g. `Teleostei`,
#'   `Actinopteri`, `Actinopterygii`, or `Osteichthyes` all work for
#'   ray-finned fish.
#' @param worms_taxa A character vector of substitute taxon names to use
#'   instead of `obis_taxa` when looking up WoRMS IDs. NA reuses
#'   `obis_taxa`. Rarely needed now that names route through [resolve_taxa()].
#' @param regional_poly A WKT POLYGON string of the form
#'   `"POLYGON ((long lat, long lat, ...))"`. Draw a region on
#'   [wktmap.com](https://wktmap.com) and copy the generated polygon.
#' @param gbif_outputname Basename (no `.csv` extension) for the output
#'   file under `datasets/`. Default `"GBIF_Species"`.
#' @param kingdom Kingdom used by [resolve_taxa()] to disambiguate the query
#'   taxa. Default `"Animalia"`; `NULL` disables the filter.
#' @param gbif_descend_to WoRMS rank to walk each query taxon down to before
#'   matching the GBIF backbone. Default `"order"`. This is what makes
#'   ray-finned fish work: GBIF's backbone has no usable class node for them
#'   (bony fish hang under phylum Chordata with no class), so a class-level
#'   match drops every fish; descending to **order** hits a rank GBIF
#'   populates. The resulting keys are unioned with each taxon's direct
#'   backbone key, so sharks/mammals/birds keep working too. Order is the
#'   default because it is **fast** (tens of keys, not hundreds -- GBIF limits
#'   concurrent downloads and large key sets are slow to prepare) and reaches
#'   ~97% of fish families: WoRMS splits the old `Perciformes` into modern
#'   orders GBIF lacks, but GBIF still files those fish under `Perciformes`,
#'   which IS matched. The remaining ~3% are picked up by `gbif_fill_families`.
#' @param gbif_fill_families If TRUE (default) and `gbif_descend_to` is above
#'   family, additionally add GBIF keys for the *individual families* whose
#'   GBIF parent at the primary rank is not in the matched set (e.g. families
#'   GBIF files with no order), lifting coverage from ~97% to ~98%. Because the
#'   whole search is a single `occ_download` (all keys in one request -- the key
#'   count does not cost extra downloads against GBIF's concurrent-download
#'   limit), completeness is the sensible default. The only cost is a second
#'   WoRMS family walk during preparation. A phylum guard prevents stray
#'   family-name collisions from adding off-target taxa. Set FALSE to skip the
#'   family walk for a faster, order-only (~97%) preparation.
#'
#' @details
#' Takes ~15 minutes to run end to end. The function prints status
#' ("PREPARING" / "RUNNING") in the console while polling the GBIF
#' download API; this is normal. GBIF coverage is uneven -- see the
#' Common troubleshooting section of the README for known issues
#' (Osteichthyes not recognized, occasional timeouts, etc.).
#'
#' @return Invisibly NULL. Writes
#'   `datasets/<gbif_outputname>.csv` as a side effect.
#'
#' @examples
#' \dontrun{
#' # Pick the gbif_outputname for your region + group; the convention
#' # GBIF_<region>_<group> is a suggestion only.
#' # Pacific NW (OCNMS) ray-finned fish + crustaceans:
#' GBIF_download(
#'   obis_taxa       = c("Osteichthyes", "Multicrustacea"),
#'   worms_taxa      = c("Actinopterygii", "Multicrustacea"),
#'   regional_poly   = "POLYGON ((-129 52, -122 52, -122 41, -129 41, -129 52))",
#'   gbif_outputname = "GBIF_OCNMS_fish_crust"
#' )
#'
#' # The same function for any other group + region, e.g. mammals over
#' # a Sonoran Desert polygon, freshwater insects over a German
#' # river-basin polygon, soil microbiota over an Antarctic polygon --
#' # just supply the class-or-broader taxon name(s) and the WKT.
#' }
#'
#' @export
GBIF_download <- function(obis_taxa,
                          worms_taxa = NA,
                          regional_poly,
                          gbif_outputname = "GBIF_Species",
                          kingdom = "Animalia",
                          gbif_descend_to = "order",
                          gbif_fill_families = TRUE
    ) {
  query_taxa <- if (is.na(worms_taxa[1])) obis_taxa else worms_taxa

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
  
  # Validate + disambiguate query taxa to unambiguous WoRMS AphiaIDs. Replaces
  # the old per-name wm_name2id() call, which errored (206 Partial Content) on
  # ambiguous names like "Vertebrata". check_gbif = TRUE also gives us each
  # taxon's direct GBIF backbone key (reliable for Mammalia, Aves, ...).
  resolved <- resolve_taxa(query_taxa, kingdom = kingdom, check_gbif = TRUE)
  message("GBIF query resolved to: ",
          paste0(resolved$valid_name, " (AphiaID ", resolved$aphia_id, ")",
                 collapse = "; "))

  # GBIF backbone keys are assembled from up to three sources, kept small and
  # fast while still complete:
  #  (1) DIRECT: each resolved taxon's own backbone key (resolve_taxa's
  #      gbif_key) -- reliable for taxa GBIF has as nodes (Mammalia, Aves,
  #      Elasmobranchii).
  #  (2) PRIMARY descent (default `gbif_descend_to = "order"`): keys for the
  #      taxon's descendants at that rank. This is the fix for ray-finned fish
  #      -- GBIF's backbone skips the class rank for them (bony fish hang under
  #      phylum Chordata with NO class node), so a class-level match drops every
  #      fish; descending to ORDER hits a rank GBIF populates. Order is fast
  #      (tens of keys) and reaches ~97% of fish families, because GBIF files
  #      most modern-order fish under its broad `Perciformes`, which is matched.
  #  (3) GAP-FILL families (`gbif_fill_families = TRUE`): for the ~3% of fish
  #      families whose GBIF parent order is NOT in the matched primary set
  #      (GBIF gives them no order, etc.), add just those family keys. This
  #      lifts coverage to ~98% while adding only a handful of keys -- far
  #      cheaper than descending everything to family (~10x the keys, a much
  #      slower download, and GBIF throttles concurrent downloads).
  direct_keys <- resolved$gbif_key[!is.na(resolved$gbif_key)]

  # Only taxa WITHOUT a usable direct GBIF key need descending. A taxon GBIF
  # has as a node -- sharks/rays under class Elasmobranchii (key 121), Mammalia
  # (359), Aves (212) -- is fully covered by that ONE key, so we skip the
  # redundant order/family walk for it. Bony fish have no usable class node, so
  # they fall through to the order descent below.
  descend_ids <- resolved$aphia_id[is.na(resolved$gbif_key)]
  if (length(descend_ids) < nrow(resolved)) {
    message(nrow(resolved) - length(descend_ids), " taxon/taxa covered by a ",
            "single direct GBIF key (no descent needed): ",
            paste(resolved$valid_name[!is.na(resolved$gbif_key)], collapse = ", "))
  }

  # Walk the taxa that need descending down to a rank; return unique names.
  .walk <- function(rank) {
    acc <- c()
    for (id in descend_ids) {
      d <- tryCatch(
        taxize::worms_downstream(id = id, downto = rank),
        error = function(e) NULL)   # tolerate taxize rank-walk failures per taxon
      if (!is.null(d) && nrow(d) > 0) acc <- rbind(acc, d)
    }
    if (is.null(acc)) character(0) else unique(acc$name)
  }

  # (2) Primary descent.
  print(paste0("Getting WoRMS downstream ", gbif_descend_to, "s"))
  primary_names <- .walk(gbif_descend_to)
  classnameno   <- length(primary_names)
  primary_keys  <- integer(0)
  valid_phyla   <- character(0)
  if (length(primary_names) > 0) {
    bbp          <- rgbif::name_backbone_checklist(primary_names)
    keep         <- bbp$matchType != "NONE"
    primary_keys <- bbp$usageKey[keep]
    if ("phylum" %in% names(bbp)) valid_phyla <- unique(stats::na.omit(bbp$phylum[keep]))
  }

  # (3) Gap-fill families not reached by the primary nodes (skip if already
  # descending to family). Guarded by phylum so a stray family-name collision
  # (e.g. a WoRMS fish-family name that matches a GBIF mollusc family) can't
  # drag in off-target taxa.
  fill_keys <- integer(0)
  if (isTRUE(gbif_fill_families) && gbif_descend_to != "family") {
    fam_names <- .walk("family")
    if (length(fam_names) > 0) {
      bbf        <- rgbif::name_backbone_checklist(fam_names)
      bbf        <- bbf[bbf$matchType != "NONE", , drop = FALSE]
      parent_col <- paste0(gbif_descend_to, "Key")
      covered    <- if (parent_col %in% names(bbf)) bbf[[parent_col]] %in% primary_keys else FALSE
      covered[is.na(covered)] <- FALSE
      in_phylum  <- if (length(valid_phyla) && "phylum" %in% names(bbf))
                      bbf$phylum %in% valid_phyla else TRUE
      fill_keys  <- bbf$usageKey[!covered & in_phylum]
      if (length(fill_keys) > 0) {
        message("Gap-fill: added ", length(fill_keys), " GBIF family key(s) for taxa ",
                "not under a matched ", gbif_descend_to,
                " (e.g. families GBIF gives no order). Set gbif_fill_families = ",
                "FALSE to skip, or gbif_descend_to = \"family\" for full family descent.")
      }
    }
  }

  backbone_keys <- unique(c(direct_keys, primary_keys, fill_keys))
  print("GBIF backbone keys:")
  print(backbone_keys)
  backbonekeyno <- length(backbone_keys)
  if (backbonekeyno == 0) {
    stop("No usable GBIF backbone keys for your query taxa (",
         paste(query_taxa, collapse = ", "), "). GBIF's backbone lacks keys ",
         "for some broad classes (notably ray-finned fish). Use OBIS for these ",
         "taxa, or pass narrower class names GBIF recognizes.")
  }

  print("Downloading Database from GBIF")
  #### Download Checklist from GBIF ####
  download <- rgbif::occ_download(
    rgbif::pred_within(regional_poly),
    rgbif::pred_in("taxonKey", backbone_keys), # important to use pred_in
    rgbif::pred("hasCoordinate", TRUE),
    rgbif::pred("hasGeospatialIssue", FALSE),
    format = "SPECIES_LIST"
  )

  # Get download ID
  download_output <- utils::capture.output(download)
  print(download_output)
  download_id <- substr(download_output[15], 17, 41)
  print(download_id)
  
  # Source - https://stackoverflow.com/a/55851721
  # Posted by sckott
  # Retrieved 2026-02-03, License - CC BY-SA 4.0
  # Theoretically pauses the code here until the download is done but it is not working
  still_running <- TRUE
  status_ping <- 9
  print(rgbif::occ_download_meta(key = download_id)) # print info about download

  while (still_running) { # starts true
    meta <- rgbif::occ_download_meta(key = download_id) # get metadata about download
    print(meta$status) # print just the status, the rest of it won't have changed since previous print 
    status <- meta$status 
    still_running <- !(status %in% c("SUCCEEDED", "KILLED")) # SUCCEEDED or KILLED are end states so they stop the loop
    Sys.sleep(status_ping) # sleep between pings
  }
  
  ##### Import and get the species list ####
  
  GBIF_list <- rgbif::occ_download_get(download_id) %>%
    rgbif::occ_download_import()
  
  print("Download Complete")
  
  # Reprint key info because the download status messages clogged the console
  print("Taxa walked downstream:")
  print(primary_names)
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
  
  utils::write.csv(GBIF_species, paste(here::here("datasets"), sep = "", "/", gbif_outputname, ".csv"), row.names = FALSE)

  print("GBIF Download & Export Complete. Check your datasets folder for the output.")
  print(paste("Output at:", sep = " ", paste(here::here("datasets"), sep = "", "/", gbif_outputname, ".csv")))
  print(paste(backbonekeyno, sep = " ", "backbone keys found out of", classnameno, "classes in list"))
}

