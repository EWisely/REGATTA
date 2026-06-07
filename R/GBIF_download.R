# Eldridge Wisely
# This copy modified for function use by Ella Crotty

#' Download a GBIF species list inside a WKT polygon
#'
#' Resolves the given high-level taxa via WoRMS, finds GBIF backbone
#' keys for the resulting classes, and submits a GBIF occurrence
#' download restricted to the polygon. **Returns** the species data frame;
#' writes it to disk only if you supply `output_dir`.
#'
#' @section Time and credentials:
#' This submits an asynchronous GBIF `occ_download` and **waits for GBIF to
#' assemble it server-side, which typically takes several minutes** (longer for
#' large polygons or broad taxa) -- the call blocks while it polls. GBIF
#' requires a free account, and `rgbif` authenticates with your account
#' **credentials** (there is no separate "API key" for downloads). Register at
#' \url{https://www.gbif.org/user/profile}, then store three values in your
#' `~/.Renviron` (via `usethis::edit_r_environ()`, then restart R) so nothing is
#' hardcoded: `GBIF_USER`, `GBIF_PWD`, and `GBIF_EMAIL` (the email you
#' registered with).
#'
#' @param obis_taxa A character vector of taxon names at any level. They are
#'   validated, disambiguated, and resolved to GBIF backbone keys by
#'   [resolve_taxa()], which walks broad taxa down to the orders/families GBIF
#'   populates -- so subphyla GBIF has no usable key for (e.g. `"Vertebrata"`,
#'   `"Crustacea"`, or the ray-finned-fish classes) still work, as do common
#'   synonyms (`Teleostei`, `Actinopteri`, `Osteichthyes`) and the `"fish"` /
#'   `"vertebrates"` shorthands.
#' @param worms_taxa A character vector of substitute taxon names to use
#'   instead of `obis_taxa` when looking up WoRMS IDs. NA reuses
#'   `obis_taxa`. Rarely needed now that names route through [resolve_taxa()].
#' @param regional_poly A WKT POLYGON string of the form
#'   `"POLYGON ((long lat, long lat, ...))"`. Draw a region on
#'   [wktmap.com](https://wktmap.com) and copy the generated polygon.
#' @param output_dir Optional directory to write the result into. Default
#'   `NULL` writes nothing (the function just returns the data); supply a
#'   directory to also save `<output_dir>/<gbif_outputname>.csv` and
#'   `<output_dir>/<gbif_outputname>_download_info.txt` (the download key +
#'   GBIF citation).
#' @param gbif_outputname Basename (no `.csv` extension) for the output file,
#'   used only when `output_dir` is supplied. Default `"GBIF_Species"`.
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
#' @param existing_download Reuse an already-finished GBIF download instead of
#'   submitting a new one: a GBIF download **key** string (e.g.
#'   `"0012345-230101000000000"`) or an `occ_download` **object**. When set,
#'   the function skips taxon resolution and submission and just imports that
#'   download, so `obis_taxa` / `regional_poly` are not needed. Default `NULL`.
#'
#' @details
#' Takes ~15 minutes to run end to end. The function prints status
#' ("PREPARING" / "RUNNING") in the console while polling the GBIF
#' download API; this is normal. GBIF coverage is uneven -- see the
#' Common troubleshooting section of the README for known issues
#' (Osteichthyes not recognized, occasional timeouts, etc.).
#'
#' @return Invisibly, a data.frame with `Species` and `Source` columns. The
#'   GBIF download key and citation are attached as `attr(result,
#'   "gbif_download")` (a one-row data.frame: `download_key`, `doi`, `created`,
#'   `format`, `citation`) -- **keep this for your methods section and to reuse
#'   the download** (pass `download_key` back via `existing_download =`). It is
#'   also printed to the console and, when `output_dir` is supplied, written to
#'   `<gbif_outputname>_download_info.txt`.
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
GBIF_download <- function(obis_taxa = NULL,
                          worms_taxa = NA,
                          regional_poly = NULL,
                          gbif_outputname = "GBIF_Species",
                          kingdom = "Animalia",
                          gbif_descend_to = "order",
                          gbif_fill_families = TRUE,
                          output_dir = NULL,
                          existing_download = NULL
    ) {
 if (is.null(existing_download)) {
  if (is.null(obis_taxa) || is.null(regional_poly))
    stop("GBIF_download() needs `obis_taxa` and `regional_poly` to submit a ",
         "new download. To reuse a finished one, pass `existing_download` ",
         "(a GBIF download key or occ_download object) instead.")
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
  
  # Validate + disambiguate query taxa AND resolve their GBIF backbone keys in
  # one call: resolve_taxa() returns, per taxon, an unambiguous WoRMS AphiaID
  # plus `gbif_keys` -- the taxon's own backbone key when usable, otherwise the
  # keys reached by walking WoRMS down to orders/families (the fix for broad
  # nodes GBIF lacks a usable key for: ray-finned fish, Crustacea, Vertebrata).
  resolved <- resolve_taxa(query_taxa, kingdom = kingdom, check_gbif = TRUE,
                           gbif_descend_to = gbif_descend_to,
                           gbif_fill_families = gbif_fill_families)
  message("GBIF query resolved to: ",
          paste0(resolved$valid_name, " (AphiaID ", resolved$aphia_id, ")",
                 collapse = "; "))

  ntaxa         <- nrow(resolved)
  backbone_keys <- unique(unlist(resolved$gbif_keys, use.names = FALSE))
  backbone_keys <- backbone_keys[!is.na(backbone_keys)]
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
  print("GBIF backbone keys:")
  print(backbone_keys)

  # --- Capture the citation/provenance the user must keep ------------------
  # rgbif stamps the full GBIF citation onto the occ_download object; the DOI
  # and timestamp come from the final metadata. This is the reproducibility
  # record (and the reuse key), so we return it rather than only printing it.
  final_meta  <- rgbif::occ_download_meta(key = download_id)
  # The DOI is often blank on the submission object but populated in the final
  # metadata; take the first non-empty of the two. (rgbif stamps doi/citation/
  # created as attributes on the occ_download object -- see print.occ_download.)
  pick <- function(...) {
    for (v in list(...))
      if (!is.null(v) && length(v) == 1 && !is.na(v) && nzchar(v)) return(v)
    NA_character_
  }
  dl_doi      <- pick(final_meta$doi, attr(download, "doi"))
  dl_created  <- pick(as.character(final_meta$created), attr(download, "created"))
  dl_citation <- attr(download, "citation")
  if (is.null(dl_citation) || !nzchar(dl_citation))
    dl_citation <- paste0("GBIF Occurrence Download https://www.gbif.org/",
                          "occurrence/download/", download_id,
                          " accessed from R via rgbif")
  gbif_info <- data.frame(
    download_key = download_id, doi = dl_doi, created = dl_created,
    format = "SPECIES_LIST", citation = dl_citation, stringsAsFactors = FALSE)
 } else {
   message("Reusing existing GBIF occ_download (no new request submitted).")
   GBIF_list <- rgbif::occ_download_get(existing_download) %>%
     rgbif::occ_download_import()
   reuse_key <- as.character(existing_download)[1]
   rm_meta   <- tryCatch(rgbif::occ_download_meta(key = reuse_key),
                         error = function(e) NULL)
   gbif_info <- data.frame(
     download_key = reuse_key,
     doi     = if (!is.null(rm_meta$doi)) rm_meta$doi else NA_character_,
     created = if (!is.null(rm_meta$created)) as.character(rm_meta$created) else NA_character_,
     format  = "SPECIES_LIST",
     citation = if (!is.null(rm_meta$doi))
                  paste0("GBIF Occurrence Download ", rm_meta$doi,
                         " accessed from R via rgbif")
                else NA_character_,
     stringsAsFactors = FALSE)
 }
  
  
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
  
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    out_path <- file.path(output_dir, paste0(gbif_outputname, ".csv"))
    utils::write.csv(GBIF_species, out_path, row.names = FALSE)
    message("GBIF species list written to ", out_path)
    # Persist the citation/provenance alongside the data (opt-in via output_dir).
    info_path <- file.path(output_dir, paste0(gbif_outputname, "_download_info.txt"))
    writeLines(c(
      "GBIF download -- cite this in your paper.",
      "GBIF citation guidelines: https://www.gbif.org/citation-guidelines",
      "",
      paste0("Download key: ", gbif_info$download_key),
      paste0("DOI: ",          gbif_info$doi),
      paste0("Created: ",      gbif_info$created),
      "Citation:",
      gbif_info$citation), info_path)
    message("GBIF download info (key + citation) written to ", info_path)
  }

  # The backbone-key summary only exists on the fresh-download path.
  if (is.null(existing_download))
    message(backbonekeyno, " GBIF backbone key(s) resolved for ", ntaxa, " query taxon/taxa")

  # Surface the citation/key prominently -- it scrolls past in the status pings,
  # and the user needs it for their methods section and to reuse the download.
  message("")
  message("=== GBIF download info -- SAVE THIS (needed for your paper) ===")
  message("  Download key : ", gbif_info$download_key)
  message("                 reuse without re-downloading via ",
          'GBIF = "', gbif_info$download_key, '" (or existing_download=).')
  if (!is.na(gbif_info$doi)) message("  DOI          : ", gbif_info$doi)
  message("  Citation     : ", gbif_info$citation)
  message('  Also on the result: attr(<result>, "gbif_download")')
  message("==============================================================")
  message("")

  message("GBIF download complete (", nrow(GBIF_species), " species).")
  attr(GBIF_species, "gbif_download") <- gbif_info
  invisible(GBIF_species)
}

