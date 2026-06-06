# Eldridge Wisely
# Originally `dataset_combine()`, modified for function use by Ella Crotty.
# Renamed to `build_regional_checklist()` to make the actual job legible:
# merge the GBIF, OBIS, and any local CSV sources into one regional list.

#' Build a regional species checklist from GBIF, OBIS, and local sources
#'
#' Reads the per-source CSVs written by [GBIF_download()] and
#' [OBIS_download()], plus any local species CSVs you supply directly,
#' stacks them, deduplicates, drops any literal `"NA"` species rows, and
#' writes the result as a tab-delimited `.txt` file under
#' `local_database_checklist/`. The output is the input to
#' [taxonomize_checklist()].
#'
#' Run this function **once per taxonomic group** (fish, crustaceans,
#' marine mammals, ...). Don't mix groups in one call -- the
#' off-target check downstream is per-group, and a fish-plus-crustacean
#' megalist defeats the filter.
#'
#' @param comb_inputnames A character vector of basenames (no `.csv`
#'   extension) of the per-source CSVs in `checklist_sources/` -- i.e. the
#'   outputs of [GBIF_download()] and [OBIS_download()]. Each entry is
#'   optional; supply whichever of the two web sources you have. Default
#'   matches the default output names of [GBIF_download()] and
#'   [OBIS_download()].
#' @param local_csvs Optional character vector of local checklist CSVs to
#'   fold in. Each entry may be a path (absolute, or relative to your
#'   working directory) to a file anywhere on disk, **or** a bare filename
#'   for a CSV in `checklist_sources/`. Each must have columns named exactly
#'   `Genus` and `Species`; rows are normalized to `"Genus species"` and
#'   tagged `Source = "Local_csv"`, and rows with `Species == "sp."` are
#'   dropped. Point this straight at a checklist kept in another project --
#'   nothing is copied.
#' @param comb_outputname Basename (no `.txt` extension) of the
#'   combined output file under `local_database_checklist/`. Default
#'   `"Comprehensive_species_list"`. Choose something descriptive --
#'   you will pass this path into [taxonomize_checklist()] next.
#'
#' @return Invisibly NULL. Writes
#'   `local_database_checklist/<comb_outputname>.txt` as a side effect. The
#'   file is tab-delimited with one species per row plus a `Source` column
#'   recording which input each species came from.
#'
#' @examples
#' \dontrun{
#' # GBIF/OBIS basenames are files in checklist_sources/ that the two
#' # downloaders wrote; local_csvs are raw Genus+Species CSVs anywhere on
#' # disk (no copying). comb_outputname is yours to pick.
#' build_regional_checklist(
#'   comb_inputnames = c("GBIF_my_region_my_group",
#'                       "OBIS_my_region_my_group"),
#'   local_csvs      = "~/other_project/data/my_local_checklist.csv",
#'   comb_outputname = "Comprehensive_my_region_my_group_list"
#' )
#' }
#'
#' @importFrom utils read.csv
#' @export
build_regional_checklist <- function(comb_inputnames = c("GBIF_Species",
                                                         "OBIS_Species"),
                                     local_csvs = NULL,
                                     comb_outputname = "Comprehensive_species_list") {
  if (!requireNamespace("here",  quietly = TRUE)) stop("`here` is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("`dplyr` is required.")
  if (!requireNamespace("readr", quietly = TRUE)) stop("`readr` is required.")

  Species_local_comb <- data.frame()

  # GBIF / OBIS per-source CSVs, read by basename from checklist_sources/.
  for (i in comb_inputnames) {
    db <- utils::read.csv(paste(here::here("checklist_sources", i), sep = "", ".csv"),
                          fileEncoding = "latin1",
                          check.names = FALSE)
    Species_local_comb <- unique(rbind(Species_local_comb, db))
  }

  # Local checklists, read by path (anywhere on disk; bare names fall back to
  # checklist_sources/) and normalized inline to the {Species, Source} shape
  # the web sources already use. Genus + Species -> "Genus species".
  for (lc in local_csvs) {
    path <- if (file.exists(lc)) lc else here::here("checklist_sources", lc)
    if (!file.exists(path))
      stop("Cannot find local CSV '", lc, "'. Give a path to an existing ",
           "file, or place the file in checklist_sources/ and pass its name.")
    db <- utils::read.csv(path, fileEncoding = "latin1", check.names = FALSE)
    if (!"Genus" %in% colnames(db))
      stop('Local CSV "', lc, '" needs a column named exactly "Genus".')
    if (!"Species" %in% colnames(db))
      stop('Local CSV "', lc, '" needs a column named exactly "Species".')
    keep <- db[db$Species != "sp.", , drop = FALSE]
    sp <- data.frame(
      Species = trimws(paste(keep$Genus, keep$Species)),
      Source  = "Local_csv",
      stringsAsFactors = FALSE
    )
    Species_local_comb <- unique(rbind(Species_local_comb, sp))
  }

  Species_local_comb <- dplyr::filter(Species_local_comb, Species != "NA")

  message(nrow(Species_local_comb), " unique species found")

  dir.create(here::here("local_database_checklist"), recursive = TRUE, showWarnings = FALSE)
  readr::write_delim(
    Species_local_comb,
    paste(here::here("local_database_checklist"), sep = "", "/", comb_outputname, ".txt"),
    delim = "\t"
  )

  message("Wrote combined checklist to ",
          paste(here::here("local_database_checklist"), sep = "", "/", comb_outputname, ".txt"))
  invisible(NULL)
}
