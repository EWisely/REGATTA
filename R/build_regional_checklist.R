# Eldridge Wisely
# Originally `dataset_combine()`, modified for function use by Ella Crotty.
# Renamed to `build_regional_checklist()` to make the actual job legible:
# merge the GBIF, OBIS, and any local CSV sources into one regional list.

#' Build a regional species checklist from GBIF, OBIS, and local sources
#'
#' Reads the per-source CSVs written by [GBIF_download()] and
#' [OBIS_download()], plus any local species CSVs you supply directly,
#' stacks them, deduplicates, drops any literal `"NA"` species rows, and
#' writes **two** tab-delimited `.txt` files under
#' `local_database_checklist/`, both named from your `region` and `taxa`:
#' \itemize{
#'   \item `comprehensive_<region>_<taxa>_list_for_making_localdb.txt` --
#'     species binomials only, for building a regional reference sequence
#'     database.
#'   \item `comprehensive_<region>_<taxa>_list_for_LCA.txt` -- the same
#'     binomials **plus** any retained genus-level entries, for
#'     [taxonomize_checklist()] -> [reconcile_checklist()].
#' }
#'
#' Genus-level entries (from local `Genus / sp.` rows) belong in the LCA
#' checklist -- they let a classifier call downgrade to a genus the region
#' actually has -- but **not** in the reference-DB list, where a bare genus
#' would pull in every species of that genus worldwide and break the
#' regional scope. That is why the two files differ by exactly those rows.
#'
#' Run this function **once per taxonomic group** (fish, crustaceans,
#' marine mammals, ...). Don't mix groups in one call -- the
#' off-target check downstream is per-group, and a fish-plus-crustacean
#' megalist defeats the filter.
#'
#' @param region Short region label for the output filenames, e.g.
#'   `"galapagos"`, `"california_current"`. Required.
#' @param taxa Short taxon-group label for the output filenames, e.g.
#'   `"fish"`, `"vertebrates"`, `"crustaceans"`. Required.
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
#'   tagged `Source = "Local_csv"`. A `Species` value of `"sp."` keeps the
#'   genus alone (no epithet); that genus-level entry goes into the `_for_LCA`
#'   file only, never the `_for_making_localdb` file. Point this straight at a
#'   checklist kept in another project -- nothing is copied.
#'
#' @return Invisibly NULL. Writes
#'   `local_database_checklist/comprehensive_<region>_<taxa>_list_for_making_localdb.txt`
#'   (species binomials only) and
#'   `local_database_checklist/comprehensive_<region>_<taxa>_list_for_LCA.txt`
#'   (binomials plus retained genus-level entries) as a side effect. Each is
#'   tab-delimited with a `Species` and a `Source` column.
#'
#' @examples
#' \dontrun{
#' # GBIF/OBIS basenames are files in checklist_sources/ that the two
#' # downloaders wrote; local_csvs are raw Genus+Species CSVs anywhere on
#' # disk (no copying). region + taxa name the two output files.
#' build_regional_checklist(
#'   region          = "galapagos",
#'   taxa            = "fish",
#'   comb_inputnames = c("GBIF_galapagos_fish", "OBIS_galapagos_fish"),
#'   local_csvs      = "~/other_project/data/my_local_checklist.csv"
#' )
#' # -> taxonomize the _for_LCA.txt file for reconcile_checklist().
#' }
#'
#' @importFrom utils read.csv
#' @export
build_regional_checklist <- function(region,
                                     taxa,
                                     comb_inputnames = c("GBIF_Species",
                                                         "OBIS_Species"),
                                     local_csvs = NULL) {
  if (!requireNamespace("here",  quietly = TRUE)) stop("`here` is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("`dplyr` is required.")
  if (!requireNamespace("readr", quietly = TRUE)) stop("`readr` is required.")
  if (missing(region) || missing(taxa) ||
      !is.character(region) || !is.character(taxa) ||
      length(region) != 1 || length(taxa) != 1 ||
      !nzchar(trimws(region)) || !nzchar(trimws(taxa)))
    stop("Provide a non-empty `region` and `taxa` -- they name the output ",
         'files, e.g. region = "galapagos", taxa = "fish".')

  stem <- paste0("comprehensive_", trimws(region), "_", trimws(taxa), "_list")

  empty        <- data.frame(Species = character(0), Source = character(0),
                             stringsAsFactors = FALSE)
  species_rows <- empty   # binomials: GBIF/OBIS + local non-"sp." rows
  genus_rows   <- empty   # genus-only: local "Genus / sp." rows

  # GBIF / OBIS per-source CSVs, read by basename from checklist_sources/.
  for (i in comb_inputnames) {
    db <- utils::read.csv(paste(here::here("checklist_sources", i), sep = "", ".csv"),
                          fileEncoding = "latin1",
                          check.names = FALSE)
    species_rows <- unique(rbind(species_rows, db))
  }

  # Local checklists, read by path (anywhere on disk; bare names fall back to
  # checklist_sources/) and normalized to the {Species, Source} shape the web
  # sources use. Binomials join species_rows; "Genus sp." rows become a bare
  # genus and join genus_rows (LCA only).
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

    is_sp <- trimws(db$Species) == "sp."
    bin   <- db[!is_sp, , drop = FALSE]
    gen   <- db[ is_sp, , drop = FALSE]
    if (nrow(bin))
      species_rows <- unique(rbind(species_rows, data.frame(
        Species = trimws(paste(bin$Genus, bin$Species)),
        Source  = "Local_csv", stringsAsFactors = FALSE)))
    if (nrow(gen))
      genus_rows <- unique(rbind(genus_rows, data.frame(
        Species = trimws(gen$Genus),
        Source  = "Local_csv", stringsAsFactors = FALSE)))
  }

  species_rows <- dplyr::filter(species_rows, Species != "NA")
  lca_rows     <- dplyr::filter(unique(rbind(species_rows, genus_rows)),
                                Species != "NA")

  dir.create(here::here("local_database_checklist"), recursive = TRUE, showWarnings = FALSE)
  db_path  <- paste0(here::here("local_database_checklist"), "/",
                     stem, "_for_making_localdb.txt")
  lca_path <- paste0(here::here("local_database_checklist"), "/",
                     stem, "_for_LCA.txt")
  readr::write_delim(species_rows, db_path,  delim = "\t")
  readr::write_delim(lca_rows,     lca_path, delim = "\t")

  message(nrow(species_rows), " species-level entries -> ", db_path,
          " (for building the local reference DB)")
  message(nrow(lca_rows), " entries (incl. ", nrow(genus_rows),
          " genus-level) -> ", lca_path, " (for reconcile_checklist)")
  invisible(NULL)
}
