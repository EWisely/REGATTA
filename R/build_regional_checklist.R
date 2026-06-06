# Eldridge Wisely
# Originally `dataset_combine()`, modified for function use by Ella Crotty.
# Renamed to `build_regional_checklist()` to make the actual job legible:
# merge the GBIF, OBIS, and Local CSV outputs into one regional list.

#' Build a regional species checklist from GBIF, OBIS, and Local sources
#'
#' Reads the per-source CSVs written by [GBIF_download()],
#' [OBIS_download()], and [Local_csv_download()], stacks them,
#' deduplicates, drops any literal `"NA"` species rows, and writes the
#' result as a tab-delimited `.txt` file under `custom_db/`. The
#' output is the input to [taxonomize_checklist()].
#'
#' Run this function **once per taxonomic group** (fish, crustaceans,
#' marine mammals, ...). Don't mix groups in one call -- the
#' off-target check downstream is per-group, and a fish-plus-crustacean
#' megalist defeats the filter.
#'
#' @param comb_inputnames A character vector of basenames (no `.csv`
#'   extension) of the per-source CSVs in `datasets/`. Each entry is
#'   optional: you can supply just one source or all three. Default
#'   matches the default output names of [GBIF_download()],
#'   [Local_csv_download()], and [OBIS_download()].
#' @param comb_outputname Basename (no `.txt` extension) of the
#'   combined output file under `custom_db/`. Default
#'   `"Comprehensive_species_list"`. Choose something descriptive --
#'   you will pass this path into [taxonomize_checklist()] next.
#'
#' @return Invisibly NULL. Writes
#'   `custom_db/<comb_outputname>.txt` as a side effect. The file is
#'   tab-delimited with one species per row plus a `Source` column
#'   recording which input each species came from.
#'
#' @examples
#' \dontrun{
#' # The names below are the basenames of files in datasets/ that
#' # GBIF_download(), OBIS_download(), and Local_csv_download() wrote.
#' # comb_outputname is yours to pick; suggested convention is
#' # Comprehensive_<region>_<group>_list but it's not required.
#' build_regional_checklist(
#'   comb_inputnames = c("GBIF_my_region_my_group",
#'                       "OBIS_my_region_my_group",
#'                       "Local_my_region_my_group"),
#'   comb_outputname = "Comprehensive_my_region_my_group_list"
#' )
#' }
#'
#' @importFrom utils read.csv
#' @export
build_regional_checklist <- function(comb_inputnames = c("GBIF_Species",
                                                         "Local_Species",
                                                         "OBIS_Species"),
                                     comb_outputname = "Comprehensive_species_list") {
  if (!requireNamespace("here",  quietly = TRUE)) stop("`here` is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("`dplyr` is required.")
  if (!requireNamespace("readr", quietly = TRUE)) stop("`readr` is required.")

  Species_local_comb <- data.frame()

  for (i in comb_inputnames) {
    db <- utils::read.csv(paste(here::here("datasets", i), sep = "", ".csv"),
                          fileEncoding = "latin1",
                          check.names = FALSE)
    Species_local_comb <- unique(rbind(Species_local_comb, db))
  }

  Species_local_comb <- dplyr::filter(Species_local_comb, Species != "NA")

  message(nrow(Species_local_comb), " unique species found")

  dir.create(here::here("custom_db"), recursive = TRUE, showWarnings = FALSE)
  readr::write_delim(
    Species_local_comb,
    paste(here::here("custom_db"), sep = "", "/", comb_outputname, ".txt"),
    delim = "\t"
  )

  message("Wrote combined checklist to ",
          paste(here::here("custom_db"), sep = "", "/", comb_outputname, ".txt"))
  invisible(NULL)
}
