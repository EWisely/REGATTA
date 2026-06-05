# Eldridge Wisely
# Originally `dataset_combine()`, modified for function use by Ella Crotty.
# Renamed to `build_regional_checklist()` to make the actual job legible:
# merge the GBIF, OBIS, and Local CSV outputs into one regional list.

#' Build a regional species checklist by combining GBIF, OBIS, and local sources
#'
#' Reads the per-source CSVs written by [GBIF_download()], [OBIS_download()],
#' and [Local_csv_download()] (one row per species, with columns
#' `Species` and `Source`), stacks them, deduplicates, drops any literal
#' `"NA"` species rows, and writes the result as a tab-delimited `.txt`
#' file under `custom_db/`. The output is the input to [taxonomize_checklist()].
#'
#' Run this function **once per taxonomic group** (fish, crustaceans, ...).
#' Don't mix groups in one call — the off-target check is per-group.
#'
#' @param comb_inputnames Character vector of basenames (no extension) of the
#'   per-source CSVs in the `datasets/` directory. Defaults to the basenames
#'   the upstream download functions write.
#' @param comb_outputname Basename (no extension) of the combined output
#'   file under `custom_db/`.
#'
#' @return Invisibly NULL; writes
#'   `custom_db/<comb_outputname>.txt` as a side effect.
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

  readr::write_delim(
    Species_local_comb,
    paste(here::here("custom_db"), sep = "", "/", comb_outputname, ".txt"),
    delim = "\t"
  )

  message("Wrote combined checklist to ",
          paste(here::here("custom_db"), sep = "", "/", comb_outputname, ".txt"))
  invisible(NULL)
}
