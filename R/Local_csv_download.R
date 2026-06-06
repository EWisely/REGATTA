# Eldridge Wisely
# This copy modified for function use by Ella Crotty

#' Read user-supplied local checklist CSVs
#'
#' Reads one or more CSV files of locally-known species, drops rows
#' where `Species` is `"sp."`, unites `Genus` and `Species` into a
#' single `Genus_species` column, and writes the deduplicated result
#' to `checklist_sources/<loc_outputname>.csv` along with a `Source` column
#' marking each row as `"Local_csv"`.
#'
#' @param loc_csvs Character vector of CSV files. Each entry may be a path
#'   (absolute, or relative to your working directory) to a file anywhere on
#'   disk, **or** a bare filename (including the `.csv` extension) for a CSV
#'   in the `checklist_sources/` folder. Each CSV must have columns named
#'   exactly `Genus` and `Species` (capitalization matters).
#' @param loc_outputname Basename (no `.csv` extension) for the output
#'   file under `checklist_sources/`. Default `"Local_Species"`. Choose a
#'   distinctive name -- you will pass it into
#'   [build_regional_checklist()] later.
#'
#' @details
#' Each CSV is read from the path you give if that file exists, otherwise
#' from `checklist_sources/<name>` (alongside the outputs of
#' [GBIF_download()] and [OBIS_download()]) -- so you can point at a
#' checklist kept in another project without copying it first. Rows with
#' `Species == "sp."` are dropped because they cannot be matched against
#' species-level NCBI taxonomy downstream.
#'
#' @return Invisibly NULL. Writes
#'   `checklist_sources/<loc_outputname>.csv` as a side effect.
#'
#' @examples
#' \dontrun{
#' # Whatever local checklists you have in checklist_sources/ that carry Genus
#' # and Species columns -- fish, mammals, plants, microbes, whatever.
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
  message("Pass either a path to each CSV (absolute or relative to your ",
          "working directory), or a bare filename for a CSV kept in the ",
          "checklist_sources/ folder of your working directory.")

  Local_species <- data.frame()

  for (i in loc_csvs) {
    # Use the path as given if it points at an existing file (absolute or
    # relative to the working directory); otherwise treat it as a bare
    # filename in checklist_sources/.
    path <- if (file.exists(i)) i else here::here("checklist_sources", i)
    if (!file.exists(path))
      stop("Cannot find local CSV '", i, "'. Give a path to an existing ",
           "file, or place the file in checklist_sources/ and pass its name.")
    message("Processing ", path)

    db <- utils::read.csv(path, fileEncoding = "latin1", check.names = FALSE)

    if (!"Genus" %in% colnames(db))
      stop('"Genus" column must exist and be named exactly that')
    if (!"Species" %in% colnames(db))
      stop('"Species" column must exist and be named exactly that')

    # Combine Genus + Species into "Genus species", dropping placeholder "sp."
    keep <- db[db$Species != "sp.", , drop = FALSE]
    sp <- data.frame(
      Species = trimws(paste(keep$Genus, keep$Species)),
      Source  = "Local_csv",
      stringsAsFactors = FALSE
    )
    Local_species <- rbind(Local_species, sp)
    message(i, " processed")
  }

  Local_species <- unique(Local_species)
  dir.create(here::here("checklist_sources"), recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(here::here("checklist_sources"), paste0(loc_outputname, ".csv"))
  utils::write.csv(Local_species, out_path, row.names = FALSE)
  message("Local species CSV written to ", out_path)
  invisible(NULL)
}

