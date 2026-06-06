# Eldridge Wisely
# Originally `dataset_combine()`, modified for function use by Ella Crotty.
# Now the orchestrator for the checklist-building stage: it can run the
# OBIS/GBIF downloads (or reuse pre-made sources), fold in local CSVs, write
# the two regional lists, and taxonomize the LCA list in the background.

#' Build a regional species checklist (orchestrator)
#'
#' The basic entry point for the checklist-building stage. Pulls regional
#' species records from OBIS and/or GBIF (or reuses pre-made sources), folds
#' in any local `Genus`+`Species` CSVs, deduplicates, and writes **two**
#' tab-delimited lists under `local_database_checklist/`:
#' \itemize{
#'   \item `comprehensive_<region>_<label>_list_for_making_localdb.txt` --
#'     species binomials only, for a reference-database builder (e.g. CRABS).
#'   \item `comprehensive_<region>_<label>_list_for_LCA.txt` -- those binomials
#'     **plus** any retained genus-level entries, for the checklist LCA.
#' }
#' The LCA list is then **taxonomized in the background** (when an
#' `accessionTaxa.sql` database is available), writing
#' `comprehensive_<region>_<label>_list_for_LCA_taxonomized.rds` ready for
#' [reconcile_checklist()] / [run_regatta()]. If no database is available the
#' taxonomize step is skipped here and runs (with a warning) at the LCA step
#' instead.
#'
#' [OBIS_download()], [GBIF_download()], and [taxonomize_checklist()] remain
#' exported for direct/advanced use (caching, inspection), but the typical
#' user only calls this wrapper.
#'
#' @param region Region label for the output filenames, e.g. `"galapagos"`.
#'   Required.
#' @param label Short taxon-group label for the output filenames, e.g.
#'   `"fish"`. Kept separate from `taxa` so a long list of query names can
#'   still produce a short filename. Required.
#' @param taxa The query: taxon name(s) or shortcut(s) (e.g. `"fish"`,
#'   `"Mammalia"`) passed to the downloaders. Required only when a fresh
#'   OBIS or GBIF download will actually run.
#' @param regional_poly WKT `POLYGON((...))` string. Required only when a
#'   fresh download will run.
#' @param OBIS Source control for OBIS: `TRUE` (default) runs
#'   [OBIS_download()]; `FALSE` skips OBIS; a path to a pre-made OBIS source
#'   CSV (with `Species`/`Source` columns) feeds that instead of downloading.
#' @param GBIF Source control for GBIF: `FALSE` (default) skips GBIF; `TRUE`
#'   submits and runs a fresh `occ_download`; a GBIF download **key** string
#'   or an `occ_download` **object** reuses an already-finished download
#'   (no resubmission); a path to a pre-made GBIF source CSV feeds that.
#' @param CSV Local checklist input: `FALSE` (default) for none, or a
#'   character vector of paths to local `Genus`+`Species` CSV(s) (anywhere on
#'   disk). A `Species` value of `"sp."` keeps the genus alone (no epithet),
#'   and that genus-level entry goes into the `_for_LCA` file only.
#' @param marine,freshwater,terrestrial,brackish OBIS habitat filters, passed
#'   to [OBIS_download()]. Defaults (`marine = TRUE`, `terrestrial = FALSE`,
#'   the others `NA`) keep marine plus anadromous/catadromous species (which
#'   carry a freshwater flag) and drop only land contaminants. Set explicitly
#'   to override.
#' @param sql_path Path to the local `accessionTaxa.sql` taxonomizr DB used
#'   for the background taxonomize step. If the file is absent, taxonomizing
#'   is deferred to the LCA step.
#' @param kingdom Kingdom used by [resolve_taxa()] to disambiguate query
#'   taxa in the downloaders. Default `"Animalia"`.
#' @param gbif_fill_families Passed to [GBIF_download()] for a fresh GBIF run.
#'
#' @return Invisibly a list with `db_path`, `lca_path`, and (when
#'   taxonomized here) `rds_path`. Writes the two `.txt` lists and, when a DB
#'   is available, the taxonomized `.rds`.
#'
#' @examples
#' \dontrun{
#' # Fish in the Galapagos: OBIS (default) + a local checklist, no GBIF.
#' build_regional_checklist(
#'   region        = "galapagos",
#'   label         = "fish",
#'   taxa          = "fish",
#'   regional_poly = "POLYGON ((-92 2, -89 2, -89 -2, -92 -2, -92 2))",
#'   CSV           = "~/other_project/galapagos_fish_checklist.csv"
#' )
#'
#' # Add GBIF, reusing an occ_download you already ran:
#' build_regional_checklist(
#'   region = "galapagos", label = "fish", taxa = "fish",
#'   regional_poly = "POLYGON ((-92 2, -89 2, -89 -2, -92 -2, -92 2))",
#'   GBIF   = "0012345-230101000000000"   # an existing GBIF download key
#' )
#' }
#'
#' @importFrom utils read.csv
#' @export
build_regional_checklist <- function(region,
                                     label,
                                     taxa = NULL,
                                     regional_poly = NULL,
                                     OBIS = TRUE,
                                     GBIF = FALSE,
                                     CSV  = FALSE,
                                     marine = TRUE,
                                     freshwater = NA,
                                     terrestrial = FALSE,
                                     brackish = NA,
                                     sql_path = "accessionTaxa.sql",
                                     kingdom = "Animalia",
                                     gbif_fill_families = TRUE) {
  if (!requireNamespace("here",  quietly = TRUE)) stop("`here` is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("`dplyr` is required.")
  if (!requireNamespace("readr", quietly = TRUE)) stop("`readr` is required.")
  if (missing(region) || missing(label) ||
      !is.character(region) || !is.character(label) ||
      length(region) != 1 || length(label) != 1 ||
      !nzchar(trimws(region)) || !nzchar(trimws(label)))
    stop("Provide a non-empty `region` and `label` -- they name the output ",
         'files, e.g. region = "galapagos", label = "fish".')

  stem      <- paste0("comprehensive_", trimws(region), "_", trimws(label), "_list")
  src_dir   <- here::here("checklist_sources")
  needs_dl  <- isTRUE(OBIS) || isTRUE(GBIF)
  if (needs_dl && (is.null(taxa) || is.null(regional_poly)))
    stop("Running a fresh OBIS/GBIF download needs both `taxa` and ",
         "`regional_poly`. Set OBIS/GBIF to FALSE or to a pre-made source ",
         "to skip downloading.")

  empty        <- data.frame(Species = character(0), Source = character(0),
                             stringsAsFactors = FALSE)
  species_rows <- empty   # binomials: OBIS/GBIF + local non-"sp." rows
  genus_rows   <- empty   # genus-only: local "Genus / sp." rows

  # --- OBIS source --------------------------------------------------------
  if (isTRUE(OBIS)) {
    out <- paste0("OBIS_", trimws(region), "_", trimws(label))
    OBIS_download(obis_taxa = taxa, regional_poly = regional_poly,
                  obis_outputname = out, kingdom = kingdom,
                  marine = marine, freshwater = freshwater,
                  terrestrial = terrestrial, brackish = brackish)
    species_rows <- unique(rbind(species_rows,
      .regatta_read_source_csv(file.path(src_dir, paste0(out, ".csv")), "OBIS")))
  } else if (is.character(OBIS)) {
    species_rows <- unique(rbind(species_rows,
      .regatta_read_source_csv(OBIS, "OBIS")))
  }

  # --- GBIF source --------------------------------------------------------
  if (!.regatta_is_off(GBIF)) {
    out <- paste0("GBIF_", trimws(region), "_", trimws(label))
    gbif_csv <- file.path(src_dir, paste0(out, ".csv"))
    if (isTRUE(GBIF)) {
      GBIF_download(obis_taxa = taxa, regional_poly = regional_poly,
                    gbif_outputname = out, kingdom = kingdom,
                    gbif_fill_families = gbif_fill_families)
      species_rows <- unique(rbind(species_rows,
        .regatta_read_source_csv(gbif_csv, "GBIF")))
    } else if (is.character(GBIF) && length(GBIF) == 1 &&
               (grepl("\\.csv$", GBIF, ignore.case = TRUE) || file.exists(GBIF))) {
      species_rows <- unique(rbind(species_rows,
        .regatta_read_source_csv(GBIF, "GBIF")))            # pre-made CSV
    } else {
      # a download key (character) or an occ_download object -> reuse it
      GBIF_download(existing_download = GBIF, gbif_outputname = out)
      species_rows <- unique(rbind(species_rows,
        .regatta_read_source_csv(gbif_csv, "GBIF")))
    }
  }

  # --- Local CSV(s) -------------------------------------------------------
  if (!.regatta_is_off(CSV)) {
    for (lc in CSV) {
      path <- if (file.exists(lc)) lc else file.path(src_dir, lc)
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
  }

  if (nrow(species_rows) == 0 && nrow(genus_rows) == 0)
    stop("No species found from any source. Enable OBIS/GBIF or supply a CSV.")

  species_rows <- dplyr::filter(species_rows, Species != "NA")
  lca_rows     <- dplyr::filter(unique(rbind(species_rows, genus_rows)),
                                Species != "NA")

  out_dir  <- here::here("local_database_checklist")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  db_path  <- file.path(out_dir, paste0(stem, "_for_making_localdb.txt"))
  lca_path <- file.path(out_dir, paste0(stem, "_for_LCA.txt"))
  readr::write_delim(species_rows, db_path,  delim = "\t")
  readr::write_delim(lca_rows,     lca_path, delim = "\t")
  message(nrow(species_rows), " species-level entries -> ", db_path,
          " (for building the local reference DB)")
  message(nrow(lca_rows), " entries (incl. ", nrow(genus_rows),
          " genus-level) -> ", lca_path, " (for reconcile_checklist)")

  # --- Background taxonomize of the LCA list ------------------------------
  rds_path <- NULL
  if (file.exists(sql_path)) {
    message("Taxonomizing the LCA checklist in the background ...")
    tx <- taxonomize_checklist(lca_path, sql_path = sql_path)
    rds_path <- file.path(out_dir, paste0(stem, "_for_LCA_taxonomized.rds"))
    saveRDS(tx, rds_path)
    message("Taxonomized checklist -> ", rds_path,
            " (ready for reconcile_checklist / run_regatta)")
  } else {
    message("No taxonomizr DB at '", sql_path, "' -- skipping taxonomize here. ",
            "It will run (with a warning) at the LCA step. Pass `sql_path` to ",
            "a built accessionTaxa.sql to taxonomize now.")
  }

  invisible(list(db_path = db_path, lca_path = lca_path, rds_path = rds_path))
}

# FALSE / NULL / NA all mean "this source is off". A non-FALSE value (TRUE, a
# path, a key, an object) means "use it".
.regatta_is_off <- function(x) is.null(x) || (is.logical(x) && length(x) == 1 && !isTRUE(x))

# At the LCA step, accept a checklist that may not be taxonomized yet: if it
# lacks the 7 rank columns, warn and taxonomize it on the fly. A path, a
# character vector of names, or a names-only data.frame all work.
# build_regional_checklist() normally taxonomizes ahead of time so this is a
# no-op, but a hand-fed raw checklist still runs.
.regatta_ensure_taxonomized <- function(checklist, sql_path = "accessionTaxa.sql") {
  ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  if (is.data.frame(checklist) && all(ranks %in% names(checklist)))
    return(checklist)
  warning("checklist was not pre-taxonomized; taxonomizing it now via ",
          "taxonomize_checklist(). Build it with build_regional_checklist() ",
          "to do this ahead of time and silence this warning.", call. = FALSE)
  taxonomize_checklist(checklist, sql_path = sql_path)
}

# Read a per-source checklist CSV (OBIS/GBIF output or a fed-in equivalent)
# down to the {Species, Source} shape the combiner stacks.
.regatta_read_source_csv <- function(path, kind) {
  if (!file.exists(path))
    stop(kind, " source CSV not found: ", path)
  db <- utils::read.csv(path, fileEncoding = "latin1", check.names = FALSE)
  if (!"Species" %in% names(db))
    stop(kind, " source CSV '", path, "' needs a 'Species' column.")
  if (!"Source" %in% names(db)) db$Source <- kind
  db[, c("Species", "Source"), drop = FALSE]
}
