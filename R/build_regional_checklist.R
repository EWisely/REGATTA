# Eldridge Wisely
# Originally `dataset_combine()`, modified for function use by Ella Crotty.
# Now the orchestrator for the checklist-building stage: it can run the
# OBIS/GBIF downloads (or reuse pre-made sources), fold in local CSVs,
# taxonomize the LCA list, and RETURN everything. It writes files only when
# given an output_dir.

#' Build a regional species checklist (orchestrator)
#'
#' The basic entry point for the checklist-building stage. Pulls regional
#' species records from OBIS and/or GBIF (or reuses pre-made sources), folds
#' in any local `Genus`+`Species` CSVs, deduplicates, taxonomizes the LCA
#' list (when a `taxonomizr` database is available), and **returns** the
#' results. Three outputs are produced:
#' \itemize{
#'   \item `for_making_localdb` -- a bare, sorted character vector of the unique
#'     species binomials (no columns). Write it one-name-per-line and feed it to
#'     a reference-database builder (e.g. CRABS).
#'   \item `checklist_summary` -- the full taxonomized table, including how each
#'     name resolved (scientific name / synonym / unresolved). The audit view.
#'     `NULL` when no `taxonomizr` database was available.
#'   \item `for_LCA` -- the `taxID` column plus the 7 ranks only (resolved rows,
#'     genus-level entries retained), cleaned and ready to hand straight to
#'     [run_regatta()] / [reconcile_checklist()]. When no database is available
#'     this is the un-taxonomized name list, deferred to the reconcile step.
#' }
#' Nothing is written to disk unless you pass `output_dir`; then the outputs are
#' saved there, named from `region` and `label`.
#'
#' [OBIS_download()], [GBIF_download()], and [taxonomize_checklist()] remain
#' exported for direct/advanced use, but the typical user only calls this
#' wrapper.
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
#'   runs a fresh `occ_download`; a GBIF download **key** string or an
#'   `occ_download` **object** reuses an already-finished download; a path to
#'   a pre-made GBIF source CSV feeds that.
#' @param CSV Local checklist input: `FALSE` (default) for none, or a
#'   character vector of paths to local `Genus`+`Species` CSV(s). A `Species`
#'   value of `"sp."` keeps the genus alone (no epithet); that genus-level
#'   entry goes into `for_LCA` only, never `for_making_localdb`.
#' @param marine,freshwater,terrestrial,brackish OBIS habitat filters, passed
#'   to [OBIS_download()]. Defaults keep marine plus anadromous/catadromous
#'   species and drop only land contaminants.
#' @param sql_path Path to the local `accessionTaxa.sql` taxonomizr DB used to
#'   taxonomize the LCA list. If absent, the LCA list is returned
#'   un-taxonomized and taxonomizing is deferred to the reconcile step.
#' @param output_dir Optional directory to write outputs into. Default `NULL`
#'   writes nothing (everything is returned); supply a directory to also save
#'   `comprehensive_<region>_<label>_list_for_making_localdb.txt` (one name per
#'   line), `..._checklist_summary.csv`, and `..._for_LCA.rds` (or `.txt` when
#'   no database was available).
#' @param kingdom Kingdom used by [resolve_taxa()] to disambiguate query
#'   taxa in the downloaders. Default `"Animalia"`.
#' @param gbif_fill_families Passed to [GBIF_download()] for a fresh GBIF run.
#'
#' @return Invisibly, a list with `for_making_localdb` (a bare character vector
#'   of unique species binomials, for a reference-database builder),
#'   `checklist_summary` (the full taxonomized table with per-name resolution
#'   status, or `NULL` if no DB was available), and `for_LCA` (`taxID` + the 7
#'   ranks, ready for the LCA step). Pass `for_LCA` straight to [run_regatta()]
#'   / [reconcile_checklist()].
#'
#' @examples
#' \dontrun{
#' # OBIS (default) + a local checklist, no GBIF, no taxonomizr DB needed yet.
#' cl <- build_regional_checklist(
#'   region        = "galapagos",
#'   label         = "fish",
#'   taxa          = "fish",
#'   regional_poly = "POLYGON ((-92 2, -89 2, -89 -2, -92 -2, -92 2))",
#'   CSV           = "~/other_project/galapagos_fish_checklist.csv"
#' )
#' run_regatta(input = "MiFish_obitools.tab", checklist = cl$for_LCA,
#'             sql_path = "/path/to/accessionTaxa.sql")
#'
#' # Persist the lists to a directory of your choosing:
#' build_regional_checklist(region = "galapagos", label = "fish", taxa = "fish",
#'   regional_poly = "POLYGON ((-92 2, -89 2, -89 -2, -92 -2, -92 2))",
#'   sql_path = "/path/to/accessionTaxa.sql", output_dir = "my_checklists")
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
                                     output_dir = NULL,
                                     kingdom = "Animalia",
                                     gbif_fill_families = TRUE) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("`dplyr` is required.")
  if (!requireNamespace("readr", quietly = TRUE)) stop("`readr` is required.")
  if (missing(region) || missing(label) ||
      !is.character(region) || !is.character(label) ||
      length(region) != 1 || length(label) != 1 ||
      !nzchar(trimws(region)) || !nzchar(trimws(label)))
    stop("Provide a non-empty `region` and `label` -- they name the output ",
         'files, e.g. region = "galapagos", label = "fish".')

  needs_dl <- isTRUE(OBIS) || isTRUE(GBIF)
  if (needs_dl && (is.null(taxa) || is.null(regional_poly)))
    stop("Running a fresh OBIS/GBIF download needs both `taxa` and ",
         "`regional_poly`. Set OBIS/GBIF to FALSE or to a pre-made source ",
         "to skip downloading.")

  cols         <- c("Species", "Source")
  empty        <- data.frame(Species = character(0), Source = character(0),
                             stringsAsFactors = FALSE)
  species_rows <- empty   # binomials: OBIS/GBIF + local non-"sp." rows
  genus_rows   <- empty   # genus-only: local "Genus / sp." rows

  # --- OBIS source (uses the function's return value, not a written file) --
  if (isTRUE(OBIS)) {
    obis_df <- OBIS_download(obis_taxa = taxa, regional_poly = regional_poly,
                            output_dir = NULL, kingdom = kingdom,
                            marine = marine, freshwater = freshwater,
                            terrestrial = terrestrial, brackish = brackish)
    species_rows <- unique(rbind(species_rows, obis_df[, cols, drop = FALSE]))
  } else if (is.character(OBIS)) {
    species_rows <- unique(rbind(species_rows, .regatta_read_source_csv(OBIS, "OBIS")))
  }

  # --- GBIF source --------------------------------------------------------
  if (!.regatta_is_off(GBIF)) {
    if (isTRUE(GBIF)) {
      gbif_df <- GBIF_download(obis_taxa = taxa, regional_poly = regional_poly,
                              output_dir = NULL, kingdom = kingdom,
                              gbif_fill_families = gbif_fill_families)
    } else if (is.character(GBIF) && length(GBIF) == 1 &&
               (grepl("\\.csv$", GBIF, ignore.case = TRUE) || file.exists(GBIF))) {
      gbif_df <- .regatta_read_source_csv(GBIF, "GBIF")          # pre-made CSV
    } else {
      gbif_df <- GBIF_download(existing_download = GBIF, output_dir = NULL)
    }
    species_rows <- unique(rbind(species_rows, gbif_df[, cols, drop = FALSE]))
  }

  # --- Local CSV(s) -------------------------------------------------------
  if (!.regatta_is_off(CSV)) {
    for (lc in CSV) {
      if (!file.exists(lc))
        stop("Cannot find local CSV '", lc, "'. Give a path to an existing file.")
      db <- utils::read.csv(lc, fileEncoding = "latin1", check.names = FALSE)
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
  message(nrow(species_rows), " species-level entries (for the reference DB); ",
          nrow(lca_rows), " entries incl. ", nrow(genus_rows),
          " genus-level (for the checklist LCA).")

  # --- Three return objects ------------------------------------------------
  # for_making_localdb : bare, sorted, unique scientific names -- nothing else.
  #   Written one-per-line, this is the input for FASTA-building tools (CRABS).
  # checklist_summary  : the full taxonomized table, including how each name
  #   resolved (scientific name / synonym / unresolved) -- the audit view.
  # for_LCA            : the taxID column + the 7 ranks only, resolved rows,
  #   cleaned and ready for the LCA step (run_regatta()/reconcile_checklist()).
  ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  for_making_localdb <- sort(unique(species_rows$Species))

  checklist_summary <- NULL
  if (file.exists(sql_path)) {
    message("Taxonomizing the LCA checklist ...")
    checklist_summary <- taxonomize_checklist(lca_rows, sql_path = sql_path)
    resolved_rows <- !is.na(checklist_summary$taxID)
    for_LCA <- checklist_summary[resolved_rows, c("taxID", ranks), drop = FALSE]
    rownames(for_LCA) <- NULL
  } else {
    message("No taxonomizr DB at '", sql_path, "' -- for_LCA is returned as a ",
            "name list and checklist_summary is NULL; reconcile_checklist()/",
            "run_regatta() will taxonomize it (with a warning) at the LCA step.")
    for_LCA <- lca_rows
  }
  lca_taxonomized <- all(ranks %in% names(for_LCA))

  # --- Write only when an output directory is given -----------------------
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    stem <- paste0("comprehensive_", trimws(region), "_", trimws(label), "_list")
    # CRABS input: one scientific name per line, no header, no columns.
    writeLines(for_making_localdb,
               file.path(output_dir, paste0(stem, "_for_making_localdb.txt")))
    if (!is.null(checklist_summary))
      utils::write.csv(checklist_summary,
                       file.path(output_dir, paste0(stem, "_checklist_summary.csv")),
                       row.names = FALSE)
    if (lca_taxonomized)
      saveRDS(for_LCA, file.path(output_dir, paste0(stem, "_for_LCA.rds")))
    else
      readr::write_delim(for_LCA,
                         file.path(output_dir, paste0(stem, "_for_LCA.txt")),
                         delim = "\t")
    message("Wrote checklist files to ", output_dir)
  }

  invisible(list(for_making_localdb = for_making_localdb,
                 checklist_summary  = checklist_summary,
                 for_LCA            = for_LCA))
}

# FALSE / NULL / NA all mean "this source is off". A non-FALSE value (TRUE, a
# path, a key, an object) means "use it".
.regatta_is_off <- function(x) is.null(x) || (is.logical(x) && length(x) == 1 && !isTRUE(x))

# At the LCA step, accept a checklist that may not be taxonomized yet: if it
# lacks the 7 rank columns, warn and taxonomize it on the fly. A path, a
# character vector of names, or a names-only data.frame all work.
.regatta_ensure_taxonomized <- function(checklist, sql_path = "accessionTaxa.sql") {
  ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  if (is.data.frame(checklist) && all(ranks %in% names(checklist)))
    return(checklist)
  warning("checklist was not pre-taxonomized; taxonomizing it now via ",
          "taxonomize_checklist(). Build it with build_regional_checklist() ",
          "to do this ahead of time and silence this warning.", call. = FALSE)
  taxonomize_checklist(checklist, sql_path = sql_path)
}

# Read a pre-made per-source checklist CSV (a fed-in OBIS/GBIF equivalent)
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
