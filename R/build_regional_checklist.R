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
#'   \item `methods` -- a ready-to-adapt methods/provenance sentence naming the
#'     sources used, the date, the polygon, and the citations to include.
#'     Package citations are pulled live from [utils::citation()] and the GBIF
#'     line uses the download's own DOI -- nothing is invented; verify against
#'     your reference manager.
#'   \item `GBIF_download_citation` -- **only when `GBIF = TRUE`** (a fresh
#'     download): the download key + GBIF citation (one-row data.frame).
#'     `NULL` otherwise. Keep it for your paper and to reuse the download via
#'     `GBIF = <download_key>`.
#' }
#' The outputs are always written to the required `output_dir` (named from
#' `region` and `label`) and also returned, so the expensive artifacts persist.
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
#' @param CSV Optional local checklist input. `NULL` (default) adds none --
#'   the public sources (OBIS/GBIF) drive the checklist. Otherwise a character
#'   vector of paths to local `Genus`+`Species` CSV(s). A `Species` value of
#'   `"sp."` keeps the genus alone (no epithet); that genus-level entry goes
#'   into `for_LCA` only, never `for_making_localdb`. If you supply one, prefer
#'   a **citable** source (e.g. a published regional checklist) so your pipeline
#'   stays reproducible -- record its citation in your methods.
#' @param marine,freshwater,terrestrial,brackish OBIS habitat filters, passed
#'   to [OBIS_download()]. Defaults keep marine plus anadromous/catadromous
#'   species and drop only land contaminants.
#' @param output_dir **Required.** Directory to write the outputs into. Each run
#'   gets its own dated subdirectory `<region>_<label>_<Date>` inside it, so
#'   successive runs don't pile up loose files (the GBIF `.zip`, when a download
#'   runs, is kept there too). The expensive artifacts are always saved (so you
#'   don't rebuild them): `comprehensive_<region>_<label>_list_for_making_localdb.txt`
#'   (one name per line), `..._checklist_summary.csv`, `..._for_LCA.rds` (or
#'   `.txt` when not taxonomized), `..._methods.txt`, and
#'   `..._GBIF_download_info.txt` (when a GBIF download ran). The same objects
#'   are also returned.
#' @param sql_path Path to the local `accessionTaxa.sql` taxonomizr DB used to
#'   taxonomize the LCA list. **Default `NULL` skips taxonomization** here
#'   (`for_LCA` is returned as a name list and `checklist_summary` is `NULL`);
#'   it is deferred to the reconcile step ([run_regatta()] /
#'   [reconcile_checklist()]), which taxonomize using their own cached DB. Pass
#'   a path to taxonomize now -- e.g. an existing DB, or the shared cache
#'   `file.path(tools::R_user_dir("REGATTA", "cache"), "accessionTaxa.sql")`. If
#'   the given path is missing, REGATTA prompts to build it (interactive) or
#'   errors with the one-line build command (non-interactive) unless
#'   `overwrite_taxonomy_files = TRUE`. Only names+nodes are needed here, so the
#'   build uses `taxonomizr::prepareDatabase(getAccessions = FALSE)`.
#' @param overwrite_taxonomy_files If `TRUE`, (re)build the taxonomy DB at
#'   `sql_path` even if one exists -- the way to refresh a stale cached snapshot
#'   in place without hunting for the file. Default `FALSE` reuses an existing
#'   DB (reporting its build date and soft-warning if it is old). The DB is a
#'   snapshot of NCBI taxonomy and is never rebuilt silently, to keep results
#'   reproducible.
#' @param kingdom Kingdom used by [resolve_taxa()] to disambiguate query
#'   taxa in the downloaders. Default `"Animalia"`.
#' @param gbif_fill_families Passed to [GBIF_download()] for a fresh GBIF run.
#'
#' @return Invisibly, a list with `for_making_localdb` (a bare character vector
#'   of unique species binomials, for a reference-database builder),
#'   `checklist_summary` (the full taxonomized table with per-name resolution
#'   status, or `NULL` if no DB was available), `for_LCA` (`taxID` + the 7
#'   ranks, ready for the LCA step), `methods` (a citations-included provenance
#'   sentence), and `GBIF_download_citation` (the GBIF download key + citation,
#'   only when `GBIF = TRUE`, else `NULL`). Pass `for_LCA` straight to
#'   [run_regatta()] / [reconcile_checklist()].
#'
#' @examples
#' \dontrun{
#' # OBIS (default) + a local checklist. output_dir is required; the default
#' # cached taxonomy DB is built on first use (prompted) or reused after.
#' cl <- build_regional_checklist(
#'   region        = "galapagos",
#'   label         = "fish",
#'   taxa          = "fish",
#'   regional_poly = "POLYGON ((-92 2, -89 2, -89 -2, -92 -2, -92 2))",
#'   CSV           = "~/other_project/galapagos_fish_checklist.csv",
#'   output_dir    = "my_checklists"
#' )
#' run_regatta(input = "MiFish_obitools.tab", checklist = cl$for_LCA)
#'
#' # Refresh a stale cached taxonomy snapshot in place:
#' build_regional_checklist(region = "galapagos", label = "fish", taxa = "fish",
#'   regional_poly = "POLYGON ((-92 2, -89 2, -89 -2, -92 -2, -92 2))",
#'   output_dir = "my_checklists", overwrite_taxonomy_files = TRUE)
#'
#' # Skip taxonomization now and defer it to the reconcile step:
#' build_regional_checklist(region = "galapagos", label = "fish", taxa = "fish",
#'   regional_poly = "POLYGON ((-92 2, -89 2, -89 -2, -92 -2, -92 2))",
#'   output_dir = "my_checklists", sql_path = NULL)
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
                                     CSV  = NULL,
                                     marine = TRUE,
                                     freshwater = NA,
                                     terrestrial = FALSE,
                                     brackish = NA,
                                     output_dir,
                                     sql_path = NULL,
                                     overwrite_taxonomy_files = FALSE,
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
  if (missing(output_dir) || is.null(output_dir) || !is.character(output_dir) ||
      length(output_dir) != 1 || !nzchar(trimws(output_dir)))
    stop("`output_dir` is required: give a directory to write the checklist ",
         'outputs into, e.g. output_dir = "my_checklists". The expensive ',
         "artifacts (taxonomized checklist, source lists) are always saved so ",
         "you don't have to rebuild them.")

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
  gbif_citation <- NULL   # set ONLY for a fresh GBIF = TRUE download
  sources_used  <- character(0)  # which sources actually contributed, for cl$methods

  # Per-run output directory, so successive runs land in their own dated folder
  # (outputs and the GBIF .zip) instead of piling up loose in output_dir.
  run_dir <- file.path(output_dir,
                       paste0(trimws(region), "_", trimws(label), "_", Sys.Date()))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

  # --- OBIS source (uses the function's return value, not a written file) --
  if (isTRUE(OBIS)) {
    obis_df <- OBIS_download(obis_taxa = taxa, regional_poly = regional_poly,
                            output_dir = NULL, kingdom = kingdom,
                            marine = marine, freshwater = freshwater,
                            terrestrial = terrestrial, brackish = brackish)
    species_rows <- unique(rbind(species_rows, obis_df[, cols, drop = FALSE]))
    sources_used <- c(sources_used, "OBIS")
  } else if (is.character(OBIS)) {
    species_rows <- unique(rbind(species_rows, .regatta_read_source_csv(OBIS, "OBIS")))
    sources_used <- c(sources_used, "OBIS")
  }

  # --- GBIF source --------------------------------------------------------
  if (!.regatta_is_off(GBIF)) {
    if (isTRUE(GBIF)) {
      gbif_df <- GBIF_download(obis_taxa = taxa, regional_poly = regional_poly,
                              output_dir = NULL, kingdom = kingdom,
                              gbif_fill_families = gbif_fill_families,
                              download_dir = run_dir)
      # The citation/key is surfaced only for a fresh download (GBIF = TRUE).
      gbif_citation <- attr(gbif_df, "gbif_download")
    } else if (is.character(GBIF) && length(GBIF) == 1 &&
               (grepl("\\.csv$", GBIF, ignore.case = TRUE) || file.exists(GBIF))) {
      gbif_df <- .regatta_read_source_csv(GBIF, "GBIF")          # pre-made CSV
    } else {
      gbif_df <- GBIF_download(existing_download = GBIF, output_dir = NULL)
    }
    species_rows <- unique(rbind(species_rows, gbif_df[, cols, drop = FALSE]))
    sources_used <- c(sources_used, "GBIF")
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
    sources_used <- c(sources_used, "a local CSV checklist")
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

  # Ensure (or deliberately skip) the taxonomy DB. sql_path = NULL skips
  # taxonomization entirely and defers it to the reconcile step.
  checklist_summary <- NULL
  db_date <- NA_character_
  have_db <- !is.null(sql_path) &&
    .regatta_ensure_taxonomy_db(sql_path, overwrite_taxonomy_files)
  if (have_db) {
    message("Taxonomizing the LCA checklist ...")
    checklist_summary <- taxonomize_checklist(lca_rows, sql_path = sql_path)
    resolved_rows <- !is.na(checklist_summary$taxID)
    for_LCA <- checklist_summary[resolved_rows, c("taxID", ranks), drop = FALSE]
    rownames(for_LCA) <- NULL
    db_date <- .regatta_taxonomy_db_date(sql_path)$date
  } else {
    message("for_LCA is returned as a name list and checklist_summary is NULL; ",
            "reconcile_checklist()/run_regatta() will taxonomize it (with a ",
            "warning) at the LCA step.")
    for_LCA <- lca_rows
  }
  lca_taxonomized <- all(ranks %in% names(for_LCA))

  # --- Methods sentence (real citations only -- pulled live, never invented) -
  methods <- .regatta_methods(region, label, sources_used, regional_poly,
                              gbif_citation, db_date)

  # --- Write the outputs into this run's dated directory ------------------
  {
    stem <- paste0("comprehensive_", trimws(region), "_", trimws(label), "_list")
    # CRABS input: one scientific name per line, no header, no columns.
    writeLines(for_making_localdb,
               file.path(run_dir, paste0(stem, "_for_making_localdb.txt")))
    if (!is.null(checklist_summary))
      utils::write.csv(checklist_summary,
                       file.path(run_dir, paste0(stem, "_checklist_summary.csv")),
                       row.names = FALSE)
    if (lca_taxonomized)
      saveRDS(for_LCA, file.path(run_dir, paste0(stem, "_for_LCA.rds")))
    else
      readr::write_delim(for_LCA,
                         file.path(run_dir, paste0(stem, "_for_LCA.txt")),
                         delim = "\t")
    writeLines(methods, file.path(run_dir, paste0(stem, "_methods.txt")))
    if (!is.null(gbif_citation))
      writeLines(c(
        "GBIF download -- cite this in your paper.",
        "GBIF citation guidelines: https://www.gbif.org/citation-guidelines",
        "",
        paste0("Download key: ", gbif_citation$download_key),
        paste0("DOI: ",          gbif_citation$doi),
        paste0("Created: ",      gbif_citation$created),
        "Citation:",
        gbif_citation$citation),
        file.path(run_dir, paste0(stem, "_GBIF_download_info.txt")))
    message("Wrote checklist files to ", run_dir)
  }

  invisible(list(for_making_localdb     = for_making_localdb,
                 checklist_summary      = checklist_summary,
                 for_LCA                = for_LCA,
                 methods                = methods,
                 GBIF_download_citation = gbif_citation))
}

# Assemble the methods/provenance sentence. CITATIONS ARE REAL: R-package
# citations are pulled live from citation() (the maintainer-supplied entry),
# the GBIF database citation is the DOI GBIF minted for this exact download,
# and nothing is hardcoded as a fabricated author/year/DOI. Where a citation
# cannot be pulled (package not installed, OBIS has no machine citation) the
# string says so and points the user to the canonical source to fill in.
.regatta_methods <- function(region, label, sources_used, regional_poly,
                             gbif_citation, db_date = NA_character_) {
  today <- as.character(Sys.Date())
  src   <- if (length(sources_used)) paste(sources_used, collapse = ", ") else "the supplied sources"

  cite_pkg <- function(p) tryCatch(
    trimws(gsub("\\s+", " ",
                paste(format(utils::citation(p), style = "text"), collapse = " "))),
    error = function(e)
      paste0("the '", p, "' R package (run citation(\"", p, "\") to get the citation)"))

  refs <- c(paste0("REGATTA: ", cite_pkg("REGATTA")))
  if ("OBIS" %in% sources_used) {
    refs <- c(refs, paste0("robis (OBIS R client): ", cite_pkg("robis")))
    refs <- c(refs, paste0("OBIS database: cite per https://obis.org/ ",
                           "(Ocean Biodiversity Information System), accessed ", today))
  }
  if ("GBIF" %in% sources_used) {
    refs <- c(refs, paste0("rgbif (GBIF R client): ", cite_pkg("rgbif")))
    gbif_ref <- if (!is.null(gbif_citation) && !is.na(gbif_citation$citation))
      gbif_citation$citation
    else
      paste0("GBIF download: cite the download DOI per ",
             "https://www.gbif.org/citation-guidelines, accessed ", today)
    refs <- c(refs, paste0("GBIF database download: ", gbif_ref))
  }

  poly_phrase <- if (is.null(regional_poly))
    "The region was defined by the polygon(s) in the supplied source files."
  else
    paste0("The region was defined by the WKT polygon encompassing ", regional_poly, ".")

  tax_phrase <- if (!is.na(db_date))
    paste0(" Taxonomy was resolved against an NCBI taxonomy snapshot built ",
           db_date, " (via the taxonomizr R package).")
  else ""

  paste0(
    "The regional species checklist for ", trimws(region), " (", trimws(label),
    ") was built from ", src, " on ", today,
    " using the R package REGATTA. ", poly_phrase, tax_phrase,
    " Please cite the following (verify against your reference manager): ",
    paste(refs, collapse = " | "))
}

# --- Taxonomy-DB lifecycle helpers --------------------------------------
# REGATTA caches a persistent taxonomizr DB (default under R_user_dir). The DB
# is a SNAPSHOT of NCBI taxonomy and goes stale as NCBI changes; we never
# rebuild it silently (that would make results irreproducible). Instead we
# report its build date, soft-warn when it is old, and rebuild only on
# explicit request (overwrite_taxonomy_files = TRUE) or interactive consent.

.regatta_db_build_cmd <- function(sql_path)
  paste0('taxonomizr::prepareDatabase("', sql_path, '", getAccessions = FALSE)')

# Package-wide default location for the taxonomizr DB: a persistent per-user
# cache, shared across projects/sessions. Used as the default sql_path
# everywhere so build/run/reconcile all look in the same place.
.regatta_default_sql_path <- function()
  file.path(tools::R_user_dir("REGATTA", "cache"), "accessionTaxa.sql")

# Like .regatta_ensure_taxonomy_db(), but for steps where a DB is mandatory
# (resolving classifier input, taxonomizing a raw checklist): error if one
# isn't available afterward instead of returning FALSE.
.regatta_require_taxonomy_db <- function(sql_path, overwrite_taxonomy_files = FALSE) {
  if (!.regatta_ensure_taxonomy_db(sql_path, overwrite_taxonomy_files))
    stop("A taxonomy database is required here but none is available. Re-run ",
         "with overwrite_taxonomy_files = TRUE, or build it first with:\n  ",
         .regatta_db_build_cmd(sql_path), call. = FALSE)
  invisible(TRUE)
}

# Build/refresh the names+nodes-only DB at sql_path and stamp its build date.
.regatta_build_taxonomy_db <- function(sql_path) {
  if (!requireNamespace("taxonomizr", quietly = TRUE))
    stop("Building the taxonomy DB needs the 'taxonomizr' package. Install it ",
         'with install.packages("taxonomizr").')
  dir.create(dirname(sql_path), recursive = TRUE, showWarnings = FALSE)
  message("Building the NCBI taxonomy DB (names/nodes only) at ", sql_path,
          " -- this downloads ~hundreds of MB and takes a few minutes ...")
  taxonomizr::prepareDatabase(sql_path, getAccessions = FALSE)
  writeLines(as.character(Sys.Date()), paste0(sql_path, ".built_on"))
  message("Taxonomy DB ready.")
  invisible(TRUE)
}

# Build date + its provenance: the sidecar we write when WE build the DB, else
# the .sql file's modification date (approximate).
.regatta_taxonomy_db_date <- function(sql_path) {
  side <- paste0(sql_path, ".built_on")
  if (file.exists(side)) {
    d <- tryCatch(trimws(readLines(side, n = 1)), error = function(e) NA_character_)
    return(list(date = d, source = "REGATTA build"))
  }
  if (file.exists(sql_path))
    return(list(date = as.character(as.Date(file.info(sql_path)$mtime)),
                source = "file date (approx)"))
  list(date = NA_character_, source = NA_character_)
}

# Make sure a usable taxonomy DB is present at sql_path. Returns TRUE if one is
# available afterward, FALSE if the user (interactively) declined to build one
# -- in which case taxonomization is deferred to the reconcile step.
.regatta_ensure_taxonomy_db <- function(sql_path, overwrite_taxonomy_files,
                                        stale_days = 180L) {
  if (isTRUE(overwrite_taxonomy_files)) {
    if (file.exists(sql_path)) {
      message("overwrite_taxonomy_files = TRUE: replacing the taxonomy DB at ", sql_path)
      unlink(c(sql_path, paste0(sql_path, ".built_on")))
    }
    return(.regatta_build_taxonomy_db(sql_path))
  }
  if (file.exists(sql_path)) {
    info <- .regatta_taxonomy_db_date(sql_path)
    message("Using NCBI taxonomy snapshot built ", info$date, " (", info$source,
            "). Pass overwrite_taxonomy_files = TRUE to refresh it.")
    bd <- suppressWarnings(as.Date(info$date))
    if (!is.na(bd) && as.numeric(Sys.Date() - bd) > stale_days)
      warning("This taxonomy snapshot is ", as.numeric(Sys.Date() - bd),
              " days old; NCBI taxonomy has likely changed since. For a current ",
              "analysis, rebuild it with overwrite_taxonomy_files = TRUE.",
              call. = FALSE)
    return(TRUE)
  }
  # DB missing, no overwrite flag.
  if (interactive()) {
    ans <- readline(paste0(
      "No taxonomy DB at '", sql_path, "'. Build one now? ",
      "(~hundreds of MB, a few minutes) [y/N]: "))
    if (tolower(trimws(ans)) %in% c("y", "yes"))
      return(.regatta_build_taxonomy_db(sql_path))
    message("Skipping the taxonomy-DB build at your request.")
    return(FALSE)
  }
  stop("No taxonomy DB at '", sql_path, "'. In a non-interactive session, ",
       "build it once with:\n  ", .regatta_db_build_cmd(sql_path),
       "\nor call build_regional_checklist(..., overwrite_taxonomy_files = TRUE) ",
       "to build it here, or point sql_path at an existing DB (or sql_path = ",
       "NULL to skip taxonomization and defer it to the reconcile step).",
       call. = FALSE)
}

# FALSE / NULL / NA all mean "this source is off". A non-FALSE value (TRUE, a
# path, a key, an object) means "use it".
.regatta_is_off <- function(x) is.null(x) || (is.logical(x) && length(x) == 1 && !isTRUE(x))

# At the LCA step, accept a checklist that may not be taxonomized yet: if it
# lacks the 7 rank columns, warn and taxonomize it on the fly. A path, a
# character vector of names, or a names-only data.frame all work.
.regatta_ensure_taxonomized <- function(checklist,
                                        sql_path = .regatta_default_sql_path(),
                                        overwrite_taxonomy_files = FALSE) {
  ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  if (is.data.frame(checklist) && all(ranks %in% names(checklist)))
    return(checklist)
  warning("checklist was not pre-taxonomized; taxonomizing it now via ",
          "taxonomize_checklist(). Build it with build_regional_checklist() ",
          "to do this ahead of time and silence this warning.", call. = FALSE)
  .regatta_require_taxonomy_db(sql_path, overwrite_taxonomy_files)
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
