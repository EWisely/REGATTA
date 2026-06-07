# run_regatta.R
# Eldridge Wisely

# High-level REGATTA wrapper. Accepts file paths (single classifier, a
# vsearch LCA + userout pair, or a named list of global + local
# classifier outputs) or a folder of input files. Auto-detects file
# formats, dispatches to the right preprocessor, runs the right
# reconcile* steps, writes 3-CSV-per-stage output triples plus a
# top-level 21-row summary, and emits a run_log.txt describing what
# was detected and run.

# ---- internal: format detection ----------------------------------

#' Detect the format of a single classifier-output file (internal)
#' @keywords internal
#' @noRd
.regatta_detect_format <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") return("rds")
  if (ext == "csv") {
    headers <- strsplit(readLines(path, n = 1), ",", fixed = TRUE)[[1]]
    if ("BestTaxon" %in% headers) return("besttaxon")
    return("unknown_csv")
  }
  if (ext %in% c("tab", "txt")) {
    first <- readLines(path, n = 1)
    headers <- strsplit(first, "\t", fixed = TRUE)[[1]]
    if ("TAXID" %in% headers && "BEST_IDENTITY" %in% headers) {
      return("obitools_tab")
    }
    # Scan a window of lines to find the first non-empty col 2 -- the
    # LCA file has many rows where unassigned ASVs leave col 2 blank.
    lines <- readLines(path, n = 30)
    for (ln in lines) {
      fields <- strsplit(ln, "\t", fixed = TRUE)[[1]]
      if (length(fields) >= 2 && nzchar(fields[2])) {
        if (length(fields) == 6 && grepl(";tax=", fields[2])) {
          return("vsearch_userout")
        }
        if (length(fields) == 2 && grepl("^[dpcofgs]:", fields[2])) {
          return("vsearch_lca")
        }
        break
      }
    }
  }
  "unknown"
}

# ---- internal: input normalization -------------------------------

#' Normalize the `input` arg into {role, classifier paths} (internal)
#' @keywords internal
#' @noRd
.regatta_normalize_input <- function(input) {
  if (is.character(input)) {
    if (length(input) == 1) {
      if (dir.exists(input)) return(.regatta_scan_folder(input))
      return(list(role = "single", classifier = list(input)))
    }
    if (length(input) == 2) {
      fmts <- vapply(input, .regatta_detect_format, character(1))
      if (all(fmts %in% c("vsearch_lca", "vsearch_userout")) &&
          length(unique(fmts)) == 2) {
        return(list(role = "single", classifier = list(input)))
      }
      stop("A vector of two paths must be a vsearch LCA + userout pair. ",
           "Detected formats: ", paste(fmts, collapse = ", "),
           ". For a two-classifier workflow, use ",
           "list(global = ..., local = ...).")
    }
    stop("Unsupported input: character vector of length ", length(input))
  }
  if (is.list(input)) {
    if (!all(c("global", "local") %in% names(input))) {
      stop("Named-list input must contain `global` and `local` entries.")
    }
    return(list(role = "two-DB", global = input$global, local = input$local))
  }
  stop("`input` must be a path, a folder, a length-2 vector, ",
       "or list(global=, local=).")
}

#' Scan a folder for REGATTA inputs (internal)
#' @keywords internal
#' @noRd
.regatta_scan_folder <- function(folder) {
  files <- list.files(folder, full.names = TRUE)
  files <- files[!file.info(files)$isdir]
  if (length(files) == 0) stop("Input folder is empty: ", folder)
  base  <- tolower(basename(files))
  fmts  <- vapply(files, .regatta_detect_format, character(1))
  is_g  <- grepl("^global", base)
  is_l  <- grepl("^local",  base)
  if (any(is_g) || any(is_l)) {
    if (!any(is_g) || !any(is_l)) {
      stop("Folder has files marked global or local but not both: ",
           "found global=", sum(is_g), " local=", sum(is_l))
    }
    return(list(role = "two-DB",
                global = files[is_g], local = files[is_l]))
  }
  has_lca     <- any(fmts == "vsearch_lca")
  has_userout <- any(fmts == "vsearch_userout")
  if (has_lca && has_userout) {
    return(list(role = "single",
                classifier = list(files[fmts %in% c("vsearch_lca","vsearch_userout")])))
  }
  if (length(files) == 1) {
    return(list(role = "single", classifier = list(files)))
  }
  stop("Folder contains multiple files without clear roles. Found:\n",
       paste0("  ", basename(files), " (", fmts, ")", collapse = "\n"),
       "\nRename to global.* and local.* in the folder, ",
       "or pass an explicit list(global=, local=).")
}

# ---- internal: preprocessor dispatch -----------------------------

#' Collapse a long (per-sample) classifier table to one row per ASV (internal)
#'
#' BestTaxon-style outputs are often long: one row per (sample, ASV), so a
#' single ASV's taxonomy is repeated across every sample it appears in (with
#' per-sample columns like read counts varying). Taxonomy reconciliation is
#' per-ASV, so we collapse to one row per `id_col`, keeping only the columns
#' that are constant within each ASV (taxonomy, sequence length, etc.) and
#' dropping the ones that vary across samples (Sample_name, nReads, ...).
#' @keywords internal
#' @noRd
.regatta_dedupe_by_id <- function(df, id_col) {
  if (!anyDuplicated(df[[id_col]])) return(df)   # already one row per ASV
  keep <- vapply(names(df), function(cn) {
    if (cn == id_col) return(TRUE)
    # keep only columns functionally determined by id_col (constant per ASV)
    !any(tapply(df[[cn]], df[[id_col]],
                function(x) length(unique(x)) > 1L))
  }, logical(1))
  unique(df[, keep, drop = FALSE])
}

#' Ensure a per-ASV id column exists (internal)
#'
#' Fallback for classifier outputs that carry no `id_col` and where the user
#' did not name one: synthesize an id from the NCBI taxID that resolve_names /
#' resolve_taxids attach (`taxID_<id>`), falling back to a row index
#' (`row_<n>`) only where the taxID is NA (unresolved). When `id_col` already
#' exists (the common case -- obitools `ID`, vsearch `ASV_id`, or a user
#' `Hash`), the table is returned unchanged.
#' @keywords internal
#' @noRd
.regatta_ensure_id_col <- function(tax, id_col = "ASV_id") {
  if (id_col %in% names(tax)) return(tax)
  new_id <- if ("taxID" %in% names(tax)) {
    ifelse(is.na(tax$taxID), paste0("row_", seq_len(nrow(tax))),
           paste0("taxID_", tax$taxID))
  } else {
    paste0("row_", seq_len(nrow(tax)))
  }
  cbind(stats::setNames(list(new_id), id_col), tax, stringsAsFactors = FALSE)
}

#' Read a classifier file or vsearch pair into a tax table (internal)
#' @keywords internal
#' @noRd
.regatta_read_input <- function(paths, sql_path, id_col = "ASV_id") {
  if (length(paths) == 1) {
    fmt <- .regatta_detect_format(paths)
    if (fmt == "obitools_tab") {
      raw <- utils::read.delim(paths, sep = "\t", check.names = FALSE,
                               stringsAsFactors = FALSE)
      tax <- resolve_taxids(raw, taxid_col = "TAXID", sql_path = sql_path)
      tax$pct_id <- suppressWarnings(as.numeric(raw$BEST_IDENTITY))
      names(tax)[names(tax) == "ID"] <- "ASV_id"
      return(list(table = .regatta_ensure_id_col(tax, id_col),
                  format = "obitools_tab"))
    }
    if (fmt == "vsearch_userout") {
      return(list(table = .regatta_ensure_id_col(parse_vsearch_userout(paths), id_col),
                  format = fmt))
    }
    if (fmt == "vsearch_lca") {
      lca <- utils::read.delim(paths, sep = "\t", header = FALSE,
                               stringsAsFactors = FALSE,
                               col.names = c("ASV_id", "sintax"))
      tax <- parse_sintax(lca, sintax_col = "sintax")
      tax$pct_id <- NA_real_
      return(list(table = .regatta_ensure_id_col(tax, id_col), format = fmt))
    }
    if (fmt == "besttaxon") {
      raw <- utils::read.csv(paths, stringsAsFactors = FALSE)
      # If the file carries the id column (e.g. a sequence Hash), collapse a
      # long sample x ASV table to one row per ASV before resolving. If it
      # carries no id column, skip dedupe and let .regatta_ensure_id_col()
      # synthesize one from the resolved taxID after lookup.
      if (id_col %in% names(raw)) {
        n_in <- nrow(raw)
        raw  <- .regatta_dedupe_by_id(raw, id_col)
        if (nrow(raw) < n_in) {
          message("Collapsed long BestTaxon table from ", n_in, " rows to ",
                  nrow(raw), " unique ASVs (by ", id_col, ").")
        }
      }
      tax <- resolve_names(raw, name_col = "BestTaxon", sql_path = sql_path)
      return(list(table = .regatta_ensure_id_col(tax, id_col), format = fmt))
    }
    stop("Could not detect format of input file: ", paths)
  }
  if (length(paths) == 2) {
    fmts <- vapply(paths, .regatta_detect_format, character(1))
    if (all(fmts %in% c("vsearch_lca", "vsearch_userout"))) {
      lca     <- paths[fmts == "vsearch_lca"][1]
      userout <- paths[fmts == "vsearch_userout"][1]
      return(list(table = .regatta_ensure_id_col(
                            parse_vsearch_results(lca_path = lca,
                                                  userout_path = userout), id_col),
                  format = "vsearch_pair"))
    }
    stop("Two-path input must be a vsearch LCA + userout pair (got ",
         paste(fmts, collapse = " / "), ").")
  }
  stop("Unsupported number of paths: ", length(paths))
}

# ---- public: run_regatta --------------------------------------------

#' Run the full REGATTA pipeline from file paths
#'
#' One-call wrapper around [reconcile_global_local()], [reconcile_checklist()],
#' and [summarize_regatta()]. Accepts file paths or a folder, auto-detects
#' classifier format (obitools `.tab`, vsearch userout, vsearch LCA + userout,
#' BestTaxon `.csv`), dispatches to the right preprocessor, runs the
#' appropriate reconcile* steps, and writes the standard 3-CSVs-per-stage
#' output triples plus a top-level 21-row summary and a `run_log.txt`
#' describing what was detected.
#'
#' @param input One of:
#'   * a single file path (single-DB workflow);
#'   * a length-2 character vector of vsearch `lca` + `userout` paths
#'     (single-DB workflow using the LCA-resolved taxonomy);
#'   * a folder path (the function scans the folder -- see Details for the
#'     convention);
#'   * `list(global = ..., local = ...)` for an explicit two-DB workflow
#'     where each entry is a single path or a length-2 vsearch pair.
#' @param checklist Path to a taxonomized regional checklist (`.rds` or
#'   `.csv`), or a data.frame with the 7 rank columns.
#' @param out_dir Output directory. Default `NULL` returns the results
#'   without writing anything; supply a directory to also write the per-stage
#'   CSV triples, the 21-row summary, and `run_log.txt` there (created if
#'   missing).
#' @param sql_path Path to the local `accessionTaxa.sql` taxonomizr DB
#'   (used for obitools / BestTaxon input shapes).
#' @param id_col Name of the per-ASV identifier column. Default
#'   `"ASV_id"` (matches the obitools and vsearch preprocessors). For
#'   BestTaxon/Kraken2-style CSVs whose identifier is named differently
#'   (e.g. `"Hash"`), set this so the wrapper can find the id and, if the
#'   table is long (one row per sample x ASV), collapse it to one row per
#'   ASV before resolving. If a classifier output has **no** id column at
#'   all and you name none, the wrapper synthesizes one from the resolved
#'   NCBI taxID (`taxID_<id>`, or `row_<n>` where unresolved), so the
#'   pipeline still runs -- though naming a real id column is preferable
#'   when one exists, since only that path collapses long tables.
#' @param Local_advantage TRUE (default): local wins ties at the
#'   `best_pctid` step in [reconcile_global_local()].
#'
#' @details
#' **Role is always declared, never inferred from filenames.** Single-input
#' calls don't need a role. Two-input workflows must use either a named
#' list or a folder whose files start with `global.` and `local.`. The
#' wrapper never guesses based on contains-"obi" / contains-"vs" heuristics.
#'
#' **Folder convention:**
#' \itemize{
#'   \item One classifier file in the folder -> single-DB workflow.
#'   \item Two vsearch files (LCA + userout) in the folder -> single-DB
#'     workflow using the joined LCA taxonomy + userout pct_id.
#'   \item Files starting with `global.*` and `local.*` -> two-DB workflow.
#'     Either side may itself be a vsearch LCA + userout pair (matching
#'     `global_lca*` + `global_userout*` etc.).
#' }
#'
#' @return Invisibly a list with the final tables and the
#'   [summarize_regatta()] data.frame. Writes CSV outputs and `run_log.txt`
#'   under `out_dir` only when `out_dir` is supplied.
#'
#' @importFrom utils read.csv read.delim write.csv
#' @export
run_regatta <- function(input,
                        checklist,
                        out_dir         = NULL,
                        sql_path        = "accessionTaxa.sql",
                        id_col          = "ASV_id",
                        Local_advantage = TRUE) {
  t_start <- Sys.time()
  norm <- .regatta_normalize_input(input)

  # out_dir = NULL: return results only, write nothing. Supply a directory to
  # also persist the per-stage CSV triples + summary + run_log there.
  if (!is.null(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  .sub <- function(name) if (is.null(out_dir)) NULL else file.path(out_dir, name)

  if (is.character(checklist) && length(checklist) == 1) {
    cl <- if (grepl("\\.rds$", checklist, ignore.case = TRUE)) readRDS(checklist)
          else utils::read.csv(checklist, stringsAsFactors = FALSE)
  } else if (is.data.frame(checklist)) {
    cl <- checklist
  } else {
    stop("checklist must be a path to .rds/.csv or a data.frame.")
  }
  # If a raw (un-taxonomized) checklist was supplied, taxonomize it now with a
  # warning; a pre-taxonomized checklist passes through unchanged.
  cl <- .regatta_ensure_taxonomized(cl, sql_path)

  dest <- if (is.null(out_dir)) "(return only; no files written)" else normalizePath(out_dir)
  log_lines <- c(
    paste0("REGATTA run started: ", format(t_start)),
    paste0("Workflow: ", norm$role),
    paste0("Output dir: ", dest)
  )
  message("REGATTA: ", norm$role, " workflow -> ", dest)

  if (norm$role == "two-DB") {
    g <- .regatta_read_input(norm$global, sql_path, id_col = id_col)
    log_lines <- c(log_lines, paste0("global (", g$format, "): ",
                                     nrow(g$table), " ASVs"))
    l <- .regatta_read_input(norm$local, sql_path, id_col = id_col)
    log_lines <- c(log_lines, paste0("local  (", l$format, "): ",
                                     nrow(l$table), " ASVs"))

    rec <- reconcile_global_local(
      g$table, l$table, id_col = id_col, Local_advantage = Local_advantage,
      output_dir    = .sub("reconcile_global_local"),
      output_prefix = "reconcile_global_local")
    log_lines <- c(log_lines, "reconcile_global_local done")

    post <- reconcile_checklist(
      rec$result, cl, id_col = id_col,
      output_dir    = .sub("reconcile_checklist"),
      output_prefix = "reconcile_checklist",
      prior_dir     = .sub("reconcile_global_local"),
      prior_prefix  = "reconcile_global_local")
    log_lines <- c(log_lines, "reconcile_checklist done")

    summ <- summarize_regatta(rec, post, g$table, l$table)
    result <- list(global_tax = g$table, local_tax = l$table,
                   reconciled = rec, post_checklist = post, summary = summ)
  } else {
    c1 <- .regatta_read_input(norm$classifier[[1]], sql_path, id_col = id_col)
    log_lines <- c(log_lines, paste0("classifier (", c1$format, "): ",
                                     nrow(c1$table), " ASVs"))
    post <- reconcile_checklist(
      c1$table, cl, id_col = id_col,
      output_dir    = .sub("reconcile_checklist"),
      output_prefix = "reconcile_checklist",
      prior_dir     = NULL)
    log_lines <- c(log_lines, "reconcile_checklist done")

    summ <- summarize_regatta(post_checklist = post)
    result <- list(input_tax = c1$table, post_checklist = post, summary = summ)
  }

  log_lines <- c(log_lines,
                 paste0("Total elapsed: ", format(Sys.time() - t_start)))
  if (!is.null(out_dir)) {
    utils::write.csv(summ, file.path(out_dir, "regatta_summary.csv"),
                     row.names = FALSE)
    writeLines(log_lines, file.path(out_dir, "run_log.txt"))
    message("Wrote ", normalizePath(file.path(out_dir, "regatta_summary.csv")))
  }

  invisible(result)
}
