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
    headers <- trimws(gsub('^"|"$', "", headers))   # strip surrounding quotes
    if ("BestTaxon" %in% headers) return("besttaxon")
    # A pre-resolved taxonomy table: already carries the rank columns, so no
    # name/taxID lookup is needed. Detect on the six unambiguous lower ranks
    # (case-insensitive); the top rank (domain/kingdom/superkingdom) and an id
    # column are normalized/synthesized downstream.
    low <- tolower(headers)
    if (all(c("phylum", "class", "order", "family", "genus", "species") %in% low))
      return("ranks_csv")
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

#' Normalize taxonomy rank column names to the canonical lowercase 7 (internal)
#'
#' Renames case-variant rank columns (e.g. `Class`) to lowercase, maps a
#' `kingdom`/`superkingdom` column to `domain` when no `domain` is present, and
#' adds any of the 7 ranks that are missing as all-NA, so a pre-resolved
#' taxonomy CSV drops straight into [reconcile_checklist()].
#' @keywords internal
#' @noRd
.regatta_normalize_rank_cols <- function(df) {
  ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  low <- tolower(trimws(names(df)))
  if (!("domain" %in% low)) {
    alt <- which(low %in% c("kingdom", "superkingdom"))
    if (length(alt)) low[alt[1]] <- "domain"
  }
  for (r in ranks) {
    hit <- which(low == r)
    if (length(hit)) names(df)[hit[1]] <- r
  }
  for (r in ranks) if (!(r %in% names(df))) df[[r]] <- NA_character_
  df
}

#' Ensure a per-ASV id column exists and is populated (internal)
#'
#' Synthesizes an id when `id_col` is absent OR present-but-empty (all NA/blank
#' -- e.g. an empty `ASV_id`/`MOTU` column): an id from the NCBI taxID when one
#' is attached (`taxID_<id>`), else a row index (`row_<n>`). A populated
#' `id_col` (the common case -- obitools `ID`, vsearch `ASV_id`, a user `Hash`)
#' is returned unchanged.
#' @keywords internal
#' @noRd
.regatta_ensure_id_col <- function(tax, id_col = "ASV_id") {
  present  <- id_col %in% names(tax)
  populated <- present &&
    !all(is.na(tax[[id_col]]) | !nzchar(trimws(as.character(tax[[id_col]]))))
  if (populated) return(tax)
  new_id <- if ("taxID" %in% names(tax)) {
    ifelse(is.na(tax$taxID), paste0("row_", seq_len(nrow(tax))),
           paste0("taxID_", tax$taxID))
  } else {
    paste0("row_", seq_len(nrow(tax)))
  }
  message("No usable '", id_col, "' values found; synthesized per-ASV ids ",
          "(", if ("taxID" %in% names(tax)) "taxID_*/row_*" else "row_*", ").")
  if (present) { tax[[id_col]] <- new_id; tax }   # fill the empty column in place
  else cbind(stats::setNames(list(new_id), id_col), tax, stringsAsFactors = FALSE)
}

#' Read a classifier file or vsearch pair into a tax table (internal)
#' @keywords internal
#' @noRd
.regatta_read_input <- function(paths, sql_path, id_col = "ASV_id",
                                overwrite_taxonomy_files = FALSE) {
  if (length(paths) == 1) {
    fmt <- .regatta_detect_format(paths)
    if (fmt == "obitools_tab") {
      raw <- utils::read.delim(paths, sep = "\t", check.names = FALSE,
                               stringsAsFactors = FALSE)
      .regatta_require_taxonomy_db(sql_path, overwrite_taxonomy_files)
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
      .regatta_require_taxonomy_db(sql_path, overwrite_taxonomy_files)
      tax <- resolve_names(raw, name_col = "BestTaxon", sql_path = sql_path)
      return(list(table = .regatta_ensure_id_col(tax, id_col), format = fmt))
    }
    if (fmt == "ranks_csv") {
      # Already-resolved taxonomy table: normalize the rank column names to the
      # canonical lowercase 7, keep any pct_id/passthrough columns, and ensure a
      # per-ASV id (synthesized if the id column is missing or empty). No DB.
      raw <- utils::read.csv(paths, stringsAsFactors = FALSE, check.names = FALSE)
      # A leading unnamed column (e.g. a written row-index, which read.csv names
      # ""): drop it if it's empty, or recover it as the id if it carries values
      # and no id column is named -- so it doesn't ride along as a junk column.
      blank <- !nzchar(trimws(names(raw)))
      empty <- vapply(raw, function(v)
        all(is.na(v) | !nzchar(trimws(as.character(v)))), logical(1))
      if (any(blank & empty)) raw <- raw[, !(blank & empty), drop = FALSE]
      blank <- !nzchar(trimws(names(raw)))
      if (!(id_col %in% names(raw)) && any(blank))
        names(raw)[which(blank)[1]] <- id_col
      tax <- .regatta_normalize_rank_cols(raw)
      return(list(table = .regatta_ensure_id_col(tax, id_col),
                  format = "ranks_csv"))
    }
    # Unrecognized: tell the user exactly what we saw and what we accept.
    hdr <- tryCatch(
      trimws(gsub('^"|"$', "",
                  strsplit(readLines(paths, n = 1),
                           if (tolower(tools::file_ext(paths)) == "csv") "," else "\t",
                           fixed = TRUE)[[1]])),
      error = function(e) character(0))
    stop("Could not detect the format of '", paths, "'.\n",
         if (length(hdr)) paste0("  Columns found: ", paste(hdr, collapse = ", "), "\n") else "",
         "  REGATTA accepts:\n",
         "    - obitools .tab (a TAXID + BEST_IDENTITY column),\n",
         "    - a vsearch SINTAX lca file, or an lca + userout pair,\n",
         "    - a Kraken2/BestTaxon CSV (a 'BestTaxon' column),\n",
         "    - a pre-resolved taxonomy CSV with the 7 rank columns ",
         "(domain, phylum, class, order, family, genus, species; ",
         "case-insensitive, plus an optional id column).\n",
         "  This looks like a rank table missing one or more of phylum/class/",
         "order/family/genus/species -- check those column names.",
         call. = FALSE)
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
#'   * a single file path (single-DB workflow). The format is auto-detected:
#'     obitools `.tab` (TAXID + BEST_IDENTITY), a vsearch SINTAX `lca` file, a
#'     Kraken2/BestTaxon CSV (a `BestTaxon` column), or a **pre-resolved
#'     taxonomy CSV** that already has the 7 rank columns (domain, phylum,
#'     class, order, family, genus, species; case-insensitive, with an optional
#'     id column -- a missing or empty one is synthesized). An unrecognized file
#'     errors with the columns it found and the formats accepted;
#'   * a length-2 character vector of vsearch `lca` + `userout` paths
#'     (single-DB workflow using the LCA-resolved taxonomy);
#'   * a folder path (the function scans the folder -- see Details for the
#'     convention);
#'   * `list(global = ..., local = ...)` for an explicit two-DB workflow
#'     where each entry is a single path or a length-2 vsearch pair.
#' @param checklist Path to a taxonomized regional checklist (`.rds` or
#'   `.csv`), or a data.frame with the 7 rank columns.
#' @param out_dir **Required.** Output directory. Each run writes into its own
#'   dated `<region>_<label>_<Date>` subfolder of it, so successive runs don't
#'   pile up. That subfolder gets the per-stage CSV triples, the 21-row summary,
#'   and `run_log.txt`. The results list is also returned.
#' @param region,label **Required.** Region and taxon-group labels for this run.
#'   REGATTA filters/adjusts the taxonomy against the regional checklist built
#'   for that region x group, so the run is labeled by them -- pass the same
#'   values you gave [build_regional_checklist()]. They name the dated run
#'   subfolder `<region>_<label>_<Date>`.
#' @param sql_path Path to the local `accessionTaxa.sql` taxonomizr DB, used to
#'   resolve obitools / BestTaxon input and to taxonomize a raw `checklist`.
#'   Defaults to the same persistent per-user cache as
#'   [build_regional_checklist()] (`tools::R_user_dir("REGATTA", "cache")`). If
#'   it is missing when needed, REGATTA prompts to build it (interactive) or
#'   errors with the build command (non-interactive) unless
#'   `overwrite_taxonomy_files = TRUE`. Pre-resolved input (vsearch) and a
#'   pre-taxonomized `checklist` need no DB. (Note: accession-based input would
#'   need the full `taxonomizr` build with `accession2taxid`; the auto-build is
#'   names+nodes only.)
#' @param overwrite_taxonomy_files If `TRUE`, (re)build the taxonomy DB at
#'   `sql_path` even if one exists. Default `FALSE`. See
#'   [build_regional_checklist()].
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
#'   [summarize_regatta()] data.frame. Always also writes the CSV outputs and
#'   `run_log.txt` into the dated run subfolder of `out_dir`.
#'
#' @importFrom utils read.csv read.delim write.csv
#' @export
run_regatta <- function(input,
                        checklist,
                        out_dir,
                        region,
                        label,
                        sql_path        = .regatta_default_sql_path(),
                        overwrite_taxonomy_files = FALSE,
                        id_col          = "ASV_id",
                        Local_advantage = TRUE) {
  t_start <- Sys.time()
  if (missing(region) || missing(label) ||
      !is.character(region) || !is.character(label) ||
      length(region) != 1 || length(label) != 1 ||
      !nzchar(trimws(region)) || !nzchar(trimws(label)))
    stop("Provide a non-empty `region` and `label`. REGATTA filters/adjusts ",
         "the taxonomy against the regional checklist built for that ",
         "region x group, so the run is labeled by them. e.g. ",
         'region = "galapagos", label = "fish".')
  if (missing(out_dir) || is.null(out_dir) || !is.character(out_dir) ||
      length(out_dir) != 1 || !nzchar(trimws(out_dir)))
    stop("`out_dir` is required: give a directory to write the run outputs ",
         'into, e.g. out_dir = "regatta_out". Outputs land in a dated ',
         "<region>_<label>_<Date> subfolder of it.")
  norm <- .regatta_normalize_input(input)

  # Each run lands in its own dated <region>_<label>_<Date> subfolder.
  run_dir <- file.path(out_dir,
                       paste0(trimws(region), "_", trimws(label), "_", Sys.Date()))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  .sub <- function(name) file.path(run_dir, name)

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
  cl <- .regatta_ensure_taxonomized(cl, sql_path, overwrite_taxonomy_files)

  dest <- normalizePath(run_dir)
  log_lines <- c(
    paste0("REGATTA run started: ", format(t_start)),
    paste0("Workflow: ", norm$role),
    paste0("Output dir: ", dest)
  )
  message("REGATTA: ", norm$role, " workflow -> ", dest)

  if (norm$role == "two-DB") {
    g <- .regatta_read_input(norm$global, sql_path, id_col = id_col,
                             overwrite_taxonomy_files = overwrite_taxonomy_files)
    log_lines <- c(log_lines, paste0("global (", g$format, "): ",
                                     nrow(g$table), " ASVs"))
    l <- .regatta_read_input(norm$local, sql_path, id_col = id_col,
                             overwrite_taxonomy_files = overwrite_taxonomy_files)
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
    c1 <- .regatta_read_input(norm$classifier[[1]], sql_path, id_col = id_col,
                              overwrite_taxonomy_files = overwrite_taxonomy_files)
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
  utils::write.csv(summ, file.path(run_dir, "regatta_summary.csv"),
                   row.names = FALSE)
  writeLines(log_lines, file.path(run_dir, "run_log.txt"))
  message("Wrote run outputs to ", normalizePath(run_dir))

  invisible(result)
}
