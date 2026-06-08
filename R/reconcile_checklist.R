# reconcile_checklist.R
# Eldridge Wisely

# reconcile_checklist(): Reconcile a taxonomy table against a regional
# species checklist by walking each row's lineage from species up
# toward domain, finding the lowest rank present in the checklist,
# and NA-ing every rank below that match. This is the REGATTA
# geographic-reconciliation / specificity-downgrade step (originally
# Pass 3 in Validate_local_assignments.R) -- it preserves specificity
# where the regional checklist supports it and downgrades where it
# doesn't, without using percent-identity heuristics.
#
# By default, output_dir = "reconcile_checklist_out" and the three
# return-list elements are written to disk there as CSVs named
#   <output_prefix>_taxonomy_table.csv  (the reconciled $result, strict 8 cols)
#   <output_prefix>_tracking.csv        (per-ASV before/after audit)
#   <output_prefix>_summary.csv         (the stats data frame)
# The directory is created if it doesn't exist. Pass output_dir = NULL
# to disable file writing entirely. The function still returns the same
# list(result, tracking, stats) regardless.
#
# prior_dir / prior_prefix: if a folder named prior_dir exists in the
# working directory at call time AND contains both
#   <prior_prefix>_tracking.csv  and  <prior_prefix>_summary.csv
# (i.e. the outputs of a previous reconcile_global_local() run),
# reconcile_checklist() will READ those files and write AUGMENTED
# versions of tracking + summary into output_dir, combining the
# reconcile_global_local stage columns/rows with the new
# reconcile_checklist columns/rows. The taxonomy_table.csv is ALWAYS
# strict 8 columns (id_col + 7 ranks), regardless of prior detection.
# Defaults: prior_dir = "reconcile_global_local_out",
#           prior_prefix = "reconcile_global_local".
#
# Accepts either:
#   - the $result output of reconcile_global_local() (8 columns:
#     id_col + 7 ranks), or
#   - a standalone classifier output table from resolve_names,
#     resolve_taxids, or parse_sintax (ASV_id + 7 ranks + optional
#     extras -- extras are ignored for the LCA itself).
#
# Output: a list with three elements.
#
#   $result    The reconciled per-ASV taxonomy table. EXACTLY 8
#              columns: id_col + the 7 lowercase rank columns. Same
#              shape as reconcile_global_local()$result so REGATTA
#              functions chain cleanly. Drop-in to phyloseq
#              tax_table() or to a MetabaR MOTU table after joining
#              read counts back.
#
#   $tracking  Per-ASV before-vs-after audit: before_<rank> and
#              after_<rank> at every rank, regatta_match_rank (the
#              rank at which the input matched the checklist, or NA
#              when nothing matched at any rank -- no regional record),
#              and any input metadata columns
#              passed through. One row per ASV in the input.
#
#   $stats     A compact per-step transition headline: total ASVs,
#              assigned before/after, and how many ASVs were unchanged
#              vs downgraded (specificity reduced) vs dropped (no
#              regional record). The per-rank before/after distribution
#              is NOT here -- it is the input vs regatta_result columns
#              of summarize_regatta()'s report.

#' Reconcile a taxonomy table against a regional species checklist
#'
#' The core REGATTA step. Walks each ASV's lineage species -> domain,
#' finds the lowest rank present in the regional checklist, and NAs every
#' rank below that match. Preserves taxonomic specificity where the
#' checklist supports it; downgrades to a higher rank where it does not.
#' Accepts either the `$result` output of [reconcile_global_local()] or a
#' standalone classifier-output taxonomy table.
#'
#' Returns a strict 8-column `$result` (REGATTA exchange format), a per-ASV
#' `$tracking` before/after audit, and a `$stats` summary; optionally
#' writes those as three CSVs. If a `reconcile_global_local` output folder
#' is present, the `$tracking` and `$stats` written to disk are *augmented*
#' versions combining both stages.
#'
#' @param taxonomy_table A data.frame with `id_col` + 7 lowercase rank
#'   columns. Additional metadata columns are passed through to `$tracking`
#'   unless they match `tracking_drop_pattern`.
#' @param checklist A taxonomized regional checklist (output of
#'   [taxonomize_checklist()]) -- a data.frame with the 7 rank columns. A
#'   not-yet-taxonomized checklist (a names-only data.frame, a character
#'   vector, or a path) is accepted too: it is taxonomized on the fly with a
#'   warning, using `sql_path`.
#' @param id_col ASV identifier column name.
#' @param sql_path Path to the local `accessionTaxa.sql` taxonomizr DB, used
#'   only if `checklist` still needs taxonomizing. Defaults to the same
#'   persistent per-user cache as [build_regional_checklist()]
#'   (`tools::R_user_dir("REGATTA", "cache")`). If it is missing when a raw
#'   checklist must be taxonomized, REGATTA prompts to build it (interactive)
#'   or errors with the build command (non-interactive) unless
#'   `overwrite_taxonomy_files = TRUE`.
#' @param overwrite_taxonomy_files If `TRUE`, (re)build the taxonomy DB at
#'   `sql_path` even if one exists. Only consulted when a raw checklist needs
#'   taxonomizing. Default `FALSE`. See [build_regional_checklist()].
#' @param warn_pct_id If `TRUE` (default), warn when `taxonomy_table` has no
#'   usable `pct_id` column -- percent-identity is optional for the geographic
#'   LCA but enables percent-ID filtering of low-confidence (likelier-wrong)
#'   calls, and is required for the two-DB [reconcile_global_local()] comparison.
#'   [run_regatta()] sets this `FALSE` for its internal two-DB step, whose
#'   intermediate result intentionally carries no `pct_id`.
#' @param output_dir Directory path; default `NULL` writes nothing (the
#'   `result`/`tracking`/`stats` list is returned). Supply a directory to also
#'   write the 3 CSVs there.
#' @param output_prefix Filename prefix for the output CSVs. Default
#'   `"reconcile_checklist"`.
#' @param write_summary If `TRUE` (default), also write the `_summary.csv`
#'   alongside the taxonomy-table and tracking CSVs. [run_regatta()] sets this
#'   `FALSE` so a run writes only one top-level summary (`regatta_summary.csv`).
#' @param prior_dir,prior_prefix Optional directory of a prior
#'   [reconcile_global_local()] output to augment the tracking/summary CSVs
#'   with. Default `prior_dir = NULL` (no augmentation); supply the directory
#'   to enable it (only relevant when `output_dir` is also set).
#' @param tracking_drop_pattern Regex matched against input column names;
#'   matching columns are dropped before they enter `$tracking`.
#'
#' @return A list with `result`, `tracking`, and `stats` elements.
#'
#' @importFrom utils write.csv read.csv
#' @export
reconcile_checklist <- function(taxonomy_table,
                                checklist,
                                id_col        = "ASV_id",
                                sql_path      = .regatta_default_sql_path(),
                                overwrite_taxonomy_files = FALSE,
                                warn_pct_id   = TRUE,
                                output_dir    = NULL,
                                output_prefix = "reconcile_checklist",
                                write_summary = TRUE,
                                prior_dir     = NULL,
                                prior_prefix  = "reconcile_global_local",
                                tracking_drop_pattern =
                                  "^(MERGED_sample:|obiclean_|seq_rank|ID_STATUS|DEFINITION)") {
  # tracking_drop_pattern is a regex matched against column names in
  # the input taxonomy_table BEFORE they get carried into $tracking.
  # Default strips obitools per-sample read-count matrices
  # (MERGED_sample:*), the obiclean_* family, seq_rank, ID_STATUS, and
  # DEFINITION. These are not specific to the checklist-LCA decision
  # REGATTA records, and they balloon the tracking CSV to hundreds of
  # irrelevant columns. BEST_MATCH_IDS / BEST_MATCH_TAXIDS are NOT
  # dropped by default because they record the obitools winning
  # accession and taxID. Pass tracking_drop_pattern = NULL (or "")
  # to keep everything.

  ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

  missing_t <- setdiff(c(id_col, ranks), names(taxonomy_table))
  if (length(missing_t) > 0) {
    stop("taxonomy_table is missing required columns: ",
         paste(missing_t, collapse = ", "))
  }
  # pct_id is optional for the geographic LCA, but its absence forecloses two
  # useful things -- warn (suppressed by run_regatta on its two-DB rec$result,
  # which legitimately carries no pct_id).
  if (isTRUE(warn_pct_id) &&
      (!("pct_id" %in% names(taxonomy_table)) || all(is.na(taxonomy_table$pct_id)))) {
    warning(
      "No usable `pct_id` (percent-identity) column in the taxonomy table. ",
      "The geographic reconciliation still runs, but without percent identity ",
      "you cannot apply percent-ID filtering -- which removes low-confidence ",
      "assignments that are more likely to be the wrong species -- and the ",
      "two-database global-vs-local comparison (reconcile_global_local() / the ",
      "run_regatta() two-DB workflow) cannot run. Add a `pct_id` column to your ",
      "taxonomy table where you can.", call. = FALSE)
  }
  # Accept a not-yet-taxonomized checklist: taxonomize on the fly (with a
  # warning) so the LCA still runs. A pre-taxonomized checklist is unchanged.
  checklist <- .regatta_ensure_taxonomized(checklist, sql_path,
                                           overwrite_taxonomy_files)
  missing_c <- setdiff(ranks, names(checklist))
  if (length(missing_c) > 0) {
    stop("checklist is missing required rank columns after taxonomizing: ",
         paste(missing_c, collapse = ", "))
  }

  # Drop obitools / per-sample bookkeeping from the input before it
  # flows into tracking. Essential columns (id, ranks) are preserved.
  if (!is.null(tracking_drop_pattern) && nzchar(tracking_drop_pattern)) {
    essential  <- c(id_col, ranks)
    candidates <- setdiff(names(taxonomy_table), essential)
    to_drop    <- candidates[grepl(tracking_drop_pattern, candidates)]
    if (length(to_drop) > 0) {
      taxonomy_table <- taxonomy_table[, setdiff(names(taxonomy_table), to_drop),
                                       drop = FALSE]
    }
  }

  checklist_sets <- lapply(ranks, function(r) unique(stats::na.omit(checklist[[r]])))
  names(checklist_sets) <- ranks

  before    <- taxonomy_table
  corrected <- taxonomy_table

  # Find the lowest rank at which each row matches the checklist.
  # Walk species -> domain so a more specific match wins.
  match_idx <- rep(NA_integer_, nrow(taxonomy_table))
  for (p in seq(length(ranks), 1, by = -1)) {
    vals <- taxonomy_table[[ranks[p]]]
    hit <- !is.na(vals) & vals %in% checklist_sets[[ranks[p]]] & is.na(match_idx)
    match_idx[hit] <- p
  }

  # NA out ranks below the matched rank; rows with no match anywhere
  # get every rank NA'd (no regional record at any rank).
  for (p in seq_along(ranks)) {
    blank <- is.na(match_idx) | p > match_idx
    corrected[[ranks[p]]][blank] <- NA
  }

  match_rank <- ifelse(is.na(match_idx), NA_character_, ranks[match_idx])

  # --- Build $result: ASV_id + 7 ranks, nothing else ---
  result <- data.frame(placeholder = corrected[[id_col]], stringsAsFactors = FALSE)
  names(result)[1] <- id_col
  for (r in ranks) result[[r]] <- corrected[[r]]
  rownames(result) <- NULL

  # --- Build $tracking: per-ASV before/after audit ---
  # Start from the full input (preserves all metadata columns), then
  # add before/after rank columns and the match rank.
  tracking <- before
  for (r in ranks) {
    tracking[[paste0("before_", r)]] <- before[[r]]
    tracking[[paste0("after_",  r)]] <- corrected[[r]]
  }
  tracking$regatta_match_rank <- match_rank
  # Re-derive scientific_name from corrected lineage (lowest non-NA rank)
  tracking$after_scientific_name <- vapply(seq_len(nrow(corrected)), function(i) {
    vals <- unlist(corrected[i, ranks])
    non_na <- !is.na(vals)
    if (!any(non_na)) return(NA_character_)
    vals[max(which(non_na))]
  }, character(1))
  rownames(tracking) <- NULL

  # --- Build $stats: the before -> after specificity story ---
  # Per-ASV "depth" = the index of its lowest (most specific) non-NA rank, NA
  # if it has no rank at all. "Assigned" = any non-NA rank, so the counts are
  # right even when a pre-resolved input has no `domain` column. The checklist
  # LCA only ever removes specificity, so an ASV is either unchanged, downgraded
  # to a coarser rank, or dropped entirely (no regional record).
  depth_of <- function(df) {
    m <- !is.na(as.matrix(df[, ranks, drop = FALSE]))
    vapply(seq_len(nrow(m)),
           function(i) if (any(m[i, ])) max(which(m[i, ])) else NA_integer_,
           integer(1))
  }
  before_d <- depth_of(before)
  after_d  <- depth_of(corrected)
  asg_b <- !is.na(before_d)
  asg_a <- !is.na(after_d)
  n_before    <- sum(asg_b)
  n_after     <- sum(asg_a)
  n_unchanged <- sum(asg_b & asg_a & after_d == before_d)
  n_downgrade <- sum(asg_b & asg_a & after_d <  before_d)
  n_dropped   <- sum(asg_b & !asg_a)   # was assigned; no regional record at all

  # A compact per-step transition headline. The per-rank before/after
  # *distribution* is NOT duplicated here -- it is read across the input vs
  # regatta_result columns of summarize_regatta()'s report (e.g. the
  # "ID'ed to species" row). This keeps one authoritative report.
  stats <- data.frame(
    metric = c("total ASVs",
               "assigned before checklist-LCA",
               "assigned after checklist-LCA",
               "ASVs unchanged (specificity kept)",
               "ASVs downgraded (specificity reduced)",
               "no regional record (call dropped)"),
    count  = c(nrow(corrected), n_before, n_after,
               n_unchanged, n_downgrade, n_dropped),
    stringsAsFactors = FALSE
  )

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    # taxonomy_table is ALWAYS the strict 8-column result.
    utils::write.csv(result,
                     file.path(output_dir, paste0(output_prefix, "_taxonomy_table.csv")),
                     row.names = FALSE)

    # Default behavior for tracking + summary: just write the reconcile_checklist
    # versions. But if a prior reconcile_global_local output folder is present,
    # read its tracking + summary and AUGMENT before writing.
    tracking_to_write <- tracking
    summary_to_write  <- stats

    prior_tracking_path <- if (!is.null(prior_dir)) file.path(prior_dir, paste0(prior_prefix, "_tracking.csv")) else NULL
    prior_summary_path  <- if (!is.null(prior_dir)) file.path(prior_dir, paste0(prior_prefix, "_summary.csv"))  else NULL
    prior_detected <- !is.null(prior_dir) && dir.exists(prior_dir) &&
                      !is.null(prior_tracking_path) && file.exists(prior_tracking_path) &&
                      !is.null(prior_summary_path)  && file.exists(prior_summary_path)

    if (prior_detected) {
      message("Detected prior reconcile_global_local output at ", normalizePath(prior_dir),
              "; augmenting tracking + summary.")
      prior_tracking <- utils::read.csv(prior_tracking_path, stringsAsFactors = FALSE)
      prior_summary  <- utils::read.csv(prior_summary_path,  stringsAsFactors = FALSE)

      # Pull the new (post-checklist) columns from the freshly-built tracking.
      # Only the columns that are NEW relative to the prior stage:
      new_cols <- c(id_col,
                    paste0("after_",  ranks),
                    "regatta_match_rank",
                    "after_scientific_name")
      new_cols <- intersect(new_cols, names(tracking))
      addons   <- tracking[, new_cols, drop = FALSE]
      # Rename "after_*" to "post_checklist_*" in the augmented view so the
      # two stages are unambiguous when read alongside preferred_* columns
      # from reconcile_global_local.
      rename_map <- c(paste0("after_", ranks),
                      "after_scientific_name")
      new_names  <- c(paste0("post_checklist_", ranks),
                      "post_checklist_scientific_name")
      for (i in seq_along(rename_map)) {
        hit <- which(names(addons) == rename_map[i])
        if (length(hit) == 1) names(addons)[hit] <- new_names[i]
      }
      tracking_to_write <- merge(prior_tracking, addons, by = id_col, all.x = TRUE)

      # Augmented summary: stack reconcile_global_local stats rows on top of
      # reconcile_checklist stats rows, with a stage label so each row is
      # clearly attributed.
      if (all(c("metric", "count") %in% names(prior_summary)) &&
          all(c("metric", "count") %in% names(stats))) {
        summary_to_write <- rbind(
          data.frame(stage = "reconcile_global_local",
                     metric = prior_summary$metric,
                     count  = prior_summary$count,
                     stringsAsFactors = FALSE),
          data.frame(stage = "reconcile_checklist",
                     metric = stats$metric,
                     count  = stats$count,
                     stringsAsFactors = FALSE)
        )
      }
    }

    utils::write.csv(tracking_to_write,
                     file.path(output_dir, paste0(output_prefix, "_tracking.csv")),
                     row.names = FALSE)
    if (isTRUE(write_summary)) {
      utils::write.csv(summary_to_write,
                       file.path(output_dir, paste0(output_prefix, "_summary.csv")),
                       row.names = FALSE)
    }
    message("Wrote ", if (isTRUE(write_summary)) "3" else "2",
            " CSVs to ", normalizePath(output_dir))
  }

  list(result = result, tracking = tracking, stats = stats)
}
