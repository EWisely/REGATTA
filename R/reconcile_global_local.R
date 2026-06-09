# reconcile_global_local.R
# Eldridge Wisely

# Reconcile two taxonomy tables that classify the same set of ASVs
# against a GLOBAL reference database (e.g. obitools + EMBL/NCBI) and
# a LOCAL reference database (e.g. vsearch + a regionally-curated DB).
#
# The methodological core of REGATTA: when both DBs return an
# assignment, the one with HIGHER percent identity to its reference is
# treated as the better-supported call. When the global DB wins by
# %ID and the local DB also has any assignment, the two lineages are
# LCA'd so that disagreements at lower ranks (e.g. global says
# Sebastes mystinus, local says S. paucispinis) downgrade to the
# shared ancestor (Sebastes) rather than silently taking global's call.
#
# Two reconciliation steps, named descriptively:
#
#   best_pctid              Per ASV, pick the database whose pct_id
#                           to its reference sequence is higher. Ties
#                           are broken by Local_advantage.
#
#   global_lca_to_local     For ASVs where global won best_pctid AND
#                           the local DB also returned an assignment,
#                           replace the preferred lineage with the LCA
#                           of the two and label the database column
#                           "global_lca_to_local". $stats reports how many
#                           were triggered, and -- of those that actually
#                           reduced rank -- a "downgraded: <from> -> <to>"
#                           breakdown (the rest agreed to global's depth).
#
# (The original Validate_local_assignments.R called these Pass 1 and
# Pass 2. The descriptive names make the methods write-up cleaner.)
#
# pct_id scale: the global classifier's pct_id may be on a 0-1 scale
# (e.g. obitools BEST_IDENTITY) while the local classifier's is
# typically 0-100 (e.g. vsearch pctid). The function auto-detects each
# input's scale (if max <= 1, treats as 0-1 and rescales to 0-100) or
# accepts an explicit scale per input. Comparison is on the common
# 0-100 scale.

# Inputs:
#   global_table, local_table  Taxonomy tables: id_col + 7 lowercase
#                              rank columns + a percent-identity
#                              column. Other columns flow through to
#                              $tracking unchanged.
#   id_col                     ASV identifier column name (default
#                              "ASV_id"), same in both tables.
#   global_pct_id_col          Percent-identity column name in
#                              global_table (default "pct_id").
#   local_pct_id_col           Percent-identity column name in
#                              local_table (default "pct_id").
#   global_pct_id_scale        "auto" (default), "0-1", or "0-100".
#   local_pct_id_scale         Same options for local_table.
#   Local_advantage            TRUE (default): local wins ties. FALSE:
#                              global wins ties. Matches the original
#                              Validate_local_assignments.R default.
#   output_dir                 Directory path. Defaults to
#                              "reconcile_global_local_out" (relative to
#                              the working directory). Writes 3 CSVs
#                              there: <output_prefix>_taxonomy_table.csv
#                              (the $result table),
#                              <output_prefix>_tracking.csv (the
#                              $tracking audit), and
#                              <output_prefix>_summary.csv (the $stats
#                              counts). The directory is created if it
#                              does not already exist. Pass
#                              output_dir = NULL to disable file writing
#                              entirely.
#   output_prefix              Filename prefix for the 3 CSVs above.
#                              Default "reconcile_global_local".

# Output: a list with three elements.
#
#   $result    The reconciled per-ASV taxonomy table: id_col, then
#              `pct_id` (the winning database's percent identity, carried
#              right after the id so downstream percent-ID filtering still
#              has it), then the 7 lowercase rank columns. Drop-in to
#              any downstream REGATTA function (regatta_checklist_lca),
#              into phyloseq tax_table(), or into a MetabaR MOTU
#              table after joining read counts back. No bookkeeping
#              columns. By design.
#
#   $tracking  best_ID_combined-style per-ASV decision record, columns
#              grouped per field as <field>_global / <field>_local /
#              <field>_regatta (the reconciled/preferred value), so each
#              of the 7 ranks and pct_id read global / local / regatta
#              side by side. Then the decision columns -- best_pctid
#              winner, preferred_database, preferred_scientific_name,
#              global_lca_to_local_triggered -- followed by any remaining
#              classifier columns from either input (also _global /
#              _local suffixed). One row per ASV in the union of both
#              inputs.
#
#   $stats    Counts of best_pctid winners, global_lca_to_local
#              triggers, and final per-database categories. A small
#              data.frame of (metric, count) rows for diagnostics.

#' Reconcile global-DB vs. local-DB classifier outputs
#'
#' Given two taxonomy tables for the same set of ASVs -- one from a global
#' reference database (e.g. obitools + NCBI/EMBL) and one from a local
#' curated database (e.g. vsearch + a regional reference) -- pick a per-ASV
#' preferred lineage. The function compares percent identity across the two
#' databases (the `best_pctid` step) and, when global wins but local also
#' has an assignment, takes the LCA of the two (the `global_lca_to_local` step);
#' `$stats` reports how many were triggered and, of those, a
#' `downgraded: <from> -> <to>` breakdown of the ones that actually reduced rank.
#' Returns a `$result` taxonomy table
#' (`id_col` + the winning `pct_id` + the 7 ranks -- REGATTA exchange format), a
#' per-ASV `$tracking` audit, and a step-level `$stats` summary; optionally
#' writes those as three CSVs.
#'
#' @param global_table Taxonomy table from the global-DB classifier:
#'   `id_col` + 7 lowercase rank columns + a percent-identity column.
#' @param local_table  Taxonomy table from the local-DB classifier: same shape.
#' @param id_col       ASV identifier column name, same in both inputs.
#' @param global_pct_id_col,local_pct_id_col Percent-identity column name
#'   in each input.
#' @param global_pct_id_scale,local_pct_id_scale One of `"auto"` (default),
#'   `"0-1"`, or `"0-100"`. `"auto"` detects 0-1 scale by checking if the
#'   input max is <= 1 and rescales to 0-100 with `signif(x, 3) * 100`.
#' @param Local_advantage TRUE (default): local wins ties at the
#'   `best_pctid` step. FALSE: global wins ties.
#' @param output_dir Directory path; default `NULL` writes nothing (the
#'   `result`/`tracking`/`stats` list is returned). Supply a directory to also
#'   write `<output_prefix>_taxonomy_table.csv`, `<output_prefix>_tracking.csv`,
#'   and (unless `write_summary = FALSE`) `<output_prefix>_summary.csv` there.
#' @param write_summary If `TRUE` (default), also write the `_summary.csv`.
#'   [run_regatta()] sets this `FALSE` so a run writes only one top-level
#'   summary (`regatta_summary.csv`).
#' @param output_prefix Filename prefix for the output CSVs. Default
#'   `"reconcile_global_local"`.
#' @param tracking_drop_pattern Regex matched against input column names;
#'   matching columns are dropped before they enter `$tracking`. Default
#'   strips obitools per-sample matrices and internal bookkeeping but
#'   preserves `BEST_MATCH_*`, `COUNT`, `TAXID`, `BEST_IDENTITY`,
#'   `NUC_SEQ`, and `SCIENTIFIC_NAME`.
#'
#' @return A list with three elements:
#'   * `result`: data.frame with `id_col`, then `pct_id` (the winning
#'     database's percent identity, kept right after the id for downstream
#'     filtering), then the 7 rank columns. Feeds straight into
#'     [reconcile_checklist()].
#'   * `tracking`: per-ASV decision record with columns grouped per field as
#'     `<field>_global` / `<field>_local` / `<field>_regatta` (the reconciled
#'     value) -- each of the 7 ranks and `pct_id` side by side -- then the
#'     decision columns (`best_pctid_winner`, `preferred_database`,
#'     `preferred_scientific_name`, `global_lca_to_local_triggered`) and any
#'     remaining classifier columns.
#'   * `stats`: (metric, count) data.frame.
#'
#' @importFrom utils write.csv
#' @export
reconcile_global_local <- function(global_table,
                                   local_table,
                                   id_col              = "ASV_id",
                                   global_pct_id_col   = "pct_id",
                                   local_pct_id_col    = "pct_id",
                                   global_pct_id_scale = c("auto", "0-1", "0-100"),
                                   local_pct_id_scale  = c("auto", "0-1", "0-100"),
                                   Local_advantage     = TRUE,
                                   output_dir          = NULL,
                                   output_prefix       = "reconcile_global_local",
                                   write_summary       = TRUE,
                                   tracking_drop_pattern =
                                     "^(MERGED_sample:|obiclean_|seq_rank|ID_STATUS|DEFINITION)") {
  # tracking_drop_pattern is a regex matched against column names in the
  # input tables BEFORE they get joined into $tracking. The default
  # strips obitools per-sample read-count matrices (MERGED_sample:*),
  # the obitools obiclean_* family, seq_rank, ID_STATUS, and
  # DEFINITION. These are not specific to the taxonomic assignment
  # decision REGATTA records, and they balloon the tracking CSV to
  # hundreds of irrelevant columns. BEST_MATCH_IDS / BEST_MATCH_TAXIDS
  # are NOT dropped by default because they record the obitools
  # winning accession and taxID -- useful audit information. Pass
  # tracking_drop_pattern = NULL (or "") to keep everything.

  ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  global_pct_id_scale <- match.arg(global_pct_id_scale)
  local_pct_id_scale  <- match.arg(local_pct_id_scale)

  # Validation
  if (!is.data.frame(global_table)) stop("global_table must be a data.frame.")
  if (!is.data.frame(local_table))  stop("local_table must be a data.frame.")
  mg <- setdiff(c(id_col, ranks, global_pct_id_col), names(global_table))
  if (length(mg) > 0) stop("global_table is missing required columns: ", paste(mg, collapse = ", "))
  ml <- setdiff(c(id_col, ranks, local_pct_id_col),  names(local_table))
  if (length(ml) > 0) stop("local_table is missing required columns: ",  paste(ml, collapse = ", "))

  # Drop obitools / per-sample bookkeeping from tracking. Essential
  # columns (id, ranks, pct_id) are preserved regardless of pattern.
  drop_extra <- function(df, essential) {
    if (is.null(tracking_drop_pattern) || !nzchar(tracking_drop_pattern)) return(df)
    candidates <- setdiff(names(df), essential)
    to_drop    <- candidates[grepl(tracking_drop_pattern, candidates)]
    df[, setdiff(names(df), to_drop), drop = FALSE]
  }
  global_table <- drop_extra(global_table, c(id_col, ranks, global_pct_id_col))
  local_table  <- drop_extra(local_table,  c(id_col, ranks, local_pct_id_col))

  # Rescale pct_id to a common 0-100 scale per input.
  rescale_pct <- function(x, scale_arg, label) {
    # Rescaling from 0-1 truncates to 3 significant figures first so the
    # precision matches what 0-100 classifiers typically report (e.g.
    # vsearch pctid 99.4). Matches the original Validate_local_assignments.R
    # treatment (signif(BEST_IDENTITY, 3) * 100) and prevents tie-break
    # flips from spurious sub-decimal differences across classifiers.
    x_num <- suppressWarnings(as.numeric(as.character(x)))
    if (scale_arg == "0-100") return(x_num)
    if (scale_arg == "0-1")   return(signif(x_num, 3) * 100)
    # auto: peek at the non-NA max
    mx <- suppressWarnings(max(x_num, na.rm = TRUE))
    if (is.finite(mx) && mx <= 1) {
      message(label, " pct_id detected as 0-1 scale (max = ",
              signif(mx, 3), "); rescaling to 0-100 (signif(x, 3) * 100).")
      return(signif(x_num, 3) * 100)
    }
    x_num
  }

  global_table[[global_pct_id_col]] <- rescale_pct(
    global_table[[global_pct_id_col]], global_pct_id_scale, "global")
  local_table[[local_pct_id_col]]   <- rescale_pct(
    local_table[[local_pct_id_col]],  local_pct_id_scale,  "local")

  sfx <- c("_global", "_local")
  joined <- merge(global_table, local_table,
                  by = id_col, all = TRUE, suffixes = sfx)

  g_rank_cols <- paste0(ranks, sfx[1])
  l_rank_cols <- paste0(ranks, sfx[2])

  # The two pct_id columns may have collided into <pct_id>_global /
  # <pct_id>_local if they share a name (the usual case) -- or if they
  # have distinct names, may not be suffixed at all. Resolve both.
  g_pct_col <- if (paste0(global_pct_id_col, sfx[1]) %in% names(joined)) {
    paste0(global_pct_id_col, sfx[1])
  } else global_pct_id_col
  l_pct_col <- if (paste0(local_pct_id_col, sfx[2]) %in% names(joined)) {
    paste0(local_pct_id_col, sfx[2])
  } else local_pct_id_col

  g_pct <- ifelse(is.na(joined[[g_pct_col]]), 0, joined[[g_pct_col]])
  l_pct <- ifelse(is.na(joined[[l_pct_col]]), 0, joined[[l_pct_col]])

  g_has_any <- !apply(is.na(joined[, g_rank_cols, drop = FALSE]), 1, all)
  l_has_any <- !apply(is.na(joined[, l_rank_cols, drop = FALSE]), 1, all)
  both_unassigned <- !g_has_any & !l_has_any

  # Step: best_pctid
  if (Local_advantage) {
    local_wins <- l_pct >= g_pct
  } else {
    local_wins <- l_pct > g_pct
  }
  best_pctid_winner <- ifelse(local_wins, "local", "global")

  preferred_lineage <- as.data.frame(
    matrix(NA_character_, nrow = nrow(joined), ncol = length(ranks),
           dimnames = list(NULL, ranks)),
    stringsAsFactors = FALSE
  )
  for (p in seq_along(ranks)) {
    preferred_lineage[[ranks[p]]] <- ifelse(
      local_wins, joined[[l_rank_cols[p]]], joined[[g_rank_cols[p]]]
    )
  }
  preferred_pctid <- ifelse(local_wins, l_pct, g_pct)
  database <- best_pctid_winner

  # Step: global_lca_to_local -- fires where global won best_pctid AND
  # local has any assignment to LCA against.
  global_won            <- database == "global"
  global_lca_to_local_triggered <- global_won & l_has_any

  if (any(global_lca_to_local_triggered)) {
    lca_rank <- rep(NA_character_, nrow(joined))
    for (p in seq(length(ranks), 1, by = -1)) {
      g <- joined[[g_rank_cols[p]]]
      l <- joined[[l_rank_cols[p]]]
      matched <- global_lca_to_local_triggered &
                 !is.na(g) & !is.na(l) & g == l & is.na(lca_rank)
      lca_rank[matched] <- ranks[p]
    }

    p_idx <- match(lca_rank, ranks)
    for (p in seq_along(ranks)) {
      fill <- global_lca_to_local_triggered & !is.na(p_idx) & p <= p_idx
      preferred_lineage[[ranks[p]]][fill]  <- joined[[g_rank_cols[p]]][fill]
      blank <- global_lca_to_local_triggered & (is.na(p_idx) | p > p_idx)
      preferred_lineage[[ranks[p]]][blank] <- NA_character_
    }

    database[global_lca_to_local_triggered] <- "global_lca_to_local"
  }

  # An ASV whose preferred lineage came out entirely empty is not assigned --
  # either both DBs were unassigned, or a global_lca_to_local row where global
  # and local share no rank at all (empty LCA, e.g. a cross-domain
  # disagreement). Mark it unassigned so the database label, pct_id, and the
  # "assigned" count stay consistent with the (empty) lineage.
  empty_lineage <- rowSums(!is.na(as.matrix(preferred_lineage[, ranks]))) == 0
  database[empty_lineage]                      <- NA_character_
  preferred_pctid[empty_lineage]               <- NA_real_
  global_lca_to_local_triggered[empty_lineage] <- FALSE

  # Downgrade breakdown for the global_lca_to_local step: of the triggered rows,
  # which actually reduced rank, and from which rank to which (global's deepest
  # rank -> the LCA rank). The rows where local merely agreed to global's depth
  # are triggered but contribute no downgrade pair.
  .depth_of <- function(M) vapply(seq_len(nrow(M)), function(i) {
    nz <- which(!is.na(M[i, ])); if (length(nz)) max(nz) else 0L
  }, integer(1))
  g_depth <- .depth_of(as.matrix(joined[, g_rank_cols, drop = FALSE]))
  p_depth <- .depth_of(as.matrix(preferred_lineage[, ranks, drop = FALSE]))
  glc_dn  <- global_lca_to_local_triggered & p_depth > 0 & p_depth < g_depth
  glc_pairs <- if (any(glc_dn))
    paste0("downgraded: ", ranks[g_depth[glc_dn]], " -> ", ranks[p_depth[glc_dn]]) else character(0)

  # Lowest non-NA rank in the preferred lineage
  preferred_scientific_name <- vapply(seq_len(nrow(joined)), function(i) {
    vals <- unlist(preferred_lineage[i, ranks])
    non_na <- !is.na(vals)
    if (!any(non_na)) return(NA_character_)
    vals[max(which(non_na))]
  }, character(1))

  # --- Build $result: ASV_id + the winning pct_id + 7 ranks ---
  # The winning database's percent identity travels with the winning taxa (right
  # after the id), so downstream percent-ID filtering and reconcile_checklist()
  # still have it.
  result <- data.frame(placeholder = joined[[id_col]], stringsAsFactors = FALSE)
  names(result)[1] <- id_col
  result$pct_id <- preferred_pctid
  for (r in ranks) result[[r]] <- preferred_lineage[[r]]
  rownames(result) <- NULL

  # --- Build $tracking: per-ASV audit, columns grouped per field as
  # <field>_global / <field>_local / <field>_regatta (the reconciled value),
  # so each rank and the pct_id read global / local / regatta side by side ---
  tracking <- joined
  # canonical pct_id names (the merge leaves these unsuffixed if the two inputs
  # used distinct pct_id column names)
  if (!("pct_id_global" %in% names(tracking))) tracking$pct_id_global <- joined[[g_pct_col]]
  if (!("pct_id_local"  %in% names(tracking))) tracking$pct_id_local  <- joined[[l_pct_col]]
  tracking$pct_id_regatta <- preferred_pctid
  for (r in ranks) tracking[[paste0(r, "_regatta")]] <- preferred_lineage[[r]]
  # decision / audit columns
  tracking$best_pctid_winner             <- best_pctid_winner
  tracking$preferred_database            <- database
  tracking$preferred_scientific_name     <- preferred_scientific_name
  tracking$global_lca_to_local_triggered <- global_lca_to_local_triggered

  # Order: id, then each rank's global/local/regatta triple, then the pct_id
  # triple, then the decision columns, then any remaining classifier columns.
  rank_triples <- as.vector(rbind(paste0(ranks, "_global"),
                                  paste0(ranks, "_local"),
                                  paste0(ranks, "_regatta")))
  front <- c(id_col, rank_triples,
             "pct_id_global", "pct_id_local", "pct_id_regatta",
             "best_pctid_winner", "preferred_database",
             "preferred_scientific_name", "global_lca_to_local_triggered")
  front <- front[front %in% names(tracking)]
  tracking <- tracking[, c(front, setdiff(names(tracking), front)), drop = FALSE]
  rownames(tracking) <- NULL

  # --- Build $stats ---
  n_total <- nrow(joined)
  stats <- data.frame(
    metric = c("total ASVs",
               "assigned ASVs",
               "both DBs unassigned",
               "best_pctid: local won",
               "best_pctid: global won",
               "global_lca_to_local triggered",
               "Final: database = local",
               "Final: database = global",
               "Final: database = global_lca_to_local"),
    count  = c(n_total,
               sum(!is.na(database)),
               sum(both_unassigned),
               sum(local_wins & !both_unassigned),
               sum(!local_wins & !both_unassigned),
               sum(global_lca_to_local_triggered),
               sum(database == "local",                na.rm = TRUE),
               sum(database == "global",               na.rm = TRUE),
               sum(database == "global_lca_to_local",  na.rm = TRUE)),
    stringsAsFactors = FALSE
  )

  # This step's transition headline: of the assigned ASVs, how many kept the
  # winning DB's full specificity vs were downgraded by the global_lca_to_local
  # LCA. Shares row labels with reconcile_checklist()'s transition so the report
  # can show both steps' counts side by side.
  n_downgraded <- length(glc_pairs)
  n_assigned   <- sum(!is.na(database))
  stats <- rbind(stats, data.frame(
    metric = c("ASVs unchanged (specificity kept)",
               "ASVs downgraded (specificity reduced)"),
    count  = c(n_assigned - n_downgraded, n_downgraded),
    stringsAsFactors = FALSE))

  # Of the global_lca_to_local rows that actually downgraded, the rank-pair
  # breakdown (global's call rank -> the LCA rank). Sums to the downgraded count
  # above; the rest were triggered but agreed to global's depth (no downgrade).
  if (length(glc_pairs)) {
    pt <- sort(table(glc_pairs), decreasing = TRUE)
    stats <- rbind(stats, data.frame(metric = names(pt),
                                     count  = as.integer(pt),
                                     stringsAsFactors = FALSE))
  }

  # Per-rank specificity counts: how far down each ASV's preferred
  # lineage actually goes. Useful when this function is run standalone.
  stats <- rbind(
    stats,
    data.frame(
      metric = c("ID'ed to kingdom only",
                 "ID'ed to phylum only",
                 "ID'ed to class only",
                 "ID'ed to order only",
                 "ID'ed to family only",
                 "ID'ed to genus only",
                 "ID'ed to species"),
      count  = c(sum(!is.na(result$domain)  & is.na(result$phylum)),
                 sum(!is.na(result$phylum)  & is.na(result$class)),
                 sum(!is.na(result$class)   & is.na(result$order)),
                 sum(!is.na(result$order)   & is.na(result$family)),
                 sum(!is.na(result$family)  & is.na(result$genus)),
                 sum(!is.na(result$genus)   & is.na(result$species)),
                 sum(!is.na(result$species))),
      stringsAsFactors = FALSE
    )
  )

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(result,   file.path(output_dir, paste0(output_prefix, "_taxonomy_table.csv")), row.names = FALSE)
    utils::write.csv(tracking, file.path(output_dir, paste0(output_prefix, "_tracking.csv")),       row.names = FALSE)
    if (isTRUE(write_summary))
      utils::write.csv(stats,  file.path(output_dir, paste0(output_prefix, "_summary.csv")),        row.names = FALSE)
    message("Wrote ", if (isTRUE(write_summary)) "3" else "2", " CSVs to ", normalizePath(output_dir))
  }

  list(result = result, tracking = tracking, stats = stats)
}
