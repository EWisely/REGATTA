# reconcile_checklist.R
# Eldridge Wisely

# reconcile_checklist(): Reconcile a taxonomy table against a regional
# species checklist by walking each row's lineage from species up
# toward domain, finding the lowest rank present in the checklist,
# and NA-ing every rank below that match. This is the REGATTA
# "off-target downgrade" step (originally Pass 3 in
# Validate_local_assignments.R) — it preserves specificity where the
# regional checklist supports it and downgrades where it doesn't,
# without using percent-identity heuristics.
#
# By default, output_dir = "reconcile_checklist_out" and the three
# return-list elements are written to disk there as CSVs named
#   <output_prefix>_taxonomy_table.csv  (the corrected $result, strict 8 cols)
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
#     extras — extras are ignored for the LCA itself).
#
# Output: a list with three elements.
#
#   $result    The corrected per-ASV taxonomy table. EXACTLY 8
#              columns: id_col + the 7 lowercase rank columns. Same
#              shape as reconcile_global_local()$result so REGATTA
#              functions chain cleanly. Drop-in to phyloseq
#              tax_table() or to a MetabaR MOTU table after joining
#              read counts back.
#
#   $tracking  Per-ASV before-vs-after audit: before_<rank> and
#              after_<rank> at every rank, regatta_match_rank (the
#              rank at which the input matched the checklist, or NA
#              for fully off-target), and any input metadata columns
#              passed through. One row per ASV in the input.
#
#   $stats     Match-rank distribution and per-rank specificity
#              counts of the corrected output, plus the change in
#              number of ASVs assigned before vs. after.

reconcile_checklist <- function(taxonomy_table,
                                checklist,
                                id_col        = "ASV_id",
                                output_dir    = "reconcile_checklist_out",
                                output_prefix = "reconcile_checklist",
                                prior_dir     = "reconcile_global_local_out",
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
  missing_c <- setdiff(ranks, names(checklist))
  if (length(missing_c) > 0) {
    stop("checklist is missing required rank columns: ",
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

  checklist_sets <- lapply(ranks, function(r) unique(na.omit(checklist[[r]])))
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
  # get every rank NA'd (fully off-target).
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

  # --- Build $stats ---
  kin_only <- sum(!is.na(corrected$domain)  & is.na(corrected$phylum))
  phy_only <- sum(!is.na(corrected$phylum)  & is.na(corrected$class))
  cla_only <- sum(!is.na(corrected$class)   & is.na(corrected$order))
  ord_only <- sum(!is.na(corrected$order)   & is.na(corrected$family))
  fam_only <- sum(!is.na(corrected$family)  & is.na(corrected$genus))
  gen_only <- sum(!is.na(corrected$genus)   & is.na(corrected$species))
  sp_ct    <- sum(!is.na(corrected$species))

  n_before <- sum(!is.na(before$domain))
  n_after  <- sum(!is.na(corrected$domain))

  stats <- data.frame(
    metric = c("total ASVs",
               "assigned before checklist-LCA",
               "assigned after checklist-LCA",
               "change in number of ASVs assigned",
               "matched at species",
               "matched at genus",
               "matched at family",
               "matched at order",
               "matched at class",
               "matched at phylum",
               "matched at domain",
               "not matched (off-target)",
               "ID'ed to kingdom only",
               "ID'ed to phylum only",
               "ID'ed to class only",
               "ID'ed to order only",
               "ID'ed to family only",
               "ID'ed to genus only",
               "ID'ed to species"),
    count  = c(nrow(corrected),
               n_before,
               n_after,
               n_after - n_before,
               sum(match_rank == "species", na.rm = TRUE),
               sum(match_rank == "genus",   na.rm = TRUE),
               sum(match_rank == "family",  na.rm = TRUE),
               sum(match_rank == "order",   na.rm = TRUE),
               sum(match_rank == "class",   na.rm = TRUE),
               sum(match_rank == "phylum",  na.rm = TRUE),
               sum(match_rank == "domain",  na.rm = TRUE),
               sum(is.na(match_rank)),
               kin_only, phy_only, cla_only, ord_only, fam_only, gen_only, sp_ct),
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
    utils::write.csv(summary_to_write,
                     file.path(output_dir, paste0(output_prefix, "_summary.csv")),
                     row.names = FALSE)
    message("Wrote 3 CSVs to ", normalizePath(output_dir))
  }

  list(result = result, tracking = tracking, stats = stats)
}
