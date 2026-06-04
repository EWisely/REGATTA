# reconcile_global_local.R
# Eldridge Wisely

# Optional reconciliation step. Given two taxonomy tables that classify
# the same set of ASVs against a GLOBAL reference database (e.g.
# obitools + EMBL/NCBI) and a LOCAL reference database (e.g.
# vsearch + a regionally-curated DB), produce:
#
#   $reconciled  A clean per-ASV taxonomy table (ASV_id + 7 ranks +
#                reconciliation_status) suitable as a phyloseq
#                tax_table() input or a MetabaR MOTU table.
#
#   $summary     A best_ID_combined-style audit table preserving every
#                column from both inputs (suffixed _global / _local)
#                plus reconciliation_rank, reconciliation_status, and
#                the reconciled_<rank> lineage. This is the per-ASV
#                decision record — origin of each input assignment,
#                the final reconciled call, and all the data the
#                reconciliation considered.
#
# Reconciliation policy (per ASV):
#
#   - Both DBs agree at rank R          → use the agreed lineage from
#                                          R upward; NA below R.
#   - Both DBs assigned but disagree at every rank
#                                       → reconciled lineage is all NA
#                                          (flagged "disagree" — likely
#                                          contamination or a real
#                                          cross-domain mismatch).
#   - Only the global DB assigned       → use the global lineage as-is.
#   - Only the local DB assigned        → use the local lineage as-is.
#   - Both unassigned                   → all NA.
#
# This differs from the original Validate_local_assignments.R three-pass
# %ID-based cascade: there is no winner-takes-all by percent identity,
# and no implicit preference for one DB over the other. Where the two
# DBs agree, the agreement wins; where only one has a call, that call
# stands; where they truly conflict, the row is flagged for review.
# Percent identities (if present in either input) are preserved through
# the summary table so the user can post-hoc filter by them.

# Pure base R; no taxonomizr dependency.

reconcile_global_local <- function(global_table,
                                   local_table,
                                   id_col = "ASV_id") {
  ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

  for (nm in c("global_table", "local_table")) {
    t <- get(nm)
    if (!is.data.frame(t)) {
      stop(nm, " must be a data.frame.")
    }
    missing_cols <- setdiff(c(id_col, ranks), names(t))
    if (length(missing_cols) > 0) {
      stop(nm, " is missing required columns: ",
           paste(missing_cols, collapse = ", "))
    }
  }

  sfx <- c("_global", "_local")
  joined <- merge(global_table, local_table,
                  by = id_col, all = TRUE, suffixes = sfx)

  g_rank_cols <- paste0(ranks, sfx[1])
  l_rank_cols <- paste0(ranks, sfx[2])

  # Lowest rank at which global and local agree. Walk species -> domain
  # so a more specific agreement wins. Vectorized across rows.
  reconciliation_rank <- rep(NA_character_, nrow(joined))
  for (p in seq(length(ranks), 1, by = -1)) {
    g <- joined[[g_rank_cols[p]]]
    l <- joined[[l_rank_cols[p]]]
    matched <- !is.na(g) & !is.na(l) & g == l & is.na(reconciliation_rank)
    reconciliation_rank[matched] <- ranks[p]
  }

  g_has_any <- !apply(is.na(joined[, g_rank_cols, drop = FALSE]), 1, all)
  l_has_any <- !apply(is.na(joined[, l_rank_cols, drop = FALSE]), 1, all)

  status <- rep(NA_character_, nrow(joined))
  status[!g_has_any & !l_has_any] <- "both_unassigned"
  status[ g_has_any & !l_has_any] <- "only_global"
  status[!g_has_any &  l_has_any] <- "only_local"
  both <- g_has_any & l_has_any
  status[both & is.na(reconciliation_rank)] <- "disagree"
  status[both & !is.na(reconciliation_rank)] <-
    paste0("agree_at_", reconciliation_rank[both & !is.na(reconciliation_rank)])

  # Build the reconciled lineage per row, following the policy above.
  reconciled_ranks <- as.data.frame(
    matrix(NA_character_, nrow = nrow(joined), ncol = length(ranks),
           dimnames = list(NULL, paste0("reconciled_", ranks))),
    stringsAsFactors = FALSE
  )

  # Case 1: both agree at some rank — fill the LCA from that rank up.
  p_idx <- match(reconciliation_rank, ranks)
  for (p in seq_along(ranks)) {
    fill <- !is.na(p_idx) & p <= p_idx
    reconciled_ranks[[paste0("reconciled_", ranks[p])]][fill] <-
      joined[[g_rank_cols[p]]][fill]
  }

  # Case 2: only the global DB assigned — use its full lineage.
  only_g <- status == "only_global"
  for (p in seq_along(ranks)) {
    reconciled_ranks[[paste0("reconciled_", ranks[p])]][only_g] <-
      joined[[g_rank_cols[p]]][only_g]
  }

  # Case 3: only the local DB assigned — use its full lineage.
  only_l <- status == "only_local"
  for (p in seq_along(ranks)) {
    reconciled_ranks[[paste0("reconciled_", ranks[p])]][only_l] <-
      joined[[l_rank_cols[p]]][only_l]
  }

  # Cases 4 & 5 (disagree, both_unassigned): leave reconciled all NA.

  # Audit table (full): everything that went into the decision
  summary <- cbind(joined,
                   reconciliation_rank   = reconciliation_rank,
                   reconciliation_status = status,
                   reconciled_ranks)
  rownames(summary) <- NULL

  # Reconciled table (clean): ASV-id + 7 ranks + status — ready for
  # phyloseq tax_table() or MetabaR MOTU input.
  reconciled <- data.frame(
    placeholder = joined[[id_col]],
    domain  = reconciled_ranks$reconciled_domain,
    phylum  = reconciled_ranks$reconciled_phylum,
    class   = reconciled_ranks$reconciled_class,
    order   = reconciled_ranks$reconciled_order,
    family  = reconciled_ranks$reconciled_family,
    genus   = reconciled_ranks$reconciled_genus,
    species = reconciled_ranks$reconciled_species,
    reconciliation_status = status,
    stringsAsFactors = FALSE
  )
  names(reconciled)[1] <- id_col
  rownames(reconciled) <- NULL

  list(reconciled = reconciled, summary = summary)
}
