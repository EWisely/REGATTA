# regatta_compare_assignments.R
# Eldridge Wisely

# Optional comparison feature. Given two taxonomy tables that classify
# the same ASVs via different reference databases (e.g. obitools+global
# DB vs. vsearch+local DB), produce a per-ASV side-by-side comparison
# plus a "merged" taxonomy table containing the lowest common ancestor
# between the two assignments. Comparison is symmetric — no notion of
# which database is "preferred".

# Inputs are full taxonomy tables in the shape regatta_checklist_lca()
# expects: an ASV-id column plus the 7 lowercase rank columns. All other
# columns (pct_id, count, sequence, etc.) are preserved through the
# join with input-specific suffixes so the user can audit both sides.

# Output: a single data.frame with one row per ASV in the UNION of both
# inputs. Columns:
#   <id_col>            the ASV identifier
#   <all input cols from A, suffixed _A>
#   <all input cols from B, suffixed _B>
#   agreement_rank      lowest taxonomic rank where A and B agree
#                       (one of species/genus/.../domain, or NA)
#   agreement_category  one of:
#                         "agree_at_<rank>"      A and B reach the same
#                                                 taxon at that rank
#                         "disagree"             both have assignments
#                                                 but no common ancestor
#                                                 in the input
#                         "only_<label_A>"       only A has an assignment
#                         "only_<label_B>"       only B has an assignment
#                         "both_unassigned"      both A and B all-NA
#   merged_<rank>       LCA of A and B: agreed values from the agreement
#                       rank up to domain; NA at all lower ranks. If no
#                       agreement, all NA. This is feedable directly
#                       into regatta_checklist_lca() as a taxonomy_table.

# Pure base R; no taxonomizr dependency.

regatta_compare_assignments <- function(table_A,
                                        table_B,
                                        id_col   = "ASV_id",
                                        label_A  = "A",
                                        label_B  = "B") {
  ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

  for (nm in c("table_A", "table_B")) {
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
  if (label_A == label_B) {
    stop("label_A and label_B must differ so column suffixes are distinct.")
  }

  sfx <- c(paste0("_", label_A), paste0("_", label_B))

  joined <- merge(table_A, table_B, by = id_col, all = TRUE, suffixes = sfx)

  a_rank_cols <- paste0(ranks, sfx[1])
  b_rank_cols <- paste0(ranks, sfx[2])

  # Lowest rank at which A and B agree. Walk species -> domain so a more
  # specific agreement wins. Vectorized across all rows for performance.
  agreement_rank <- rep(NA_character_, nrow(joined))
  for (p in seq(length(ranks), 1, by = -1)) {
    a_vals <- joined[[a_rank_cols[p]]]
    b_vals <- joined[[b_rank_cols[p]]]
    matched <- !is.na(a_vals) & !is.na(b_vals) &
               a_vals == b_vals & is.na(agreement_rank)
    agreement_rank[matched] <- ranks[p]
  }

  # Build merged_* columns: for each row, ranks at or ABOVE the
  # agreement rank get the agreed value; ranks BELOW stay NA. (NCBI
  # lineages are internally consistent, so ranks above the agreement
  # also necessarily match between A and B.)
  merged <- as.data.frame(
    matrix(NA_character_, nrow = nrow(joined), ncol = length(ranks),
           dimnames = list(NULL, paste0("merged_", ranks))),
    stringsAsFactors = FALSE
  )
  p_idx <- match(agreement_rank, ranks)
  for (p in seq_along(ranks)) {
    fill <- !is.na(p_idx) & p <= p_idx
    merged[[paste0("merged_", ranks[p])]][fill] <- joined[[a_rank_cols[p]]][fill]
  }

  a_has_any <- !apply(is.na(joined[, a_rank_cols, drop = FALSE]), 1, all)
  b_has_any <- !apply(is.na(joined[, b_rank_cols, drop = FALSE]), 1, all)

  category <- rep(NA_character_, nrow(joined))
  category[!a_has_any & !b_has_any] <- "both_unassigned"
  category[ a_has_any & !b_has_any] <- paste0("only_", label_A)
  category[!a_has_any &  b_has_any] <- paste0("only_", label_B)
  both <- a_has_any & b_has_any
  category[both & is.na(agreement_rank)] <- "disagree"
  category[both & !is.na(agreement_rank)] <-
    paste0("agree_at_", agreement_rank[both & !is.na(agreement_rank)])

  joined$agreement_rank     <- agreement_rank
  joined$agreement_category <- category
  out <- cbind(joined, merged)
  rownames(out) <- NULL
  out
}
