# regatta_checklist_lca.R
# Eldridge Wisely

# Reconcile a taxonomy table against a regional species checklist by walking
# each row's lineage from species up toward domain, finding the lowest rank
# present in the checklist, and NA-ing every rank below that match. Preserves
# all non-rank columns (pct_id, count, sequence, etc.) untouched. Returns a
# list with before/corrected/changes for snapshotting and downstream summaries.

# Inputs are required to be in "shape 1" — already-resolved 7-rank columns.
# Use parse_sintax() or a taxID/accession/name resolver upstream if needed.
# This function has no taxonomizr dependency by design.

regatta_checklist_lca <- function(taxonomy_table,
                                  checklist,
                                  id_col = "ASV_id") {
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

  checklist_sets <- lapply(ranks, function(r) unique(na.omit(checklist[[r]])))
  names(checklist_sets) <- ranks

  before <- taxonomy_table
  corrected <- taxonomy_table

  # Find the lowest rank at which each row matches the checklist.
  # Walk from species (most specific, index 7) up to domain (index 1).
  match_idx <- rep(NA_integer_, nrow(taxonomy_table))
  for (p in seq(length(ranks), 1, by = -1)) {
    vals <- taxonomy_table[[ranks[p]]]
    hit <- !is.na(vals) & vals %in% checklist_sets[[ranks[p]]] & is.na(match_idx)
    match_idx[hit] <- p
  }

  # NA out ranks below each row's matched rank.
  # Rows with no match anywhere get every rank NA'd.
  for (p in seq_along(ranks)) {
    blank <- is.na(match_idx) | p > match_idx
    corrected[[ranks[p]]][blank] <- NA
  }

  match_rank <- ifelse(is.na(match_idx), NA_character_, ranks[match_idx])
  corrected$regatta_match_rank <- match_rank

  # Build the changes table: one row per ASV whose taxonomy was modified.
  na_safe_diff <- function(a, b) {
    xor(is.na(a), is.na(b)) | (!is.na(a) & !is.na(b) & a != b)
  }
  any_changed <- Reduce(
    `|`,
    lapply(ranks, function(r) na_safe_diff(before[[r]], corrected[[r]]))
  )

  changes <- before[any_changed, id_col, drop = FALSE]
  for (r in ranks) {
    changes[[paste0("before_", r)]] <- before[[r]][any_changed]
    changes[[paste0("after_", r)]] <- corrected[[r]][any_changed]
  }
  changes$match_rank <- match_rank[any_changed]

  list(before = before, corrected = corrected, changes = changes)
}
