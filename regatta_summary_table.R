# regatta_summary_table.R
# Eldridge Wisely
#
# Summary-table layout designed by Ella Crotty.

# Build a per-stage stats summary from one or more taxonomy tables.
# Output is a data.frame with exactly the 21 rows below (one column per
# stage in the order they appear in `stages`) so users can drop it into
# supplementary tables without renaming anything.
#
# Row groups:
#   1-3   counts (total, assigned, % assigned)
#   4-7   assignment-source breakdown (populated only when `comparison`
#         and the two labels are supplied; NA otherwise)
#   8     change in number of ASVs assigned (populated only when
#         `lca_result` is supplied; NA otherwise)
#   9-15  ID'ed-to-<rank>-only specificity counts
#   16-21 diversity counts (number of distinct phyla, classes, ...,
#         species)

# Inputs:
#   stages       Named list of taxonomy tables (each: ASV-id + 7 rank
#                columns; other columns ignored for stats). Names become
#                column headers in the output. Use whatever labels make
#                sense for the user's pipeline ("global", "local",
#                "post", or just "before"/"after").
#   comparison   Optional output of regatta_compare_assignments(). When
#                supplied with both labels, the source-breakdown rows
#                fill in (counts and percentages of only_<label> cases).
#   global_label, local_label
#                Character labels matching the `label_A`/`label_B` used
#                in the comparison call. Required together with
#                `comparison` to fill the source-breakdown rows.
#   lca_result   Optional output of regatta_checklist_lca(). When
#                supplied, row 8 is populated as nrow(corrected with
#                preferred_name) minus nrow(before with preferred_name).
#
# The "assignment-source breakdown" semantics differ slightly from the
# original Validate_local_assignments.R, which counted cases where one
# DB was *picked over* the other based on percent identity. In the
# redesigned symmetric comparison there is no winner-takes-all step;
# only_<label> counts the cases where one classifier had an assignment
# and the other did not — a clean, well-defined version of the same
# "which DB contributed this ASV" question.

regatta_summary_table <- function(stages,
                                  comparison   = NULL,
                                  global_label = NULL,
                                  local_label  = NULL,
                                  lca_result   = NULL) {
  ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

  if (!is.list(stages) || is.null(names(stages)) || any(!nzchar(names(stages)))) {
    stop("`stages` must be a named list of taxonomy tables.")
  }
  for (nm in names(stages)) {
    missing_cols <- setdiff(ranks, names(stages[[nm]]))
    if (length(missing_cols) > 0) {
      stop("Stage '", nm, "' is missing required rank columns: ",
           paste(missing_cols, collapse = ", "))
    }
  }

  per_stage <- function(t) {
    n_total <- nrow(t)
    n_assigned <- sum(!is.na(t$domain))
    pct_assigned <- if (n_total > 0) 100 * n_assigned / n_total else NA_real_

    # ID'ed-to-<rank>-only: row counts where that rank is the most
    # specific populated rank. Uses NA-pattern below the rank.
    kin_only <- sum(!is.na(t$domain)  & is.na(t$phylum))
    phy_only <- sum(!is.na(t$phylum)  & is.na(t$class))
    cla_only <- sum(!is.na(t$class)   & is.na(t$order))
    ord_only <- sum(!is.na(t$order)   & is.na(t$family))
    fam_only <- sum(!is.na(t$family)  & is.na(t$genus))
    gen_only <- sum(!is.na(t$genus)   & is.na(t$species))
    sp_count <- sum(!is.na(t$species))

    # Diversity: unique non-NA values per rank.
    n_distinct <- function(x) length(unique(x[!is.na(x)]))
    div <- vapply(ranks[-1], function(r) n_distinct(t[[r]]), integer(1))

    c(n_total, n_assigned, pct_assigned,
      NA, NA, NA, NA, NA,                                # rows 4-8 fill below
      kin_only, phy_only, cla_only, ord_only, fam_only, gen_only, sp_count,
      div["phylum"], div["class"], div["order"],
      div["family"], div["genus"], div["species"])
  }

  cols <- lapply(stages, per_stage)
  # Assemble: 21 rows × N stages
  out <- as.data.frame(cols, stringsAsFactors = FALSE)
  names(out) <- names(stages)

  # Source-breakdown rows (4-7) when comparison + both labels supplied
  if (!is.null(comparison) && !is.null(global_label) && !is.null(local_label)) {
    if (!"agreement_category" %in% names(comparison)) {
      stop("`comparison` must be the output of regatta_compare_assignments().")
    }
    n_local  <- sum(comparison$agreement_category == paste0("only_", local_label))
    n_global <- sum(comparison$agreement_category == paste0("only_", global_label))
    n_cmp    <- nrow(comparison)
    # Fill the source-breakdown rows for EVERY stage column with the same
    # comparison-derived numbers. They are properties of the comparison,
    # not of any single stage.
    for (j in seq_along(out)) {
      out[[j]][4] <- n_local
      out[[j]][5] <- if (n_cmp > 0) 100 * n_local / n_cmp else NA_real_
      out[[j]][6] <- n_global
      out[[j]][7] <- if (n_cmp > 0) 100 * n_global / n_cmp else NA_real_
    }
  }

  # Row 8: change in number of ASVs assigned, from lca_result.
  if (!is.null(lca_result)) {
    if (!all(c("before", "corrected") %in% names(lca_result))) {
      stop("`lca_result` must be the output of regatta_checklist_lca().")
    }
    # Use any rank: an ASV is "assigned" if any rank is non-NA. Use the
    # most permissive (domain) on both sides.
    n_before <- sum(!is.na(lca_result$before$domain))
    n_after  <- sum(!is.na(lca_result$corrected$domain))
    for (j in seq_along(out)) out[[j]][8] <- n_after - n_before
  }

  row_labels <- c(
    "total ASVs",
    "assigned ASVs",
    "Percentage of ASVs assigned",
    "count of local assignment preferred",
    "percent local assignments",
    "count of global assignment preferred",
    "percent global assignments",
    "change in number of ASVs assigned",
    "ID'ed to kingdom only",
    "ID'ed to phylum only",
    "ID'ed to class only",
    "ID'ed to order only",
    "ID'ed to family only",
    "ID'ed to genus only",
    "ID'ed to species",
    "Number of phyla",
    "Number of classes",
    "Number of orders",
    "Number of families",
    "Number of genera",
    "Number of species"
  )

  cbind(row_names = row_labels, out, stringsAsFactors = FALSE)
}
