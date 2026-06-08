# summarize_regatta.R
# Eldridge Wisely
#
# Summary-table layout designed by Ella Crotty.

# Compare the input(s) and output(s) of a REGATTA run and produce
# Ella's 24-row stats summary -- one column per stage, capturing how
# taxonomic specificity and diversity shifted as the data flowed
# through reconciliation and the regional checklist filter.
#
# Inputs are optional and combinable. The function fills whatever
# rows it can from what you supply:
#
#   reconciled     The output list from reconcile_global_local() --
#                  uses $result for the "reconciled" stage column and
#                  $tracking for the source-breakdown rows 4-7
#                  (count of local/global assignment preferred).
#   post_checklist The output list from reconcile_checklist() -- uses
#                  $result for the "regatta" stage column and the
#                  before/after change-in-assignment count for row 8.
#   global_input   The raw global-DB classifier output taxonomy table
#                  (e.g. obitools-derived). Used to populate the
#                  "global" stage column.
#   local_input    The raw local-DB classifier output taxonomy table
#                  (e.g. vsearch+SINTAX-derived). Used to populate
#                  the "local" stage column.
#
# Each supplied input becomes its own column. Stage column names are
# fixed: global, local, reconciled, regatta.
#
# Row groups (Ella's layout, plus the checklist-membership pair):
#   1-3   counts (total, assigned, % assigned)
#   4-7   source breakdown (local-preferred and global-preferred from
#         the best_pctid step -- only populated when `reconciled` is
#         supplied)
#   8     change in number of ASVs assigned through the checklist
#         step (only populated when `post_checklist` is supplied)
#   9-10  data -> checklist membership -- % of a stage's calls on the regional
#         checklist, per ASV and per distinct taxon (only when `checklist`
#         is supplied; global < 100%, local & regatta = 100%)
#   11    checklist -> data recovery -- % of checklist species detected in the
#         stage's data (survey completeness; only when `checklist` is supplied)
#   12-18 ID'ed-to-<rank>-only specificity counts
#   19-24 diversity counts (distinct phyla, classes, ..., species)

#' Summarize a REGATTA run as a 24-row stats table
#'
#' Compares inputs and outputs across stages and produces Ella's summary:
#' counts, source breakdown, checklist membership, ID'ed-to-rank specificity,
#' and diversity. Each supplied input becomes its own column in the output.
#' Stage columns are fixed: `global`, `local`, `reconciled`, `regatta`.
#'
#' @param reconciled Output list of [reconcile_global_local()].
#' @param post_checklist Output list of [reconcile_checklist()].
#' @param global_input Raw global-DB classifier output taxonomy table.
#' @param local_input Raw local-DB classifier output taxonomy table.
#' @param checklist The taxonomized regional checklist (a data.frame with the
#'   7 rank columns -- e.g. `build_regional_checklist()$for_LCA`). When supplied,
#'   the checklist-comparison rows are filled per stage: the two "percent ... on
#'   regional checklist" rows (how much of the *data* is on the checklist;
#'   global typically < 100%, local and regatta = 100%) and "percent of
#'   checklist species detected" (the reverse -- how much of the *checklist* the
#'   data recovered). Otherwise those rows are `NA`. [run_regatta()] passes the
#'   run's checklist automatically.
#'
#' @return A 24-row data.frame with one column per supplied stage.
#'
#' @export
summarize_regatta <- function(reconciled     = NULL,
                              post_checklist = NULL,
                              global_input   = NULL,
                              local_input    = NULL,
                              checklist      = NULL) {
  # Argument order is tuned for the common call patterns:
  #   summarize_regatta(post)                 single-DB workflow
  #   summarize_regatta(rec, post)            two-DB workflow
  #   summarize_regatta(rec, post, g, l)      full audit including raw inputs
  # All four arguments are still accepted by name in any order.
  ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  # "Assigned" = any non-NA rank, so the count is right even for inputs that
  # carry no `domain` column (e.g. a pre-resolved phylum..species table).
  n_assigned_any <- function(t) sum(rowSums(!is.na(t[, ranks, drop = FALSE])) > 0)

  # Checklist membership: the fraction of a stage's calls that are on the
  # regional checklist (global < 100%; local and regatta = 100% by
  # construction). Needs the taxonomized checklist; NA when it isn't supplied
  # or lacks rank columns. Computed two ways: per ASV, and per distinct taxon.
  cl_sets <- NULL
  if (!is.null(checklist) && is.data.frame(checklist) &&
      all(ranks %in% names(checklist))) {
    cl_sets <- lapply(ranks, function(r) unique(stats::na.omit(checklist[[r]])))
    names(cl_sets) <- ranks
  }
  pct_on_checklist <- function(t) {
    if (is.null(cl_sets)) return(c(asv = NA_real_, taxon = NA_real_))
    M    <- as.matrix(t[, ranks, drop = FALSE])
    low  <- apply(!is.na(M), 1, function(r) if (any(r)) max(which(r)) else NA_integer_)
    keep <- !is.na(low)
    if (!any(keep)) return(c(asv = NA_real_, taxon = NA_real_))
    low <- low[keep]; Mk <- M[keep, , drop = FALSE]
    vals <- vapply(seq_along(low), function(i) Mk[i, low[i]], character(1))
    on   <- vapply(seq_along(low), function(i) vals[i] %in% cl_sets[[low[i]]], logical(1))
    dup  <- !duplicated(paste(low, vals, sep = "\r"))   # distinct (rank, taxon)
    c(asv   = 100 * sum(on)      / length(on),
      taxon = 100 * sum(on[dup]) / sum(dup))
  }
  # The other direction: of the checklist's species, how many did this stage's
  # data detect? (Survey recovery / completeness against the regional list.)
  recovery_species <- function(t) {
    if (is.null(cl_sets) || !length(cl_sets$species)) return(NA_real_)
    data_sp <- unique(stats::na.omit(t$species))
    100 * length(intersect(data_sp, cl_sets$species)) / length(cl_sets$species)
  }

  # Collect stage tables in fixed order
  stages <- list()
  if (!is.null(global_input))   stages[["global"]]     <- global_input
  if (!is.null(local_input))    stages[["local"]]      <- local_input
  if (!is.null(reconciled)) {
    if (!"result" %in% names(reconciled)) {
      stop("`reconciled` must be the output of reconcile_global_local() (a list with $result).")
    }
    stages[["reconciled"]] <- reconciled$result
  }
  if (!is.null(post_checklist)) {
    if (!"result" %in% names(post_checklist)) {
      stop("`post_checklist` must be the output of regatta_checklist_lca() (a list with $result).")
    }
    stages[["regatta"]] <- post_checklist$result
  }

  if (length(stages) == 0) {
    stop("Supply at least one of global_input, local_input, reconciled, post_checklist.")
  }
  for (nm in names(stages)) {
    missing_cols <- setdiff(ranks, names(stages[[nm]]))
    if (length(missing_cols) > 0) {
      stop("Stage '", nm, "' is missing required rank columns: ",
           paste(missing_cols, collapse = ", "))
    }
  }

  per_stage <- function(t) {
    n_total      <- nrow(t)
    n_assigned   <- n_assigned_any(t)
    pct_assigned <- if (n_total > 0) 100 * n_assigned / n_total else NA_real_

    kin_only <- sum(!is.na(t$domain)  & is.na(t$phylum))
    phy_only <- sum(!is.na(t$phylum)  & is.na(t$class))
    cla_only <- sum(!is.na(t$class)   & is.na(t$order))
    ord_only <- sum(!is.na(t$order)   & is.na(t$family))
    fam_only <- sum(!is.na(t$family)  & is.na(t$genus))
    gen_only <- sum(!is.na(t$genus)   & is.na(t$species))
    sp_count <- sum(!is.na(t$species))

    n_distinct <- function(x) length(unique(x[!is.na(x)]))
    div <- vapply(ranks[-1], function(r) n_distinct(t[[r]]), integer(1))

    pc <- pct_on_checklist(t)

    c(n_total, n_assigned, pct_assigned,
      NA, NA, NA, NA, NA,                                # rows 4-8 fill below
      pc[["asv"]], pc[["taxon"]],                        # rows 9-10 data -> checklist
      recovery_species(t),                              # row 11 checklist -> data
      kin_only, phy_only, cla_only, ord_only, fam_only, gen_only, sp_count,
      div["phylum"], div["class"], div["order"],
      div["family"], div["genus"], div["species"])
  }

  cols <- lapply(stages, per_stage)
  out  <- as.data.frame(cols, stringsAsFactors = FALSE)
  names(out) <- names(stages)

  # Source-breakdown rows (4-7) from $reconciled$tracking.
  # Matches original Validate_local_assignments.R semantics:
  #  - "local preferred"  = best_pctid local won AND global made any
  #                         assignment (i.e., local OVERRODE a
  #                         non-empty global call)
  #  - "global preferred" = final preferred_database is "global"
  #                         (best_pctid global won AND local empty;
  #                         global_lca_to_local cases are separate)
  #  - denominator        = ASVs where global made any assignment
  if (!is.null(reconciled) && "tracking" %in% names(reconciled)) {
    tr <- reconciled$tracking
    g_rank_cols <- paste0(ranks, "_global")
    have_g <- if (all(g_rank_cols %in% names(tr))) {
      !apply(is.na(tr[, g_rank_cols, drop = FALSE]), 1, all)
    } else rep(TRUE, nrow(tr))

    # Both counts use the best_pctid step's winner so the local and
    # global counts are symmetric and match the original
    # Validate_local_assignments.R stats. Global wins that later get
    # downgraded by global_lca_to_local still count here as "global
    # preferred at the best_pctid step."
    n_local  <- sum(tr$best_pctid_winner == "local"  & have_g)
    n_global <- sum(tr$best_pctid_winner == "global" & have_g)
    denom    <- sum(have_g)

    for (j in seq_along(out)) {
      out[[j]][4] <- n_local
      out[[j]][5] <- if (denom > 0) 100 * n_local  / denom else NA_real_
      out[[j]][6] <- n_global
      out[[j]][7] <- if (denom > 0) 100 * n_global / denom else NA_real_
    }
  }

  # Row 8: change in number of ASVs assigned by adding local-DB
  # information to global-only. Computed as nrow(assigned in
  # reconciled) - nrow(assigned in global_input). This is the
  # "local-DB rescue count": how many additional ASVs got an
  # assignment because the local DB contributed something the global
  # DB missed. Matches original Validate_local_assignments.R semantics
  # (post_asg - obi_assigned).
  if (!is.null(reconciled) && !is.null(global_input)) {
    n_global_assigned     <- n_assigned_any(global_input)
    n_reconciled_assigned <- n_assigned_any(reconciled$result)
    for (j in seq_along(out)) out[[j]][8] <- n_reconciled_assigned - n_global_assigned
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
    "percent of ASVs on regional checklist",
    "percent of distinct taxa on regional checklist",
    "percent of checklist species detected",
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
