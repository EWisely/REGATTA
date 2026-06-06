# parse_vsearch_results.R
# Eldridge Wisely

# Canonical vsearch preprocessor. Joins the LCA file (one row per ASV,
# SINTAX-format LCA taxonomy across the top-N hits) with the userout
# file (multiple rows per ASV, each with pct_id and the per-hit
# taxonomy). The LCA file is the source of truth for taxonomy; the
# userout file contributes pct_id only. Matches the original
# Validate_local_assignments.R left_join(lca, userout, multiple = "any")
# semantics: first hit per ASV from userout for pct_id.

#' Parse a vsearch LCA + userout pair into a REGATTA taxonomy table
#'
#' Joins the conservative LCA call from vsearch's `--lcaout` file with
#' the percent identity from vsearch's `--userout` file (first hit per
#' ASV) and returns a clean per-ASV table: `ASV_id` + 7 lowercase rank
#' columns + `pct_id`. Drop-in as `global_table` or `local_table` for
#' [reconcile_global_local()] or as the input to [reconcile_checklist()].
#'
#' Use this in preference to [parse_vsearch_userout()] when you have
#' both files. The userout's column 2 carries the *best hit's* taxonomy
#' (potentially species-level even when the top-N hits disagree); the
#' LCA file carries the *consensus* across the top-N (often more
#' conservative). Calling reconciliation on the LCA taxonomy matches
#' the original code semantics and is methodologically more honest.
#'
#' @param lca_path     Path to the vsearch LCA file (2 cols: ID, sintax).
#' @param userout_path Path to the vsearch userout file (6 cols: ID,
#'   sequence;tax, pct_id, alnlen, mism, opens).
#' @param id_col Column name to use for the ASV identifier in the
#'   output (default `"ASV_id"`).
#'
#' @return A data.frame with `id_col` + 7 rank columns + `pct_id`.
#'   ASVs in the LCA file with no userout hit get `pct_id = NA`.
#'
#' @importFrom utils read.delim
#' @export
parse_vsearch_results <- function(lca_path,
                                  userout_path,
                                  id_col = "ASV_id") {
  if (!file.exists(lca_path)) {
    stop("vsearch LCA file not found at ", lca_path)
  }
  if (!file.exists(userout_path)) {
    stop("vsearch userout file not found at ", userout_path)
  }

  lca <- utils::read.delim(lca_path, header = FALSE, sep = "\t",
                           stringsAsFactors = FALSE,
                           col.names = c("ASV_id_raw", "sintax"))
  ranks_df <- parse_sintax(lca$sintax)

  user <- utils::read.delim(userout_path, header = FALSE, sep = "\t",
                            stringsAsFactors = FALSE,
                            col.names = c("ASV_id_raw", "seqtax",
                                          "pct_id_raw", "alnlen",
                                          "mism", "opens"))
  user_first <- user[!duplicated(user$ASV_id_raw), c("ASV_id_raw", "pct_id_raw"),
                     drop = FALSE]
  names(user_first) <- c("ASV_id_raw", "pct_id")
  user_first$pct_id <- suppressWarnings(as.numeric(user_first$pct_id))

  out <- data.frame(placeholder = lca$ASV_id_raw, stringsAsFactors = FALSE)
  names(out)[1] <- id_col
  out <- cbind(out, ranks_df, stringsAsFactors = FALSE)
  out$pct_id <- user_first$pct_id[match(lca$ASV_id_raw, user_first$ASV_id_raw)]
  rownames(out) <- NULL
  out
}
