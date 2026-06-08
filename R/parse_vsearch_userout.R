# parse_vsearch_userout.R
# Eldridge Wisely

# Parse a vsearch --userout file into a taxonomy table ready for the
# REGATTA pipeline. Useout files have one row per ASV-hit (top-N hits
# per ASV) with this column layout (no header):
#
#   col 1: ASV identifier
#   col 2: sequence;tax=<sintax string>
#   col 3: percent identity (0-100 scale)
#   col 4: alignment length
#   col 5: mismatches
#   col 6: gap opens
#
# This helper handles the common preprocessing pattern that the
# original Validate_local_assignments.R did inline:
#
#   - Take the first hit per ASV (matches the original code's
#     `left_join(..., multiple = "any")` semantics; the top-N hits in
#     a userout typically share the highest identity, so first-hit is
#     a reasonable canonical pick).
#   - Strip the "<sequence>;tax=" prefix off column 2 to isolate the
#     SINTAX taxonomy string.
#   - Parse the SINTAX string into the 7 lowercase rank columns via
#     parse_sintax().
#   - Carry the percent identity through as `pct_id`.
#
# Output: a data.frame with ASV_id, the 7 rank columns
# (domain..species), and pct_id. Drop-in as the `local_table` input to
# reconcile_global_local() or as a standalone taxonomy table to
# reconcile_checklist().

# Inputs:
#   path        Path to the vsearch userout file.
#   id_col      Column name to use for the ASV identifier in the
#               output (default "ASV_id").
#   first_hit_only
#               TRUE (default) -- take only the first row per ASV.
#               FALSE returns every hit, suffixed with `_hit_rank`
#               for the row's position within the per-ASV block.

#' Parse a vsearch `--userout` file into a REGATTA taxonomy table
#'
#' Reads a vsearch `--userout` file, takes the first hit per ASV,
#' strips the sequence prefix off column 2, parses the SINTAX taxonomy
#' via [parse_sintax()], and returns a clean table with `ASV_id` + the
#' 7 rank columns + `pct_id`.
#'
#' **Use [parse_vsearch_results()] in preference to this function**
#' when you also have the vsearch LCA file. The userout column-2
#' taxonomy is the *best hit's* taxonomy per ASV -- which can be
#' species-level even when vsearch's top-N hits disagree at species.
#' The LCA file carries the more conservative consensus across the
#' top-N. This function exists for the rarer case where only the
#' userout file is available.
#'
#' @param path Path to the vsearch userout file.
#' @param id_col Column name to use for the ASV identifier in the output
#'   (default `"ASV_id"`).
#' @param first_hit_only If TRUE (default), take only the first row per
#'   ASV. If FALSE, returns every hit.
#'
#' @return A data.frame: `id_col` + 7 rank columns + `pct_id`.
#'
#' @importFrom utils read.delim
#' @export
parse_vsearch_userout <- function(path,
                                  id_col         = "ASV_id",
                                  first_hit_only = TRUE) {
  message("Note: parse_vsearch_userout() returns the BEST-HIT taxonomy ",
          "per ASV, which can be more specific than the LCA across the ",
          "top-N hits. If you also have the vsearch LCA file, call ",
          "parse_vsearch_results(lca_path, userout_path) instead for the ",
          "more conservative consensus taxonomy.")
  if (!file.exists(path)) {
    stop("vsearch userout file not found at ", path)
  }

  raw <- utils::read.delim(path, header = FALSE, sep = "\t",
                           stringsAsFactors = FALSE,
                           col.names = c("ASV_id_raw", "seqtax",
                                         "pct_id_raw", "alnlen",
                                         "mism", "opens"))

  if (first_hit_only) {
    raw <- raw[!duplicated(raw$ASV_id_raw), , drop = FALSE]
  }

  # Strip the "<sequence>;tax=" prefix off column 2.
  sintax <- sub(".*;tax=", "", raw$seqtax)

  ranks_df <- parse_sintax(sintax)

  out <- data.frame(
    placeholder = raw$ASV_id_raw,
    stringsAsFactors = FALSE
  )
  names(out)[1] <- id_col
  out <- cbind(out, ranks_df,
               pct_id = suppressWarnings(as.numeric(raw$pct_id_raw)),
               stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}
