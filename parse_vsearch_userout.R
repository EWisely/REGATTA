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
#               TRUE (default) — take only the first row per ASV.
#               FALSE returns every hit, suffixed with `_hit_rank`
#               for the row's position within the per-ASV block.

parse_vsearch_userout <- function(path,
                                  id_col         = "ASV_id",
                                  first_hit_only = TRUE) {
  if (!file.exists(path)) {
    stop("vsearch userout file not found at ", path)
  }
  if (!exists("parse_sintax", mode = "function")) {
    source("parse_sintax.R")
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
