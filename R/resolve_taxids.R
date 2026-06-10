# resolve_taxids.R
# Eldridge Wisely

# Resolve a column of NCBI taxIDs into a full taxonomy table in the
# shape regatta_checklist_lca() expects. NA taxIDs become all-NA rows.

# Accepts either:
#   - a numeric/character vector of taxIDs (returns df with taxID + 7 ranks)
#   - a data.frame plus taxid_col (returns input cols + 7 rank cols)
# If output_prefix is non-NULL, also writes
#   <output_dir>/<output_prefix>_full_tax_table.csv
# for use in downstream phyloseq / MetabaR pipelines.

# This is the preprocessor for shape-2 inputs (e.g. obitools3 output
# where each ASV has a TAXID column).

#' Convert NCBI taxIDs to a 7-rank taxonomy table
#'
#' Resolves a column of NCBI taxIDs into a full 7-rank lineage using
#' `taxonomizr`. Accepts either a vector or a data.frame plus `taxid_col`.
#' NA taxIDs become all-NA rows. Optionally writes the augmented table to
#' disk as `<output_prefix>_full_tax_table.csv`.
#'
#' @param input A numeric/character vector of taxIDs, or a data.frame.
#' @param taxid_col Required when `input` is a data.frame.
#' @param sql_path Path to the local `accessionTaxa.sql` taxonomizr DB.
#'   Defaults to the persistent per-user cache shared across REGATTA
#'   (`tools::R_user_dir("REGATTA", "cache")`).
#' @param output_prefix If non-NULL, writes a CSV with this prefix.
#' @param output_dir Directory to write the CSV into.
#'
#' @return The input augmented with the 7 lowercase rank columns
#'   (`domain` through `species`); all original columns are preserved.
#'
#' @examples
#' \dontrun{
#' # Resolve NCBI taxIDs to a 7-rank lineage. Needs the local taxonomizr
#' # database (built on first use, or pass an existing one via sql_path).
#' resolve_taxids(c(8030, 8022))
#'
#' # From a data.frame column, keeping the other columns:
#' resolve_taxids(data.frame(ASV_id = c("ASV_1", "ASV_2"), TAXID = c(8030, 8022)),
#'                taxid_col = "TAXID")
#' }
#'
#' @importFrom utils write.csv
#' @export
resolve_taxids <- function(input,
                           taxid_col     = NULL,
                           sql_path      = .regatta_default_sql_path(),
                           output_prefix = NULL,
                           output_dir    = ".") {
  if (!requireNamespace("taxonomizr", quietly = TRUE)) {
    stop("Package 'taxonomizr' is required.")
  }
  if (!file.exists(sql_path)) {
    stop("SQL DB not found at ", sql_path,
         ". Build it with taxonomizr::prepareDatabase() or pass a path to an existing DB.")
  }

  ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

  if (is.numeric(input) || is.character(input)) {
    df <- data.frame(taxID = input, stringsAsFactors = FALSE)
    taxid_col <- "taxID"
  } else if (is.data.frame(input)) {
    if (is.null(taxid_col)) {
      stop("`taxid_col` must be supplied when `input` is a data.frame.")
    }
    if (!taxid_col %in% names(input)) {
      stop("Column '", taxid_col, "' not found in input data.frame.")
    }
    df <- input
  } else {
    stop("`input` must be a vector of taxIDs or a data.frame.")
  }

  taxids <- df[[taxid_col]]
  # desiredTaxa = ranks ensures the modern "domain" rank is populated;
  # the default "superkingdom" is NA in current NCBI dumps.
  taxa_mat <- taxonomizr::getTaxonomy(taxids, sql_path, desiredTaxa = ranks)
  taxa_df  <- as.data.frame(taxa_mat, stringsAsFactors = FALSE)
  rownames(taxa_df) <- NULL

  rank_collisions <- intersect(ranks, names(df))
  if (length(rank_collisions) > 0) {
    names(taxa_df) <- paste0("resolved_", ranks)
  }

  out <- cbind(df, taxa_df)
  rownames(out) <- NULL

  if (!is.null(output_prefix)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    out_path <- file.path(output_dir, paste0(output_prefix, "_full_tax_table.csv"))
    utils::write.csv(out, out_path, row.names = FALSE)
    message("Wrote ", out_path)
  }

  out
}
