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

resolve_taxids <- function(input,
                           taxid_col     = NULL,
                           sql_path      = "accessionTaxa.sql",
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
