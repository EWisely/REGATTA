# resolve_taxids.R
# Eldridge Wisely

# Resolve a vector of NCBI taxIDs into a 7-rank taxonomy data.frame using
# taxonomizr's local SQL DB. NA taxIDs come back as all-NA rows. Output
# columns: domain, phylum, class, order, family, genus, species.

# This is the preprocessor for shape-2 inputs (e.g. obitools3 output where
# each ASV has a TAXID but no parsed ranks). Bind the result to the
# original ID column yourself before passing to regatta_checklist_lca().

resolve_taxids <- function(taxids, sql_path = "accessionTaxa.sql") {
  if (!requireNamespace("taxonomizr", quietly = TRUE)) {
    stop("Package 'taxonomizr' is required.")
  }
  if (!file.exists(sql_path)) {
    stop("SQL DB not found at ", sql_path,
         ". Build it with taxonomizr::prepareDatabase() or pass a path to an existing DB.")
  }

  ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

  # desiredTaxa = ranks ensures the modern "domain" rank is populated;
  # the default "superkingdom" is NA in current NCBI dumps.
  taxa_mat <- taxonomizr::getTaxonomy(taxids, sql_path, desiredTaxa = ranks)
  taxa_df  <- as.data.frame(taxa_mat, stringsAsFactors = FALSE)
  rownames(taxa_df) <- NULL
  taxa_df[, ranks, drop = FALSE]
}
