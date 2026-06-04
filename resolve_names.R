# resolve_names.R
# Eldridge Wisely

# Resolve a column of scientific names (mixed taxonomic ranks) into a
# full taxonomy table in the shape regatta_checklist_lca() expects.
# Each name is looked up via taxonomizr::getId() and the lineage filled
# via getTaxonomy(); the matched rank is populated and ranks below it
# come back as NA. Names that fail to resolve become all-NA rows.

# Accepts either:
#   - a character vector of names (returns a df with input_name + 7 ranks)
#   - a data.frame plus name_col (returns input cols + 7 rank cols)
# If output_prefix is non-NULL, also writes
#   <output_dir>/<output_prefix>_full_tax_table.csv
# for use in downstream phyloseq / MetabaR pipelines.

# This is the preprocessor for Kraken2- or BestTaxon-style outputs where
# each row has a single name at whatever rank the classifier could
# reach. Common abbreviations (sp., spp., cf., aff., Gen., indet.) and
# embedded quotes are stripped before lookup.

# Cleanup applied: strips sp., spp., cf., aff., Gen., indet., and quote
# characters; collapses whitespace. Word-boundary anchored so e.g.
# "Sebastes" is not corrupted to "ebastes".
clean_taxon_names <- function(x) {
  x <- gsub('"', "", x, fixed = TRUE)
  x <- gsub("\\b(spp?|cf|aff|Gen|indet)\\.\\s*", "", x, perl = TRUE)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

resolve_names <- function(input,
                          name_col      = NULL,
                          sql_path      = "accessionTaxa.sql",
                          output_prefix = NULL,
                          output_dir    = ".",
                          clean         = TRUE) {
  if (!requireNamespace("taxonomizr", quietly = TRUE)) {
    stop("Package 'taxonomizr' is required.")
  }
  if (!file.exists(sql_path)) {
    stop("SQL DB not found at ", sql_path,
         ". Build it with taxonomizr::prepareDatabase() or pass a path to an existing DB.")
  }

  ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

  # Normalize input to (df, names_vec)
  if (is.character(input)) {
    df <- data.frame(input_name = input, stringsAsFactors = FALSE)
    name_col <- "input_name"
  } else if (is.data.frame(input)) {
    if (is.null(name_col)) {
      stop("`name_col` must be supplied when `input` is a data.frame.")
    }
    if (!name_col %in% names(input)) {
      stop("Column '", name_col, "' not found in input data.frame.")
    }
    df <- input
  } else {
    stop("`input` must be a character vector or a data.frame.")
  }

  names_vec <- as.character(df[[name_col]])
  if (clean) names_vec <- clean_taxon_names(names_vec)

  taxids <- taxonomizr::getId(names_vec, sql_path)
  taxa_mat <- taxonomizr::getTaxonomy(taxids, sql_path, desiredTaxa = ranks)
  taxa_df  <- as.data.frame(taxa_mat, stringsAsFactors = FALSE)
  rownames(taxa_df) <- NULL

  # Don't clobber existing rank columns silently — if any are already
  # present, append with a "resolved_" prefix so the user can compare.
  rank_collisions <- intersect(ranks, names(df))
  if (length(rank_collisions) > 0) {
    names(taxa_df) <- paste0("resolved_", ranks)
  }

  out <- cbind(df, taxa_df[, intersect(c(ranks, paste0("resolved_", ranks)),
                                       names(taxa_df)), drop = FALSE])
  rownames(out) <- NULL

  n_unres <- sum(is.na(taxids))
  if (n_unres > 0) {
    message(n_unres, " of ", length(taxids),
            " name(s) did not resolve in NCBI — see NA rows.")
  }

  if (!is.null(output_prefix)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    out_path <- file.path(output_dir, paste0(output_prefix, "_full_tax_table.csv"))
    utils::write.csv(out, out_path, row.names = FALSE)
    message("Wrote ", out_path)
  }

  out
}
