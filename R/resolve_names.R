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

# Name lookup is synonym-aware: matches against NCBI scientific names
# AND recorded synonyms (via name_to_taxid in regatta_helpers.R), so
# inputs like "Lagenorhynchus obliquidens" resolve correctly to
# Sagmatias obliquidens after the recent cetacean revisions. Common
# names and other categories are excluded by default. Override via
# accept_types to broaden or restrict the policy.

# Cleanup applied: strips sp., spp., cf., aff., Gen., indet., and quote
# characters; collapses whitespace. Word-boundary anchored so e.g.
# "Sebastes" is not corrupted to "ebastes".
#' Strip common taxonomic abbreviations from names (internal)
#' @keywords internal
#' @noRd
clean_taxon_names <- function(x) {
  x <- gsub('"', "", x, fixed = TRUE)
  x <- gsub("\\b(spp?|cf|aff|Gen|indet)\\.\\s*", "", x, perl = TRUE)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

#' Convert scientific names (mixed ranks) to a 7-rank taxonomy table
#'
#' Resolves a column of scientific names (e.g. Kraken2 / BestTaxon-style
#' outputs where the rank varies row-to-row) into a full 7-rank lineage.
#' Strips common abbreviations (`sp.`, `spp.`, `cf.`, `aff.`, `Gen.`,
#' `indet.`, quotes) before NCBI lookup. Synonym-aware: matches the NCBI
#' `scientific name` and `synonym` types by default so older nomenclature
#' is normalized to the current canonical name.
#'
#' @param input A character vector or a data.frame.
#' @param name_col Required when `input` is a data.frame.
#' @param sql_path Path to the local `accessionTaxa.sql` taxonomizr DB.
#' @param output_prefix If non-NULL, writes a CSV with this prefix.
#' @param output_dir Directory to write the CSV into.
#' @param clean If TRUE (default), strip junk strings before lookup.
#' @param accept_types NCBI name types to accept; default is scientific
#'   name plus recorded synonyms.
#'
#' @return The input augmented with the 7 rank columns + an
#'   `name_match_type` column (`"scientific name"` / `"synonym"` / NA).
#'
#' @importFrom utils write.csv
#' @export
resolve_names <- function(input,
                          name_col      = NULL,
                          sql_path      = "accessionTaxa.sql",
                          output_prefix = NULL,
                          output_dir    = ".",
                          clean         = TRUE,
                          accept_types  = c("scientific name", "synonym")) {
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

  lookup <- name_to_taxid(names_vec, sql_path, accept_types = accept_types)
  taxa_mat <- taxonomizr::getTaxonomy(lookup$taxID, sql_path, desiredTaxa = ranks)
  taxa_df  <- as.data.frame(taxa_mat, stringsAsFactors = FALSE)
  rownames(taxa_df) <- NULL

  # Don't clobber existing rank columns silently -- if any are already
  # present, append with a "resolved_" prefix so the user can compare.
  rank_collisions <- intersect(ranks, names(df))
  if (length(rank_collisions) > 0) {
    names(taxa_df) <- paste0("resolved_", ranks)
  }

  # Attach name_match_type so downstream code (and the summary report)
  # can tell whether each row resolved via canonical scientific name or
  # via a synonym (i.e. the name was normalized to current NCBI
  # taxonomy). NA means the name did not resolve.
  out <- cbind(df,
               name_match_type = lookup$match_type,
               taxa_df[, intersect(c(ranks, paste0("resolved_", ranks)),
                                   names(taxa_df)), drop = FALSE])
  rownames(out) <- NULL

  n_unres <- sum(is.na(lookup$taxID))
  if (n_unres > 0) {
    message(n_unres, " of ", length(lookup$taxID),
            " name(s) did not resolve in NCBI -- see NA rows.")
  }
  n_syn <- sum(lookup$match_type == "synonym", na.rm = TRUE)
  if (n_syn > 0) {
    message(n_syn, " of ", length(lookup$taxID),
            " name(s) resolved via synonym (auto-updated to current canonical NCBI taxonomy).")
  }

  if (!is.null(output_prefix)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    out_path <- file.path(output_dir, paste0(output_prefix, "_full_tax_table.csv"))
    utils::write.csv(out, out_path, row.names = FALSE)
    message("Wrote ", out_path)
  }

  out
}
