# parse_sintax.R
# Eldridge Wisely

# Parse vsearch SINTAX taxonomy strings into a full taxonomy table in
# the shape regatta_checklist_lca() expects. SINTAX strings look like:
#   d:Eukaryota,p:Chordata,c:Actinopteri,o:Scombriformes,f:Scombridae,g:Auxis,s:Auxis_thazard
# Truncated strings (e.g. only down to genus) are accepted; missing
# ranks come back as NA. Empty strings and NA inputs return all-NA rows.
# Species values have underscores converted to spaces.

# Accepts either:
#   - a character vector of SINTAX strings (returns df with sintax + 7 ranks)
#   - a data.frame plus sintax_col (returns input cols + 7 rank cols)
# If output_prefix is non-NULL, also writes
#   <output_dir>/<output_prefix>_full_tax_table.csv
# for use in downstream phyloseq / MetabaR pipelines.

#' Parse a vsearch SINTAX taxonomy string into a 7-rank table
#'
#' Splits SINTAX taxonomy strings (e.g.
#' `d:Eukaryota,p:Chordata,...,s:Foo_bar`) into the 7 lowercase rank
#' columns. Accepts either a character vector or a data.frame plus
#' `sintax_col`. Truncated strings (e.g. only down to genus) are
#' accepted; missing ranks come back as NA. Species values have
#' underscores converted to spaces.
#'
#' @param input A character vector or a data.frame.
#' @param sintax_col Required when `input` is a data.frame.
#' @param output_prefix If non-NULL, writes a CSV with this prefix.
#' @param output_dir Directory to write the CSV into.
#'
#' @return The input augmented with the 7 rank columns.
#'
#' @examples
#' # Split SINTAX taxonomy strings into the 7 rank columns; the second string
#' # is truncated at genus, so its deeper ranks come back NA.
#' parse_sintax(c(
#'   "d:Eukaryota,p:Chordata,c:Actinopteri,o:Perciformes,f:Lutjanidae,g:Lutjanus,s:Lutjanus_kasmira",
#'   "d:Eukaryota,p:Chordata,c:Actinopteri,o:Mugiliformes,f:Mugilidae,g:Mugil"))
#'
#' @importFrom utils write.csv
#' @export
parse_sintax <- function(input,
                         sintax_col    = NULL,
                         output_prefix = NULL,
                         output_dir    = ".") {
  ranks    <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  prefixes <- c("d",      "p",      "c",     "o",     "f",      "g",     "s")

  if (is.character(input)) {
    df <- data.frame(sintax = input, stringsAsFactors = FALSE)
    sintax_col <- "sintax"
  } else if (is.data.frame(input)) {
    if (is.null(sintax_col)) {
      stop("`sintax_col` must be supplied when `input` is a data.frame.")
    }
    if (!sintax_col %in% names(input)) {
      stop("Column '", sintax_col, "' not found in input data.frame.")
    }
    df <- input
  } else {
    stop("`input` must be a character vector or a data.frame.")
  }

  sintax <- as.character(df[[sintax_col]])

  extract_rank <- function(s, prefix) {
    pat <- paste0("(?:^|,)", prefix, ":([^,]*)")
    out <- rep(NA_character_, length(s))
    ok <- !is.na(s) & nzchar(s)
    if (!any(ok)) return(out)
    m <- regmatches(s[ok], regexec(pat, s[ok]))
    out[ok] <- vapply(m, function(x) if (length(x) >= 2) x[2] else NA_character_,
                      character(1))
    out
  }

  taxa_df <- data.frame(
    stats::setNames(lapply(prefixes, extract_rank, s = sintax), ranks),
    stringsAsFactors = FALSE
  )
  for (r in ranks) {
    empty <- !is.na(taxa_df[[r]]) & !nzchar(taxa_df[[r]])
    taxa_df[[r]][empty] <- NA_character_
  }
  taxa_df$species <- gsub("_", " ", taxa_df$species, fixed = TRUE)

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
