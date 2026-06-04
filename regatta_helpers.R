# regatta_helpers.R
# Eldridge Wisely

# Internal helpers shared by resolve_names() and taxonomize_checklist().
# Currently houses the synonym-aware name -> taxID lookup.

# name_to_taxid() resolves a vector of scientific names to NCBI taxIDs
# using the local taxonomizr SQL DB. Unlike taxonomizr::getId(), which
# matches only the canonical scientific name when onlyScientific = TRUE
# and ALL name categories when FALSE, this helper restricts matching to
# a configurable set of name types — by default the canonical scientific
# name plus recorded synonyms. Common names, BLAST names, acronyms, and
# similar categories are excluded to keep the lookup explicit and
# defensible for downstream methods reporting.

# Returns an integer vector of taxIDs aligned to the input. Unmatched
# names get NA. When a name maps to taxIDs across multiple types (e.g.
# it is both a scientific name for one taxon and a synonym for another),
# the scientific-name match is preferred and a warning is emitted.

name_to_taxid <- function(taxa,
                          sql_path,
                          accept_types = c("scientific name", "synonym")) {
  if (!requireNamespace("RSQLite", quietly = TRUE)) {
    stop("Package 'RSQLite' is required.")
  }
  if (!file.exists(sql_path)) {
    stop("SQL DB not found at ", sql_path)
  }

  taxa <- as.character(taxa)
  unique_taxa <- unique(taxa[!is.na(taxa) & nzchar(taxa)])
  if (length(unique_taxa) == 0) {
    return(rep(NA_integer_, length(taxa)))
  }

  # ATTACH a temp DB with the query names — same pattern taxonomizr uses
  # internally. Scales to thousands of names without hitting SQLite's
  # IN-clause length limit.
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  tmp_db <- RSQLite::dbConnect(RSQLite::SQLite(), tmp)
  on.exit(RSQLite::dbDisconnect(tmp_db), add = TRUE)
  RSQLite::dbWriteTable(tmp_db, "query",
                        data.frame(name = unique_taxa, stringsAsFactors = FALSE),
                        overwrite = TRUE)

  db <- RSQLite::dbConnect(RSQLite::SQLite(), sql_path)
  on.exit(RSQLite::dbDisconnect(db), add = TRUE)
  RSQLite::dbExecute(db, sprintf("ATTACH '%s' AS tmp", tmp))

  # Build the type filter as a parenthesized string-literal list.
  # Escape single quotes defensively even though our defaults don't have any.
  escaped <- gsub("'", "''", accept_types, fixed = TRUE)
  type_list <- paste0("'", escaped, "'", collapse = ", ")

  # LEFT OUTER JOIN with the type predicate IN THE JOIN CLAUSE (not WHERE)
  # so unmatched query names still appear as rows with NULL taxID.
  q <- sprintf(
    "SELECT q.name AS query_name,
            n.id   AS taxid,
            n.type AS match_type,
            n.scientific AS is_scientific
     FROM tmp.query AS q
     LEFT OUTER JOIN names AS n
       ON q.name = n.name AND n.type IN (%s)",
    type_list
  )
  res <- RSQLite::dbGetQuery(db, q)

  # For each query name, pick the best match: prefer scientific over
  # synonym; warn on ambiguity (multiple distinct taxIDs).
  pick_best <- function(rows) {
    if (all(is.na(rows$taxid))) return(NA_integer_)
    rows <- rows[!is.na(rows$taxid), , drop = FALSE]
    if (any(rows$is_scientific == 1, na.rm = TRUE)) {
      return(as.integer(rows$taxid[which(rows$is_scientific == 1)[1]]))
    }
    as.integer(rows$taxid[1])
  }
  is_ambiguous <- function(rows) {
    rows <- rows[!is.na(rows$taxid), , drop = FALSE]
    length(unique(rows$taxid)) > 1
  }

  per_name   <- split(res, res$query_name)
  best       <- vapply(per_name, pick_best, integer(1))
  ambiguous  <- vapply(per_name, is_ambiguous, logical(1))

  if (any(ambiguous)) {
    warning("Multiple taxIDs found for: ",
            paste(names(ambiguous)[ambiguous], collapse = ", "),
            ". Picked scientific-name match where available, else first.",
            call. = FALSE)
  }

  unname(best[match(taxa, names(best))])
}
