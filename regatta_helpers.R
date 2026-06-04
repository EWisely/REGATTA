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

# Returns a data.frame aligned to the input with columns:
#   taxID      — integer NCBI taxID, NA for unmatched names
#   match_type — the NCBI name type that matched ("scientific name",
#                "synonym", or whatever else is in accept_types). NA
#                for unmatched names.
# Surfacing match_type lets downstream code distinguish "this name
# resolved via a synonym (the assignment got normalized to current
# canonical NCBI taxonomy)" from "this name is the current canonical
# name" — both for reporting (which the summary table will use to
# label the change category) and for filtering.

# When a name maps to taxIDs across multiple types (e.g. it is both a
# scientific name for one taxon and a synonym for another), the
# scientific-name match is preferred and a warning is emitted.

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
  na_result <- data.frame(
    taxID      = rep(NA_integer_,   length(taxa)),
    match_type = rep(NA_character_, length(taxa)),
    stringsAsFactors = FALSE
  )
  unique_taxa <- unique(taxa[!is.na(taxa) & nzchar(taxa)])
  if (length(unique_taxa) == 0) {
    return(na_result)
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
  # synonym; warn on ambiguity (multiple distinct taxIDs). Return both
  # the chosen taxID and the matched type so callers can label rows.
  pick_best <- function(rows) {
    if (all(is.na(rows$taxid))) {
      return(c(taxid = NA_integer_, match_type = NA_character_))
    }
    rows <- rows[!is.na(rows$taxid), , drop = FALSE]
    if (any(rows$is_scientific == 1, na.rm = TRUE)) {
      i <- which(rows$is_scientific == 1)[1]
    } else {
      i <- 1L
    }
    c(taxid = as.integer(rows$taxid[i]), match_type = rows$match_type[i])
  }
  is_ambiguous <- function(rows) {
    rows <- rows[!is.na(rows$taxid), , drop = FALSE]
    length(unique(rows$taxid)) > 1
  }

  per_name  <- split(res, res$query_name)
  picks     <- do.call(rbind, lapply(per_name, pick_best))
  picks_df  <- data.frame(
    query_name = rownames(picks),
    taxID      = as.integer(picks[, "taxid"]),
    match_type = unname(picks[, "match_type"]),
    stringsAsFactors = FALSE
  )
  ambiguous <- vapply(per_name, is_ambiguous, logical(1))
  if (any(ambiguous)) {
    warning("Multiple taxIDs found for: ",
            paste(names(ambiguous)[ambiguous], collapse = ", "),
            ". Picked scientific-name match where available, else first.",
            call. = FALSE)
  }

  # Align to original input order (preserves NA + duplicate handling)
  idx <- match(taxa, picks_df$query_name)
  out <- data.frame(
    taxID      = picks_df$taxID[idx],
    match_type = picks_df$match_type[idx],
    stringsAsFactors = FALSE
  )
  # taxa that were NA or empty in the input get NA taxID and NA match_type
  out$taxID[is.na(taxa)      | !nzchar(as.character(taxa))] <- NA_integer_
  out$match_type[is.na(taxa) | !nzchar(as.character(taxa))] <- NA_character_
  out
}
