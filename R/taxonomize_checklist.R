# taxonomize_checklist.R
# Eldridge Wisely

# Resolve a list of scientific names into a 7-rank taxonomy table using
# taxonomizr's local NCBI taxonomy SQL DB. Cleans common junk strings
# (Gen., indet., sp., cf., quotes) before lookup. Returns one row per
# unique input name with `input_name`, `resolution_status`, and the 7
# lowercase rank columns (domain -> species). Unresolved names stay in
# the output with NA ranks and resolution_status = "unresolved".

# Name lookup is synonym-aware: matches against NCBI scientific names
# AND recorded synonyms (via name_to_taxid in regatta_helpers.R), so
# regional checklists that contain older or synonymous names still
# resolve correctly. Common names and other categories are excluded by
# default. Override via accept_types to broaden or restrict the policy.

# Intended to run once per region per taxonomic group. The resulting
# data.frame is the checklist input to regatta_checklist_lca().

#' Resolve a regional species checklist to a 7-rank NCBI taxonomy table
#'
#' Cleans junk strings (`sp.`, `spp.`, `cf.`, etc.), looks each name up
#' via `taxonomizr` (synonym-aware), and returns a one-row-per-species
#' table with the 7 lowercase rank columns plus `input_name`,
#' `resolution_status`, `name_match_type`, and `taxID`. Run once per
#' region per taxonomic group; the result is the `checklist` input to
#' [reconcile_checklist()].
#'
#' @param input A path to a checklist file (one-column or with a
#'   `Species` column), a character vector of names, or a data.frame
#'   with a `Species` / `scientific_name` column.
#' @param sql_path Path to the local `accessionTaxa.sql` taxonomizr DB.
#'   Defaults to the persistent per-user cache shared across REGATTA
#'   (`tools::R_user_dir("REGATTA", "cache")`).
#' @param prepare_db Build the SQL DB if missing (names+nodes only, via
#'   `taxonomizr::prepareDatabase(getAccessions = FALSE)` -- a few hundred MB).
#' @param accept_types NCBI name types to accept; default is scientific
#'   name plus recorded synonyms.
#'
#' @return A data.frame, one row per unique input name: `input_name`,
#'   `resolution_status` (`"resolved"`/`"unresolved"`), `name_match_type`
#'   (`"scientific name"`/`"synonym"`/`NA`), `taxID`, and the 7 lowercase rank
#'   columns (`domain` through `species`, `NA` where unresolved).
#'
#' @importFrom utils read.delim
#' @export
taxonomize_checklist <- function(input,
                                 sql_path = .regatta_default_sql_path(),
                                 prepare_db = !file.exists(sql_path),
                                 accept_types = c("scientific name", "synonym")) {
  if (!requireNamespace("taxonomizr", quietly = TRUE)) {
    stop("Package 'taxonomizr' is required.")
  }

  # Normalize input to a character vector of names
  names_in <- if (is.character(input) && length(input) == 1 && file.exists(input)) {
    first_line <- readLines(input, n = 1, warn = FALSE)
    if (grepl("Species", first_line, ignore.case = TRUE)) {
      df <- utils::read.delim(input, stringsAsFactors = FALSE, fileEncoding = "latin1")
      species_col <- grep("^species$", names(df), ignore.case = TRUE, value = TRUE)[1]
      if (is.na(species_col)) stop("File has a header but no 'Species' column.")
      df[[species_col]]
    } else {
      readLines(input, warn = FALSE)
    }
  } else if (is.character(input)) {
    input
  } else if (is.data.frame(input)) {
    species_col <- grep("^(species|scientific_name)$",
                       names(input), ignore.case = TRUE, value = TRUE)[1]
    if (is.na(species_col)) {
      stop("data.frame input must have a 'Species' or 'scientific_name' column.")
    }
    input[[species_col]]
  } else {
    stop("`input` must be a file path, a character vector, or a data.frame with a Species column.")
  }

  names_in <- names_in[!is.na(names_in) & nzchar(trimws(names_in))]

  # Strip junk strings (matches existing cleaning in Validate_local_assignments.R)
  cleaned <- names_in
  cleaned <- gsub("Gen\\. ",   "", cleaned)
  cleaned <- gsub("indet\\. ", "", cleaned)
  cleaned <- gsub("sp\\. ",    "", cleaned)
  cleaned <- gsub("cf\\. ",    "", cleaned)
  cleaned <- gsub("\"",        "", cleaned)
  cleaned <- trimws(cleaned)
  cleaned <- unique(cleaned[nzchar(cleaned)])

  if (prepare_db) {
    if (file.exists(sql_path)) {
      message("`prepare_db = TRUE` but ", sql_path,
              " already exists -- skipping rebuild.")
    } else {
      # names+nodes only -- all getId()/getTaxonomy() need -- and stamps the
      # build date (shared with build_regional_checklist()).
      .regatta_build_taxonomy_db(sql_path)
    }
  } else if (!file.exists(sql_path)) {
    stop("SQL DB not found at ", sql_path,
         ". Pass `prepare_db = TRUE` to build, or provide the path to an existing DB.")
  }

  ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

  lookup <- name_to_taxid(cleaned, sql_path, accept_types = accept_types)
  # Request domain explicitly -- getTaxonomy's default asks for "superkingdom",
  # which is NA in current NCBI dumps where the top rank is named "domain".
  taxa_mat <- taxonomizr::getTaxonomy(lookup$taxID, sql_path, desiredTaxa = ranks)
  taxa_df <- as.data.frame(taxa_mat, stringsAsFactors = FALSE)

  # name_match_type carries the NCBI type each name matched on
  # ("scientific name" or "synonym" by default; NA for unresolved).
  # resolution_status is kept as a coarse pass/fail for back-compat.
  out <- cbind(
    data.frame(
      input_name        = cleaned,
      resolution_status = ifelse(is.na(lookup$taxID), "unresolved", "resolved"),
      name_match_type   = lookup$match_type,
      taxID             = lookup$taxID,
      stringsAsFactors  = FALSE
    ),
    taxa_df[, ranks, drop = FALSE]
  )
  rownames(out) <- NULL

  n_unresolved <- sum(out$resolution_status == "unresolved")
  if (n_unresolved > 0) {
    message(n_unresolved, " of ", nrow(out),
            " name(s) did not resolve -- see resolution_status column.")
  }
  n_syn <- sum(out$name_match_type == "synonym", na.rm = TRUE)
  if (n_syn > 0) {
    message(n_syn, " of ", nrow(out),
            " name(s) resolved via synonym (auto-updated to current canonical NCBI taxonomy).")
  }

  out
}
