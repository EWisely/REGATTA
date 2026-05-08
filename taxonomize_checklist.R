# taxonomize_checklist.R
# Eldridge Wisely

# Resolve a list of scientific names into a 7-rank taxonomy table using
# taxonomizr's local NCBI taxonomy SQL DB. Cleans common junk strings
# (Gen., indet., sp., cf., quotes) before lookup. Returns one row per
# unique input name with `input_name`, `resolution_status`, and the 7
# lowercase rank columns (domain → species). Unresolved names stay in
# the output with NA ranks and resolution_status = "unresolved".

# Intended to run once per region per taxonomic group. The resulting
# data.frame is the checklist input to regatta_checklist_lca().

taxonomize_checklist <- function(input,
                                 sql_path = "accessionTaxa.sql",
                                 prepare_db = !file.exists(sql_path)) {
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
              " already exists — skipping rebuild.")
    } else {
      message("Building taxonomizr SQL DB at ", sql_path,
              " — this takes ~15 min and downloads several GB.")
      taxonomizr::prepareDatabase(sql_path)
    }
  } else if (!file.exists(sql_path)) {
    stop("SQL DB not found at ", sql_path,
         ". Pass `prepare_db = TRUE` to build, or provide the path to an existing DB.")
  }

  ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

  taxids <- taxonomizr::getId(cleaned, sql_path)
  # Request domain explicitly — getTaxonomy's default asks for "superkingdom",
  # which is NA in current NCBI dumps where the top rank is named "domain".
  taxa_mat <- taxonomizr::getTaxonomy(taxids, sql_path, desiredTaxa = ranks)
  taxa_df <- as.data.frame(taxa_mat, stringsAsFactors = FALSE)

  out <- cbind(
    data.frame(
      input_name        = cleaned,
      resolution_status = ifelse(is.na(taxids), "unresolved", "resolved"),
      taxID             = taxids,
      stringsAsFactors  = FALSE
    ),
    taxa_df[, ranks, drop = FALSE]
  )
  rownames(out) <- NULL

  n_unresolved <- sum(out$resolution_status == "unresolved")
  if (n_unresolved > 0) {
    message(n_unresolved, " of ", nrow(out),
            " name(s) did not resolve — see resolution_status column.")
  }

  out
}
