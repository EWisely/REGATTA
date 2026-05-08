# parse_sintax.R
# Eldridge Wisely

# Parse a vector of vsearch SINTAX taxonomy strings into a 7-column
# data.frame of ranks. SINTAX strings look like:
#   d:Eukaryota,p:Chordata,c:Actinopteri,o:Scombriformes,f:Scombridae,g:Auxis,s:Auxis_thazard
# Truncated strings (e.g. only down to genus) are accepted; missing
# ranks come back as NA. Empty strings and NA inputs return all-NA rows.
# Species values are returned with underscores converted to spaces.

# Output is a data.frame with columns: domain, phylum, class, order,
# family, genus, species. The user is responsible for binding it to
# whatever ID column they have.

parse_sintax <- function(sintax) {
  ranks    <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  prefixes <- c("d",      "p",      "c",     "o",     "f",      "g",     "s")

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

  result <- data.frame(
    setNames(lapply(prefixes, extract_rank, s = sintax), ranks),
    stringsAsFactors = FALSE
  )

  # Empty strings become NA; species underscores become spaces
  for (r in ranks) {
    empty <- !is.na(result[[r]]) & !nzchar(result[[r]])
    result[[r]][empty] <- NA_character_
  }
  result$species <- gsub("_", " ", result$species, fixed = TRUE)

  result
}
