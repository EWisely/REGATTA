# Offline tests for summarize_regatta(): the 21-row stats summary.

mk_result <- function(species) {
  n <- length(species)
  data.frame(
    ASV_id  = paste0("A", seq_len(n)),
    domain  = "Eukaryota", phylum = "Chordata", class = "Actinopterygii",
    order   = "Scorpaeniformes", family = "Sebastidae", genus = "Sebastes",
    species = species,
    stringsAsFactors = FALSE
  )
}

test_that("summarize_regatta returns 21 rows with a column per supplied stage", {
  res <- mk_result(c("Sebastes mystinus", NA, NA))   # 1 to species, 2 to genus-only
  s <- summarize_regatta(post_checklist = list(result = res))

  expect_equal(nrow(s), 21L)
  expect_true(all(c("row_names", "regatta") %in% names(s)))
  expect_equal(s$regatta[s$row_names == "total ASVs"], 3)
  expect_equal(s$regatta[s$row_names == "assigned ASVs"], 3)
  expect_equal(s$regatta[s$row_names == "ID'ed to species"], 1)
  expect_equal(s$regatta[s$row_names == "ID'ed to genus only"], 2)
})

test_that("each supplied input becomes its own stage column", {
  res <- mk_result("Sebastes mystinus")
  s <- summarize_regatta(global_input = res, post_checklist = list(result = res))
  expect_true(all(c("global", "regatta") %in% names(s)))
})

test_that("summarize_regatta errors with no inputs and on a malformed list", {
  expect_error(summarize_regatta(), "Supply at least one")
  expect_error(summarize_regatta(post_checklist = list(foo = 1)), "result")
})
