# Offline tests for parse_sintax(): pure string -> 7-rank parsing.

ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

test_that("a full SINTAX string splits into all 7 ranks (species underscores -> spaces)", {
  s <- "d:Eukaryota,p:Chordata,c:Actinopteri,o:Scombriformes,f:Scombridae,g:Auxis,s:Auxis_thazard"
  out <- parse_sintax(s)
  expect_true(all(c("sintax", ranks) %in% names(out)))
  expect_equal(out$domain,  "Eukaryota")
  expect_equal(out$genus,   "Auxis")
  expect_equal(out$species, "Auxis thazard")
})

test_that("truncated strings leave deeper ranks NA", {
  out <- parse_sintax("d:Eukaryota,p:Chordata,c:Actinopteri,o:Scombriformes,f:Scombridae,g:Auxis")
  expect_equal(out$genus, "Auxis")
  expect_true(is.na(out$species))
})

test_that("empty and NA inputs return all-NA ranks", {
  out <- parse_sintax(c("", NA))
  expect_true(all(is.na(out[1, ranks])))
  expect_true(all(is.na(out[2, ranks])))
})

test_that("a data.frame input needs sintax_col and carries other columns through", {
  df <- data.frame(ASV_id = "A1",
                   tax = "d:Eukaryota,g:Auxis,s:Auxis_thazard",
                   stringsAsFactors = FALSE)
  out <- parse_sintax(df, sintax_col = "tax")
  expect_equal(out$ASV_id,  "A1")
  expect_equal(out$species, "Auxis thazard")
  expect_error(parse_sintax(df), "sintax_col")
})

test_that("rank-name collisions get a resolved_ prefix", {
  df <- data.frame(species = "x",
                   tax = "d:Eukaryota,s:Auxis_thazard",
                   stringsAsFactors = FALSE)
  out <- parse_sintax(df, sintax_col = "tax")
  expect_true("resolved_species" %in% names(out))
  expect_equal(out$resolved_species, "Auxis thazard")
})
