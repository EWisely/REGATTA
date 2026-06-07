# Synthetic 6-row test for reconcile_checklist().

# Tiny regional checklist: 3 species across 2 families and 2 orders.
checklist <- data.frame(
  domain  = "Eukaryota",
  phylum  = "Chordata",
  class   = "Actinopterygii",
  order   = c("Scorpaeniformes", "Scorpaeniformes", "Clupeiformes"),
  family  = c("Sebastidae",      "Sebastidae",      "Engraulidae"),
  genus   = c("Sebastes",        "Sebastes",        "Engraulis"),
  species = c("Sebastes mystinus","Sebastes paucispinis","Engraulis mordax"),
  stringsAsFactors = FALSE
)

input <- data.frame(
  ASV_id  = paste0("ASV_", 1:6),
  domain  = c("Eukaryota","Eukaryota","Eukaryota","Eukaryota","Bacteria","Eukaryota"),
  phylum  = c("Chordata","Chordata","Chordata","Chordata","Firmicutes","Chordata"),
  class   = c("Actinopterygii","Actinopterygii","Actinopterygii","Actinopterygii","Bacilli","Actinopterygii"),
  order   = c("Scorpaeniformes","Scorpaeniformes","Clupeiformes","Clupeiformes","Bacillales","Clupeiformes"),
  family  = c("Sebastidae","Sebastidae","Clupeidae","Engraulidae","Bacillaceae",NA),
  genus   = c("Sebastes","Sebastes","Sardinops","Engraulis","Bacillus",NA),
  species = c("Sebastes mystinus","Sebastes goodei","Sardinops sagax",
              "Engraulis mordax","Bacillus subtilis",NA),
  pct_id  = c(99.5, 97.2, 95.0, 100.0, 98.0, 88.0),
  stringsAsFactors = FALSE
)

result <- reconcile_checklist(input, checklist, output_dir = NULL)

test_that("$result is strict 8-column shape", {
  expect_identical(
    names(result$result),
    c("ASV_id","domain","phylum","class","order","family","genus","species")
  )
})

test_that("species-level match (ASV_1) is kept", {
  expect_equal(result$result$species[1], "Sebastes mystinus")
})

test_that("species miss + genus hit (ASV_2) downgrades to genus", {
  expect_true(is.na(result$result$species[2]))
  expect_equal(result$result$genus[2], "Sebastes")
})

test_that("multi-rank downgrade (ASV_3) lands at order", {
  expect_true(is.na(result$result$species[3]))
  expect_true(is.na(result$result$genus[3]))
  expect_true(is.na(result$result$family[3]))
  expect_equal(result$result$order[3], "Clupeiformes")
})

test_that("no-regional-record (ASV_5) ends all-NA", {
  expect_true(is.na(result$result$domain[5]))
  expect_true(is.na(result$result$species[5]))
})

test_that("$tracking carries regatta_match_rank + before/after", {
  expect_true("regatta_match_rank" %in% names(result$tracking))
  expect_true("before_species"     %in% names(result$tracking))
  expect_true("after_species"      %in% names(result$tracking))
  expect_equal(result$tracking$regatta_match_rank[1], "species")
  expect_equal(result$tracking$regatta_match_rank[2], "genus")
  expect_equal(result$tracking$regatta_match_rank[3], "order")
  expect_true(is.na(result$tracking$regatta_match_rank[5]))
})

test_that("$stats has expected counts", {
  expect_equal(result$stats$count[result$stats$metric == "total ASVs"], 6)
  expect_equal(result$stats$count[result$stats$metric == "matched at species"], 2)
  expect_equal(result$stats$count[result$stats$metric == "not matched (no regional record)"], 1)
})
