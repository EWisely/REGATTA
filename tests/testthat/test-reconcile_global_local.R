# Offline tests for reconcile_global_local(): synthetic global/local tables
# covering all four per-ASV outcomes plus pct_id rescaling and validation.

ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

# A1 global-only; A2 both agree; A3 disagree (global %ID higher);
# A4 disagree (local %ID higher).
global <- data.frame(
  ASV_id  = c("A1", "A2", "A3", "A4"),
  domain  = "Eukaryota", phylum = "Chordata", class = "Actinopterygii",
  order   = c("Scorpaeniformes", "Clupeiformes", "Scorpaeniformes", "Scorpaeniformes"),
  family  = c("Sebastidae", "Engraulidae", "Sebastidae", "Sebastidae"),
  genus   = c("Sebastes", "Engraulis", "Sebastes", "Sebastes"),
  species = c("Sebastes mystinus", "Engraulis mordax", "Sebastes mystinus", "Sebastes mystinus"),
  pct_id  = c(99, 98, 99, 95),
  stringsAsFactors = FALSE
)
# local has no call for A1
local <- data.frame(
  ASV_id  = c("A2", "A3", "A4"),
  domain  = "Eukaryota", phylum = "Chordata", class = "Actinopterygii",
  order   = c("Clupeiformes", "Scorpaeniformes", "Clupeiformes"),
  family  = c("Engraulidae", "Sebastidae", "Clupeidae"),
  genus   = c("Engraulis", "Sebastes", "Sardinops"),
  species = c("Engraulis mordax", "Sebastes paucispinis", "Sardinops sagax"),
  pct_id  = c(98, 97, 99),
  stringsAsFactors = FALSE
)

rec <- reconcile_global_local(global, local, output_dir = NULL)
r  <- rec$result
tr <- rec$tracking
ri <- function(a) which(r$ASV_id  == a)
ti <- function(a) which(tr$ASV_id == a)

test_that("$result is the strict 8-column exchange format", {
  expect_identical(names(r), c("ASV_id", ranks))
})

test_that("global-only ASV keeps the global call", {
  expect_equal(r$species[ri("A1")], "Sebastes mystinus")
  expect_equal(tr$preferred_database[ti("A1")], "global")
})

test_that("agreeing ASV keeps the species (local wins the tie by default)", {
  expect_equal(r$species[ri("A2")], "Engraulis mordax")
  expect_equal(tr$preferred_database[ti("A2")], "local")
})

test_that("disagreement with global %ID winning downgrades to the LCA (genus)", {
  expect_true(is.na(r$species[ri("A3")]))
  expect_equal(r$genus[ri("A3")], "Sebastes")
  expect_equal(tr$preferred_database[ti("A3")], "global_lca_to_local")
})

test_that("disagreement with local %ID winning keeps the local call", {
  expect_equal(r$species[ri("A4")], "Sardinops sagax")
  expect_equal(tr$preferred_database[ti("A4")], "local")
})

test_that("$stats reports the expected counts", {
  cnt <- function(m) rec$stats$count[rec$stats$metric == m]
  expect_equal(cnt("total ASVs"), 4)
  expect_equal(cnt("global_lca_to_local triggered"), 1)
  expect_equal(cnt("best_pctid: local won"), 2)   # A2, A4
  expect_equal(cnt("best_pctid: global won"), 2)   # A1, A3
})

test_that("Local_advantage = FALSE flips ties to global", {
  rec2 <- reconcile_global_local(global, local, Local_advantage = FALSE, output_dir = NULL)
  # A2 ties at 98/98; now global wins the tie -> global won, local also assigned
  # -> global_lca_to_local; the two agree at species so the species survives.
  expect_equal(rec2$tracking$preferred_database[which(rec2$tracking$ASV_id == "A2")],
               "global_lca_to_local")
  expect_equal(rec2$result$species[which(rec2$result$ASV_id == "A2")], "Engraulis mordax")
})

test_that("a 0-1 global pct_id is auto-rescaled to 0-100 before comparison", {
  g <- global[global$ASV_id == "A3", ]; g$pct_id <- 0.99   # 0-1 scale -> 99
  l <- local[local$ASV_id == "A3", ]                       # pct_id 97 (0-100)
  rec3 <- suppressMessages(reconcile_global_local(g, l, output_dir = NULL))
  expect_equal(rec3$tracking$preferred_pctid[1], 99)
  expect_equal(rec3$tracking$preferred_database[1], "global_lca_to_local")
  expect_equal(rec3$result$genus[1], "Sebastes")
})

test_that("missing required columns error", {
  bad <- data.frame(ASV_id = "A1", domain = "Eukaryota", stringsAsFactors = FALSE)
  expect_error(reconcile_global_local(bad, bad, output_dir = NULL),
               "missing required columns")
})
