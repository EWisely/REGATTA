# Tests for reconcile_global_local() on the bundled real-data fixtures: the
# same 12 Galapagos MiFish ASVs classified two ways that the two-database
# vignette runs. No synthetic taxa, no network.

ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

global <- read.csv(system.file("extdata", "mifish_obitools_example.csv", package = "REGATTA"),
                   stringsAsFactors = FALSE)
local  <- read.csv(system.file("extdata", "mifish_vsearch_example.csv",  package = "REGATTA"),
                   stringsAsFactors = FALSE)

rec <- reconcile_global_local(global, local, output_dir = NULL)
r   <- rec$result
tr  <- rec$tracking
# Look ASVs up by their (real) global species call.
asv  <- function(sp) global$ASV_id[match(sp, global$species)]
rrow <- function(sp) which(r$ASV_id  == asv(sp))
trow <- function(sp) which(tr$ASV_id == asv(sp))

test_that("$result is the strict 8-column exchange format for all 12 ASVs", {
  expect_identical(names(r), c("ASV_id", ranks))
  expect_equal(nrow(r), 12L)
})

test_that("a global-only ASV keeps the global call (local made none)", {
  expect_equal(tr$preferred_database[trow("Abudefduf taurus")], "global")
  expect_equal(r$species[rrow("Abudefduf taurus")], "Abudefduf taurus")
})

test_that("an agreeing ASV keeps the species (local wins the tie by default)", {
  expect_equal(tr$preferred_database[trow("Lutjanus kasmira")], "local")
  expect_equal(r$species[rrow("Lutjanus kasmira")], "Lutjanus kasmira")
})

test_that("global wins on %ID but disagrees -> LCA downgrade to genus (Mugil)", {
  expect_equal(tr$preferred_database[trow("Mugil curema")], "global_lca_to_local")
  expect_equal(r$genus[rrow("Mugil curema")], "Mugil")
  expect_true(is.na(r$species[rrow("Mugil curema")]))
})

test_that("local wins on %ID -> the local call is kept (Etrumeus teres)", {
  expect_equal(tr$preferred_database[trow("Etrumeus micropus")], "local")
  expect_equal(r$species[rrow("Etrumeus micropus")], "Etrumeus teres")
})

test_that("$stats counts the 12 ASVs and the two LCA downgrades", {
  cnt <- function(m) rec$stats$count[rec$stats$metric == m]
  expect_equal(cnt("total ASVs"), 12)
  expect_equal(cnt("global_lca_to_local triggered"), 2)   # Mugil, Epinephelus
})

test_that("a 0-1 global pct_id is auto-rescaled before comparison", {
  g01 <- global; g01$pct_id <- g01$pct_id / 100          # the real values, on a 0-1 scale
  rec2 <- suppressMessages(reconcile_global_local(g01, local, output_dir = NULL))
  # Rescaling restores the same decision as the 0-100 run.
  i <- which(rec2$tracking$ASV_id == asv("Mugil curema"))
  expect_equal(rec2$tracking$preferred_database[i], "global_lca_to_local")
})

test_that("missing required columns error", {
  bad <- data.frame(ASV_id = "A1", domain = "Eukaryota", stringsAsFactors = FALSE)
  expect_error(reconcile_global_local(bad, bad, output_dir = NULL),
               "missing required columns")
})
