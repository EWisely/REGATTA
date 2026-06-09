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

test_that("$result is ASV_id + the winning pct_id + 7 ranks, for all 12 ASVs", {
  expect_identical(names(r), c("ASV_id", "pct_id", ranks))
  expect_equal(nrow(r), 12L)
  expect_type(r$pct_id, "double")
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

test_that("global_lca_to_local counts all LCA applications; breakdown counts the real downgrades", {
  mk <- function(id, lin, pid) {
    d <- as.data.frame(as.list(stats::setNames(lin, ranks)), stringsAsFactors = FALSE)
    d$ASV_id <- id; d$pct_id <- pid; d
  }
  ep <- function(sp) c("Eukaryota","Chordata","Actinopteri","Perciformes",
                       "Serranidae","Epinephelus", sp)
  # AGREE: both call the same species, global higher %ID -> triggered, no downgrade
  # DOWNGRADE: disagree at species (same genus), global higher -> LCA to genus
  g <- rbind(mk("agree", ep("Epinephelus merra"),     99),
             mk("downg", ep("Epinephelus bontoides"), 99))
  l <- rbind(mk("agree", ep("Epinephelus merra"),     98),
             mk("downg", ep("Epinephelus merra"),     98))
  rec4 <- suppressMessages(reconcile_global_local(
    g, l, global_pct_id_scale = "0-100", local_pct_id_scale = "0-100"))
  tr <- rec4$tracking; s <- rec4$stats
  # BOTH are global_lca_to_local (the LCA was applied to both)
  expect_equal(tr$preferred_database[tr$ASV_id == "agree"], "global_lca_to_local")
  expect_equal(tr$preferred_database[tr$ASV_id == "downg"], "global_lca_to_local")
  expect_true(all(tr$global_lca_to_local_triggered))
  expect_equal(s$count[s$metric == "global_lca_to_local triggered"], 2)
  # the agreement keeps the species; the disagreement downgrades to genus
  expect_equal(rec4$result$species[rec4$result$ASV_id == "agree"], "Epinephelus merra")
  expect_true(is.na(rec4$result$species[rec4$result$ASV_id == "downg"]))
  # the downgrade breakdown counts ONLY the real downgrade (the agreement adds none)
  expect_equal(s$count[s$metric == "downgraded: species -> genus"], 1)
  expect_equal(sum(s$count[grepl("^downgraded:", s$metric)]), 1)
})

test_that("global_lca_to_local with no shared rank is consistently unassigned", {
  # global wins %ID but local disagrees at EVERY rank (cross-domain) -> the LCA
  # is empty. The row must be unassigned everywhere consistently, not counted as
  # assigned with an empty lineage.
  mk <- function(lin, pid) {
    d <- as.data.frame(as.list(stats::setNames(lin, ranks)), stringsAsFactors = FALSE)
    d$ASV_id <- "X"; d$pct_id <- pid; d
  }
  g <- mk(c("Eukaryota","Chordata","Actinopteri","Perciformes","Serranidae",
            "Epinephelus","Epinephelus bontoides"), 99)
  l <- mk(c("Bacteria","Firmicutes","Bacilli","Bacillales","Bacillaceae",
            "Bacillus","Bacillus subtilis"), 98)
  rec3 <- suppressMessages(reconcile_global_local(
    g, l, global_pct_id_scale = "0-100", local_pct_id_scale = "0-100"))
  expect_true(all(is.na(rec3$result[1, ranks])))         # empty preferred lineage
  expect_true(is.na(rec3$tracking$preferred_database))   # -> not "global_lca_to_local"
  expect_true(is.na(rec3$result$pct_id))                 # -> no carried pct_id
  expect_false(rec3$tracking$global_lca_to_local_triggered)
  expect_equal(rec3$stats$count[rec3$stats$metric == "assigned ASVs"], 0)
})

test_that("missing required columns error", {
  bad <- data.frame(ASV_id = "A1", domain = "Eukaryota", stringsAsFactors = FALSE)
  expect_error(reconcile_global_local(bad, bad, output_dir = NULL),
               "missing required columns")
})
