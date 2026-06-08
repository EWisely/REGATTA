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

test_that("$stats captures the before -> after specificity story", {
  s <- result$stats; g <- function(m) s$count[s$metric == m]
  expect_equal(g("total ASVs"), 6)
  expect_equal(g("assigned before checklist-LCA"), 6)
  expect_equal(g("assigned after checklist-LCA"), 5)        # ASV_5 (Bacillus) dropped
  # unchanged + downgraded + dropped partition the assigned-before set
  expect_equal(g("ASVs unchanged (specificity kept)"), 3)   # ASV_1, ASV_4, ASV_6
  expect_equal(g("ASVs downgraded (specificity reduced)"), 2) # ASV_2 -> genus, ASV_3 -> order
  expect_equal(g("no regional record (call dropped)"), 1)   # ASV_5
  expect_equal(g("ASVs unchanged (specificity kept)") +
               g("ASVs downgraded (specificity reduced)") +
               g("no regional record (call dropped)"),
               g("assigned before checklist-LCA"))
  # the per-rank before/after distribution is NOT in $stats anymore -- it lives
  # in summarize_regatta()'s input/regatta_result columns.
  expect_false(any(grepl("ID'ed to", s$metric)))
})

test_that("assigned-before/after count any non-NA rank (not just domain)", {
  # A pre-resolved table with NO domain column-equivalent (domain all NA) but
  # populated phylum..species: these ASVs are assigned and must be counted.
  nodomain <- data.frame(
    ASV_id = paste0("A", 1:3), domain = NA_character_, phylum = "Chordata",
    class = "Actinopterygii", order = "Scorpaeniformes", family = "Sebastidae",
    genus = "Sebastes",
    species = c("Sebastes mystinus", "Sebastes paucispinis", "Engraulis mordax"),
    pct_id = c(99, 98, 97), stringsAsFactors = FALSE)
  s <- suppressMessages(reconcile_checklist(nodomain, checklist, output_dir = NULL))$stats
  expect_equal(s$count[s$metric == "assigned before checklist-LCA"], 3)
  expect_equal(s$count[s$metric == "assigned after checklist-LCA"], 3)
  # an all-NA row is not counted as assigned
  with_blank <- rbind(nodomain, data.frame(
    ASV_id = "A4", domain = NA, phylum = NA, class = NA, order = NA,
    family = NA, genus = NA, species = NA, pct_id = NA))
  s2 <- suppressMessages(reconcile_checklist(with_blank, checklist, output_dir = NULL))$stats
  expect_equal(s2$count[s2$metric == "total ASVs"], 4)
  expect_equal(s2$count[s2$metric == "assigned before checklist-LCA"], 3)
})

test_that("target_group splits downgrades into non-local vs off-target", {
  ck <- data.frame(domain = "Eukaryota", phylum = "Chordata", class = "Actinopterygii",
                   order = "Scorpaeniformes", family = "Sebastidae", genus = "Sebastes",
                   species = "Sebastes mystinus", stringsAsFactors = FALSE)
  # target group = the fishes (class Actinopterygii), stamped as an attribute
  attr(ck, "target_group") <- data.frame(rank = "class", name = "Actinopterygii",
                                         stringsAsFactors = FALSE)
  tax <- data.frame(
    ASV_id = c("A1", "A2"), domain = "Eukaryota", phylum = "Chordata",
    class  = c("Actinopterygii", "Mammalia"),         # A1 fish, A2 mammal
    order  = c("Scorpaeniformes", "Cetacea"),
    family = c("Sebastidae", "Delphinidae"),
    genus  = c("Sebastes", "Tursiops"),
    species = c("Sebastes goodei", "Tursiops truncatus"),  # both off the checklist
    pct_id = c(98, 99), stringsAsFactors = FALSE)
  s <- suppressMessages(suppressWarnings(reconcile_checklist(tax, ck, output_dir = NULL)))$stats
  g <- function(m) s$count[s$metric == m]
  expect_equal(g("downgraded/dropped -- non-local (geographic)"), 1)  # the fish
  expect_equal(g("downgraded/dropped -- off-target (taxonomic)"), 1)  # the mammal
})

test_that("missing pct_id warns (mentions filtering + two-DB), suppressible via warn_pct_id", {
  no_pct <- input[, setdiff(names(input), "pct_id")]
  expect_warning(
    suppressMessages(reconcile_checklist(no_pct, checklist, output_dir = NULL)),
    "pct_id")
  # and the message points at both consequences
  w <- tryCatch(suppressMessages(reconcile_checklist(no_pct, checklist, output_dir = NULL)),
                warning = function(w) conditionMessage(w))
  expect_match(w, "filtering")
  expect_match(w, "two-database")
  # suppressed (run_regatta's two-DB step passes warn_pct_id = FALSE)
  expect_warning(
    suppressMessages(reconcile_checklist(no_pct, checklist, output_dir = NULL,
                                         warn_pct_id = FALSE)),
    regexp = NA)
})
