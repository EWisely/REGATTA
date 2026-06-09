# Tests for summarize_regatta() on real data: the bundled 12-ASV Galapagos
# MiFish obitools slice run through reconcile_checklist against the bundled
# taxonomized checklist. No synthetic taxa.

global    <- read.csv(system.file("extdata", "mifish_obitools_example.csv", package = "REGATTA"),
                      stringsAsFactors = FALSE)
checklist <- readRDS(system.file("extdata", "galapagos_fish_checklist.rds", package = "REGATTA"))
post      <- suppressMessages(reconcile_checklist(global, checklist, output_dir = NULL))

test_that("summarize_regatta produces the report with real counts + downgrade breakdown", {
  s <- summarize_regatta(post_checklist = post)
  val <- function(label) s$regatta_checklist_result[s$row_names == label]

  expect_gte(nrow(s), 27L)   # 24 base rows + transition headline + breakdown
  expect_true(all(c("row_names", "regatta_checklist_result") %in% names(s)))
  expect_equal(val("total ASVs"), 12)
  expect_equal(val("assigned ASVs"), 12)
  expect_equal(val("ID'ed to species"), 6)        # 6 of the 12 stay at species
  expect_equal(val("ID'ed to genus only"), 4)
  # transition headline + per-rank-pair downgrade breakdown
  expect_equal(val("ASVs downgraded (specificity reduced)"), 6)
  expect_equal(val("downgraded: species -> genus"), 4)
  expect_equal(val("downgraded: species -> family"), 2)
  # checklist-membership rows present but NA without a checklist
  expect_true(is.na(val("percent of ASVs on regional checklist")))
  expect_true(is.na(val("percent of distinct taxa on regional checklist")))
})

test_that("checklist membership: input < 100% on the checklist, checklist result = 100%", {
  s <- summarize_regatta(global_input = global, post_checklist = post,
                         checklist = checklist)
  asv <- function(col) s[[col]][s$row_names == "percent of ASVs on regional checklist"]
  tax <- function(col) s[[col]][s$row_names == "percent of distinct taxa on regional checklist"]
  expect_lt(asv("global"), 100)                        # some global species aren't on the checklist
  expect_equal(asv("regatta_checklist_result"), 100)   # every reconciled call is on the checklist
  expect_equal(tax("regatta_checklist_result"), 100)

  # The reverse direction: checklist recovery. Identical before/after at species
  # level (REGATTA never drops a true regional detection), and > 0 here.
  rec <- function(col) s[[col]][s$row_names == "percent of checklist species detected"]
  expect_gt(rec("global"), 0)
  expect_equal(rec("global"), rec("regatta_checklist_result"))
})

test_that("columns are the input(s) + per-step result; single-DB input is 'input_file'", {
  s2 <- summarize_regatta(global_input = global, post_checklist = post)
  expect_true(all(c("global", "regatta_checklist_result") %in% names(s2)))
  expect_false("regatta_global_local_result" %in% names(s2))   # no reconciled stage supplied
  s1 <- summarize_regatta(input_file = global, post_checklist = post)
  expect_identical(setdiff(names(s1), "row_names"), c("input_file", "regatta_checklist_result"))
})

test_that("transition rows sit in regatta_checklist_result with NA in the before column", {
  s <- summarize_regatta(input_file = global, post_checklist = post)
  tr <- function(lbl) s[s$row_names == lbl, c("input_file", "regatta_checklist_result")]
  for (lbl in c("ASVs unchanged (specificity kept)",
                "ASVs downgraded (specificity reduced)")) {
    row <- tr(lbl)
    expect_true(is.na(row$input_file))                  # NA in the "before" column
    expect_false(is.na(row$regatta_checklist_result))   # the count is in the checklist column
  }
  # nothing was dropped here (the walk reaches domain), so that row is suppressed
  expect_false(any(s$row_names == "no regional record (call dropped)"))
  # unchanged + downgraded == assigned (post-checklist), 12 here (0 dropped)
  vals <- vapply(c("ASVs unchanged (specificity kept)",
                   "ASVs downgraded (specificity reduced)"),
                 function(l) s$regatta_checklist_result[s$row_names == l], numeric(1))
  expect_equal(sum(vals), 12)
})

test_that("the 'call dropped' row appears only when a call really is dropped", {
  # an off-target ASV (Bacteria) that matches the fish checklist at no rank
  drop_in <- rbind(global, transform(global[1, ], ASV_id = "BACT_1",
    domain = "Bacteria", phylum = "Firmicutes", class = "Bacilli",
    order = "Bacillales", family = "Bacillaceae", genus = "Bacillus",
    species = "Bacillus subtilis"))
  pdrop <- suppressMessages(reconcile_checklist(drop_in, checklist, output_dir = NULL))
  s <- summarize_regatta(input_file = drop_in, post_checklist = pdrop)
  expect_true(any(s$row_names == "no regional record (call dropped)"))
  expect_equal(s$regatta_checklist_result[s$row_names == "no regional record (call dropped)"], 1)
})

# helper mirroring summarize_regatta's "assigned = any non-NA rank"
n_assigned_any_test <- function(t) {
  ranks <- c("domain","phylum","class","order","family","genus","species")
  sum(rowSums(!is.na(t[, ranks, drop = FALSE])) > 0)
}

test_that("two-DB has four columns; source-breakdown lives only in the global_local column", {
  local <- read.csv(system.file("extdata", "mifish_vsearch_example.csv", package = "REGATTA"),
                    stringsAsFactors = FALSE)
  rec  <- suppressMessages(reconcile_global_local(global, local, output_dir = NULL))
  post2 <- suppressMessages(reconcile_checklist(rec$result, checklist, output_dir = NULL))
  s <- summarize_regatta(rec, post2, global, local, checklist = checklist)
  expect_identical(setdiff(names(s), "row_names"),
                   c("global", "local", "regatta_global_local_result", "regatta_checklist_result"))
  # the source breakdown describes the global-vs-local step -> only that column
  for (lbl in c("count of local assignment preferred",
                "count of global assignment preferred",
                "change in number of ASVs assigned")) {
    row <- s[s$row_names == lbl, ]
    expect_true(is.na(row$global))
    expect_true(is.na(row$local))
    expect_false(is.na(row$regatta_global_local_result))   # value here
    expect_true(is.na(row$regatta_checklist_result))       # NA in the checklist column
  }
  # the transition headline is SHARED: both steps downgrade, so the row carries
  # a count in each step's column
  dn <- s[s$row_names == "ASVs downgraded (specificity reduced)", ]
  expect_false(is.na(dn$regatta_global_local_result))   # the LCA downgraded some
  expect_false(is.na(dn$regatta_checklist_result))       # the checklist downgraded some
  # the global-local-only headline stays in its own column
  tg <- s[s$row_names == "global_lca_to_local triggered", ]
  expect_false(is.na(tg$regatta_global_local_result))
  expect_true(is.na(tg$regatta_checklist_result))
  # local-preferred + global-preferred = ASVs the global DB assigned (denominator)
  lp <- s$regatta_global_local_result[s$row_names == "count of local assignment preferred"]
  gp <- s$regatta_global_local_result[s$row_names == "count of global assignment preferred"]
  expect_equal(lp + gp, n_assigned_any_test(global))
})

test_that("summarize_regatta errors with no inputs and on a malformed list", {
  expect_error(summarize_regatta(), "Supply at least one")
  expect_error(summarize_regatta(post_checklist = list(foo = 1)), "result")
})
