# Tests for summarize_regatta() on real data: the bundled 12-ASV Galapagos
# MiFish obitools slice run through reconcile_checklist against the bundled
# taxonomized checklist. No synthetic taxa.

global    <- read.csv(system.file("extdata", "mifish_obitools_example.csv", package = "REGATTA"),
                      stringsAsFactors = FALSE)
checklist <- readRDS(system.file("extdata", "galapagos_fish_checklist.rds", package = "REGATTA"))
post      <- suppressMessages(reconcile_checklist(global, checklist, output_dir = NULL))

test_that("summarize_regatta produces the 24-row table with real counts", {
  s <- summarize_regatta(post_checklist = post)
  val <- function(label) s$regatta_result[s$row_names == label]

  expect_equal(nrow(s), 24L)
  expect_true(all(c("row_names", "regatta_result") %in% names(s)))
  expect_equal(val("total ASVs"), 12)
  expect_equal(val("assigned ASVs"), 12)
  expect_equal(val("ID'ed to species"), 6)        # 6 of the 12 stay at species
  expect_equal(val("ID'ed to genus only"), 4)
  # checklist-membership rows present but NA without a checklist
  expect_true(is.na(val("percent of ASVs on regional checklist")))
  expect_true(is.na(val("percent of distinct taxa on regional checklist")))
})

test_that("checklist membership: input < 100% on the checklist, regatta_result = 100%", {
  s <- summarize_regatta(global_input = global, post_checklist = post,
                         checklist = checklist)
  asv <- function(col) s[[col]][s$row_names == "percent of ASVs on regional checklist"]
  tax <- function(col) s[[col]][s$row_names == "percent of distinct taxa on regional checklist"]
  expect_lt(asv("global"), 100)              # some global species aren't on the checklist
  expect_equal(asv("regatta_result"), 100)   # every reconciled call is on the checklist
  expect_equal(tax("regatta_result"), 100)

  # The reverse direction: checklist recovery. Identical before/after at species
  # level (REGATTA never drops a true regional detection), and > 0 here.
  rec <- function(col) s[[col]][s$row_names == "percent of checklist species detected"]
  expect_gt(rec("global"), 0)
  expect_equal(rec("global"), rec("regatta_result"))
})

test_that("columns are the input(s) + regatta_result; single-DB input is 'input_file'", {
  s2 <- summarize_regatta(global_input = global, post_checklist = post)
  expect_true(all(c("global", "regatta_result") %in% names(s2)))
  expect_false("reconciled" %in% names(s2))   # merge intermediate is not a column
  s1 <- summarize_regatta(input_file = global, post_checklist = post)
  expect_identical(setdiff(names(s1), "row_names"), c("input_file", "regatta_result"))
})

test_that("summarize_regatta errors with no inputs and on a malformed list", {
  expect_error(summarize_regatta(), "Supply at least one")
  expect_error(summarize_regatta(post_checklist = list(foo = 1)), "result")
})
