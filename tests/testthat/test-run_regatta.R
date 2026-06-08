# End-to-end test of the run_regatta() wrapper on a single-DB vsearch
# lca + userout pair. No network or DB: the vsearch path parses purely, and a
# pre-taxonomized checklist makes the background taxonomize a no-op.

test_that("run_regatta dispatches a vsearch lca+userout pair through to the checklist LCA", {
  tmp <- tempfile(); dir.create(tmp)
  out <- file.path(tmp, "out")

  lca <- file.path(tmp, "mifish_lca.txt")
  writeLines(c(
    paste("ASV_1", "d:Eukaryota,p:Chordata,c:Actinopterygii,o:Scorpaeniformes,f:Sebastidae,g:Sebastes,s:Sebastes_mystinus", sep = "\t"),
    paste("ASV_2", "d:Eukaryota,p:Chordata,c:Actinopterygii,o:Scorpaeniformes,f:Sebastidae,g:Sebastes,s:Sebastes_goodei",   sep = "\t")
  ), lca)

  uo <- file.path(tmp, "mifish_userout.txt")
  writeLines(c(
    paste("ASV_1", "SEQ;tax=x", "99.5", "170", "1", "0", sep = "\t"),
    paste("ASV_2", "SEQ;tax=x", "98.0", "170", "2", "0", sep = "\t")
  ), uo)

  # Pre-taxonomized regional checklist: has Sebastes mystinus (and so genus
  # Sebastes), but NOT S. goodei.
  checklist <- data.frame(
    domain = "Eukaryota", phylum = "Chordata", class = "Actinopterygii",
    order = "Scorpaeniformes", family = "Sebastidae", genus = "Sebastes",
    species = "Sebastes mystinus", stringsAsFactors = FALSE)

  res <- suppressMessages(suppressWarnings(
    run_regatta(input = c(lca, uo), checklist = checklist, out_dir = out,
                region = "testreg", label = "fish")
  ))

  # ASV_1 on the checklist -> kept at species; ASV_2 species off, genus on
  # -> downgraded to genus.
  rr <- res$post_checklist$result
  expect_equal(rr$species[rr$ASV_id == "ASV_1"], "Sebastes mystinus")
  expect_true(is.na(rr$species[rr$ASV_id == "ASV_2"]))
  expect_equal(rr$genus[rr$ASV_id == "ASV_2"], "Sebastes")

  # 21-row summary returned
  expect_equal(nrow(res$summary), 21L)

  # output bundle written under a dated run subfolder of out_dir
  run <- file.path(out, paste0("testreg_fish_", Sys.Date()))
  expect_true(dir.exists(run))
  expect_true(file.exists(file.path(run, "regatta_summary.csv")))
  expect_true(file.exists(file.path(run, "run_log.txt")))
  expect_true(file.exists(file.path(run, "reconcile_checklist",
                                    "reconcile_checklist_taxonomy_table.csv")))
})

test_that("run_regatta rejects an unsupported checklist argument", {
  expect_error(
    suppressMessages(run_regatta(input = "x.tab", checklist = 42,
                                 out_dir = tempfile(),
                                 region = "r", label = "l")),
    "checklist must be")
})

test_that("run_regatta requires out_dir, region, and label", {
  expect_error(
    suppressMessages(run_regatta(input = "x.tab", checklist = 42,
                                 region = "r", label = "l")),
    "out_dir")
  expect_error(
    suppressMessages(run_regatta(input = "x.tab", checklist = 42,
                                 out_dir = tempfile(), label = "l")),
    "region")
})

test_that("a pre-resolved 7-rank CSV (capitalized ranks, empty id) is detected and run", {
  tmp <- tempfile(); dir.create(tmp)
  csv <- file.path(tmp, "ranks.csv")
  utils::write.csv(data.frame(
    ASV_id  = c("", ""),                     # empty id column -> synthesized
    Domain = "Eukaryota", Phylum = "Chordata", Class = "Actinopterygii",
    Order = "Scorpaeniformes", Family = "Sebastidae", Genus = "Sebastes",
    Species = c("Sebastes mystinus", "Sebastes goodei"),
    pct_id = c(99.5, 97.2), stringsAsFactors = FALSE), csv, row.names = FALSE)

  expect_equal(REGATTA:::.regatta_detect_format(csv), "ranks_csv")

  checklist <- data.frame(
    domain = "Eukaryota", phylum = "Chordata", class = "Actinopterygii",
    order = "Scorpaeniformes", family = "Sebastidae", genus = "Sebastes",
    species = "Sebastes mystinus", stringsAsFactors = FALSE)

  res <- suppressMessages(suppressWarnings(run_regatta(
    input = csv, checklist = checklist, out_dir = file.path(tmp, "out"),
    region = "r", label = "l")))
  rr <- res$post_checklist$result
  expect_equal(nrow(rr), 2L)
  expect_true(all(nzchar(rr$ASV_id)))                 # ids were synthesized
  expect_false(is.na(rr$species[1]))                  # mystinus on checklist
  expect_true(is.na(rr$species[2]) && rr$genus[2] == "Sebastes")  # goodei -> genus
})

test_that("an unrecognized CSV gives an informative error listing columns + formats", {
  bad <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(foo = 1, bar = 2), bad, row.names = FALSE)
  err <- tryCatch(REGATTA:::.regatta_read_input(bad, NULL, "ASV_id"),
                  error = function(e) conditionMessage(e))
  expect_match(err, "Columns found: foo, bar")
  expect_match(err, "pre-resolved taxonomy CSV")
})
