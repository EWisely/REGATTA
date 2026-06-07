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
                                 out_dir = tempfile())),
    "checklist must be")
})

test_that("run_regatta requires out_dir", {
  expect_error(suppressMessages(run_regatta(input = "x.tab", checklist = 42)),
               "out_dir")
})
