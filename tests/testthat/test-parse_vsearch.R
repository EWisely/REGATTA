# Offline tests for parse_vsearch_userout() and parse_vsearch_results(),
# using temp files in the documented column layouts.

test_that("parse_vsearch_userout keeps the first hit per ASV with its taxonomy + pct_id", {
  uo <- tempfile()
  writeLines(c(
    paste("A1", "SEQ1;tax=d:Eukaryota,g:Auxis,s:Auxis_thazard", "99.4", "170", "1", "0", sep = "\t"),
    paste("A1", "SEQ2;tax=d:Eukaryota,g:Auxis,s:Auxis_rochei",  "99.0", "170", "2", "0", sep = "\t"),
    paste("A2", "SEQ3;tax=d:Eukaryota,g:Sardinops,s:Sardinops_sagax", "97.1", "170", "5", "0", sep = "\t")
  ), uo)

  out <- suppressMessages(parse_vsearch_userout(uo))
  expect_equal(nrow(out), 2L)                    # one row per ASV
  expect_equal(out$ASV_id, c("A1", "A2"))
  expect_equal(out$species[out$ASV_id == "A1"], "Auxis thazard")   # first hit kept
  expect_equal(out$pct_id[out$ASV_id == "A1"], 99.4)
})

test_that("parse_vsearch_results uses LCA taxonomy + userout pct_id", {
  lca <- tempfile(); uo <- tempfile()
  writeLines(c(
    paste("A1", "d:Eukaryota,p:Chordata,g:Auxis", sep = "\t"),    # LCA consensus: only to genus
    paste("A2", "d:Eukaryota,g:Sardinops,s:Sardinops_sagax", sep = "\t")
  ), lca)
  writeLines(c(
    paste("A1", "SEQ;tax=d:Eukaryota,g:Auxis,s:Auxis_thazard", "99.4", "170", "1", "0", sep = "\t"),
    paste("A2", "SEQ;tax=d:Eukaryota,g:Sardinops,s:Sardinops_sagax", "97.1", "170", "5", "0", sep = "\t")
  ), uo)

  out <- parse_vsearch_results(lca, uo)
  expect_equal(out$ASV_id, c("A1", "A2"))
  expect_equal(out$genus[out$ASV_id == "A1"], "Auxis")
  expect_true(is.na(out$species[out$ASV_id == "A1"]))   # taxonomy from the conservative LCA
  expect_equal(out$pct_id[out$ASV_id == "A1"], 99.4)    # pct_id from userout
})

test_that("parse_vsearch_results gives NA pct_id for an LCA ASV absent from userout", {
  lca <- tempfile(); uo <- tempfile()
  writeLines(paste("A9", "d:Eukaryota,g:Auxis", sep = "\t"), lca)
  writeLines(paste("A1", "SEQ;tax=d:Eukaryota,g:Auxis", "99.0", "170", "1", "0", sep = "\t"), uo)
  out <- parse_vsearch_results(lca, uo)
  expect_equal(out$ASV_id, "A9")
  expect_true(is.na(out$pct_id))
})

test_that("missing input files error clearly", {
  expect_error(suppressMessages(parse_vsearch_userout("/no/such/userout.txt")), "not found")
  expect_error(parse_vsearch_results("/no/such/lca.txt", "/no/such/uo.txt"), "not found")
})
