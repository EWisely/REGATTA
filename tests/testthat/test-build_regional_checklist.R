# Offline tests for the build_regional_checklist() orchestrator. It now
# returns the lists and writes nothing unless output_dir is given, so no
# here::here() redirection is needed. Real taxonomize (DB) tests are guarded
# by skip_on_cran() + skip-if-no-DB.

ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

test_that("CSV source splits genus-only from species in the returned lists", {
  csv <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(Genus   = c("Lutjanus", "Caranx", "Sebastes"),
               Species = c("argentiventris", "sp.", "mystinus")),
    csv, row.names = FALSE)

  res <- suppressMessages(build_regional_checklist(
    region = "testreg", label = "testlab",
    OBIS = FALSE, GBIF = FALSE, CSV = csv, sql_path = tempfile()))   # no DB

  # DB list: a bare character vector of species binomials, NO bare genus
  expect_type(res$for_making_localdb, "character")
  expect_null(names(res$for_making_localdb))   # no column names, just names
  expect_setequal(res$for_making_localdb,
                  c("Lutjanus argentiventris", "Sebastes mystinus"))
  expect_false("Caranx" %in% res$for_making_localdb)
  # LCA list: + the retained bare genus. No DB -> for_LCA is the raw name list
  # (Species/Source), not a taxonomized rank table; checklist_summary is NULL.
  expect_true("Caranx" %in% res$for_LCA$Species)
  expect_false("genus" %in% names(res$for_LCA))
  expect_null(res$checklist_summary)    # no DB -> no taxonomized summary
})

test_that("nothing is written to the working directory by default", {
  csv <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(Genus = "Lutjanus", Species = "argentiventris"),
                   csv, row.names = FALSE)
  wd <- tempfile(); dir.create(wd); old <- setwd(wd); on.exit(setwd(old))

  suppressMessages(build_regional_checklist(
    region = "r", label = "l", OBIS = FALSE, CSV = csv, sql_path = tempfile()))
  expect_equal(length(list.files(wd, recursive = TRUE)), 0L)
})

test_that("output_dir writes the two lists with region/label filenames", {
  csv <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(Genus = "Lutjanus", Species = "argentiventris"),
                   csv, row.names = FALSE)
  out <- tempfile(); dir.create(out)

  suppressMessages(build_regional_checklist(
    region = "testreg", label = "testlab",
    OBIS = FALSE, CSV = csv, sql_path = tempfile(), output_dir = out))

  expect_true(file.exists(file.path(out, "comprehensive_testreg_testlab_list_for_making_localdb.txt")))
  expect_true(file.exists(file.path(out, "comprehensive_testreg_testlab_list_for_LCA.txt")))
})

test_that("pre-made OBIS/GBIF source CSVs are fed in and stacked", {
  obis <- tempfile(fileext = ".csv"); gbif <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(Species = "Aaa bbb", Source = "OBIS"), obis, row.names = FALSE)
  utils::write.csv(data.frame(Species = "Ccc ddd", Source = "GBIF"), gbif, row.names = FALSE)

  res <- suppressMessages(build_regional_checklist(
    region = "r", label = "l", OBIS = obis, GBIF = gbif, sql_path = tempfile()))
  expect_true(all(c("Aaa bbb", "Ccc ddd") %in% res$for_LCA$Species))
})

test_that("region/label are required and at least one source must be on", {
  expect_error(build_regional_checklist(label = "l", OBIS = FALSE, CSV = FALSE), "region")
  expect_error(build_regional_checklist(region = "r", OBIS = FALSE, CSV = FALSE), "label")
  expect_error(
    build_regional_checklist(region = "r", label = "l",
                             OBIS = FALSE, GBIF = FALSE, CSV = FALSE),
    "No species")
})

test_that("a fresh download requires taxa and regional_poly", {
  expect_error(build_regional_checklist(region = "r", label = "l"),
               "taxa.*regional_poly|regional_poly")
})

test_that("for_LCA is taxID + ranks, checklist_summary holds the audit columns", {
  skip_on_cran()
  sql <- Sys.getenv("REGATTA_TAXONOMIZR_SQL", "accessionTaxa.sql")
  skip_if(!file.exists(sql), "taxonomizr DB not available")

  csv <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(Genus = "Lutjanus", Species = "argentiventris"),
                   csv, row.names = FALSE)
  res <- suppressMessages(build_regional_checklist(
    region = "r", label = "l", OBIS = FALSE, CSV = csv, sql_path = sql))

  # for_LCA: taxID + the 7 ranks only -- audit columns stripped out
  expect_equal(names(res$for_LCA), c("taxID", ranks))
  expect_false("input_name" %in% names(res$for_LCA))
  expect_equal(res$for_LCA$genus[res$for_LCA$species == "Lutjanus argentiventris"],
               "Lutjanus")
  # checklist_summary: the full table with the resolution-status columns
  expect_true(all(c("input_name", "resolution_status", "name_match_type", "taxID")
                  %in% names(res$checklist_summary)))
  expect_equal(
    res$checklist_summary$genus[res$checklist_summary$input_name == "Lutjanus argentiventris"],
    "Lutjanus")
})

test_that(".regatta_ensure_taxonomized passes a taxonomized checklist through untouched", {
  cl <- data.frame(input_name = "x", stringsAsFactors = FALSE)
  cl[ranks] <- NA
  cl$genus  <- "Caranx"
  expect_silent(out <- .regatta_ensure_taxonomized(cl, "ignored"))
  expect_identical(out, cl)
})

test_that(".regatta_ensure_taxonomized warns and taxonomizes a raw checklist", {
  skip_on_cran()
  sql <- Sys.getenv("REGATTA_TAXONOMIZR_SQL", "accessionTaxa.sql")
  skip_if(!file.exists(sql), "taxonomizr DB not available")

  raw <- data.frame(Species = "Lutjanus argentiventris", Source = "Local_csv",
                    stringsAsFactors = FALSE)
  expect_warning(out <- suppressMessages(.regatta_ensure_taxonomized(raw, sql)),
                 "not pre-taxonomized")
  expect_true(all(ranks %in% names(out)))
  expect_equal(out$genus[out$input_name == "Lutjanus argentiventris"], "Lutjanus")
})
