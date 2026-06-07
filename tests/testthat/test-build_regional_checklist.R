# Offline tests for the build_regional_checklist() orchestrator. output_dir is
# now REQUIRED, and sql_path = NULL skips taxonomization (deferring it to the
# reconcile step) -- the offline tests use that to stay free of a taxonomizr DB
# and network. Real taxonomize (DB) tests are guarded by skip_on_cran() +
# skip-if-no-DB. We write to tempfile() dirs, so no here::here() redirection.

ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

test_that("CSV source splits genus-only from species in the returned lists", {
  csv <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(Genus   = c("Lutjanus", "Caranx", "Sebastes"),
               Species = c("argentiventris", "sp.", "mystinus")),
    csv, row.names = FALSE)

  res <- suppressMessages(build_regional_checklist(
    region = "testreg", label = "testlab", output_dir = tempfile(),
    OBIS = FALSE, GBIF = FALSE, CSV = csv, sql_path = NULL))   # skip taxonomization

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

test_that("output_dir is required", {
  csv <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(Genus = "Lutjanus", Species = "argentiventris"),
                   csv, row.names = FALSE)
  expect_error(
    build_regional_checklist(region = "r", label = "l", OBIS = FALSE, CSV = csv,
                             sql_path = NULL),
    "output_dir")
})

test_that("output_dir writes the lists + methods with region/label filenames", {
  csv <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(Genus = "Lutjanus", Species = "argentiventris"),
                   csv, row.names = FALSE)
  out <- tempfile()

  suppressMessages(build_regional_checklist(
    region = "testreg", label = "testlab",
    OBIS = FALSE, CSV = csv, sql_path = NULL, output_dir = out))

  # outputs land in a dated per-run subdirectory: <region>_<label>_<Date>
  run <- file.path(out, paste0("testreg_testlab_", Sys.Date()))
  expect_true(dir.exists(run))
  expect_true(file.exists(file.path(run, "comprehensive_testreg_testlab_list_for_making_localdb.txt")))
  expect_true(file.exists(file.path(run, "comprehensive_testreg_testlab_list_for_LCA.txt")))
  expect_true(file.exists(file.path(run, "comprehensive_testreg_testlab_list_methods.txt")))
})

test_that("pre-made OBIS/GBIF source CSVs are fed in and stacked", {
  obis <- tempfile(fileext = ".csv"); gbif <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(Species = "Aaa bbb", Source = "OBIS"), obis, row.names = FALSE)
  utils::write.csv(data.frame(Species = "Ccc ddd", Source = "GBIF"), gbif, row.names = FALSE)

  res <- suppressMessages(build_regional_checklist(
    region = "r", label = "l", output_dir = tempfile(),
    OBIS = obis, GBIF = gbif, sql_path = NULL))
  expect_true(all(c("Aaa bbb", "Ccc ddd") %in% res$for_LCA$Species))
})

test_that("region/label are required and at least one source must be on", {
  expect_error(build_regional_checklist(label = "l", OBIS = FALSE, CSV = FALSE), "region")
  expect_error(build_regional_checklist(region = "r", OBIS = FALSE, CSV = FALSE), "label")
  expect_error(
    build_regional_checklist(region = "r", label = "l", output_dir = tempfile(),
                             OBIS = FALSE, GBIF = FALSE, CSV = FALSE),
    "No species")
})

test_that("a fresh download requires taxa and regional_poly", {
  expect_error(build_regional_checklist(region = "r", label = "l",
                                        output_dir = tempfile()),
               "taxa.*regional_poly|regional_poly")
})

test_that("for_LCA is taxID + ranks, checklist_summary holds the audit columns", {
  skip_on_cran()
  sql <- Sys.getenv("REGATTA_TAXONOMIZR_SQL", "accessionTaxa.sql")
  skip_if(!file.exists(sql), "taxonomizr DB not available")

  csv <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(Genus = "Lutjanus", Species = "argentiventris"),
                   csv, row.names = FALSE)
  res <- suppressWarnings(suppressMessages(build_regional_checklist(
    region = "r", label = "l", output_dir = tempfile(),
    OBIS = FALSE, CSV = csv, sql_path = sql)))   # existing DB -> used as-is

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
