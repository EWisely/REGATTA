# Offline tests for the build_regional_checklist() orchestrator and the
# background/lazy taxonomize helper. here::here() is mocked so output lands in
# a temp dir; taxonomize_checklist() is mocked so no network or taxonomizr DB
# is ever touched.

ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

test_that("CSV source splits genus-only from species, and names files from region/label", {
  tmp <- tempfile(); dir.create(tmp)
  local_mocked_bindings(here = function(...) file.path(tmp, ...), .package = "here")

  csv <- file.path(tmp, "local.csv")
  utils::write.csv(
    data.frame(Genus   = c("Lutjanus", "Caranx", "Sebastes"),
               Species = c("argentiventris", "sp.", "mystinus")),
    csv, row.names = FALSE)

  res <- build_regional_checklist(
    region = "testreg", label = "testlab",
    OBIS = FALSE, GBIF = FALSE, CSV = csv,
    sql_path = file.path(tmp, "absent.sql"))   # no DB -> skip taxonomize

  expect_match(res$db_path,  "comprehensive_testreg_testlab_list_for_making_localdb\\.txt$")
  expect_match(res$lca_path, "comprehensive_testreg_testlab_list_for_LCA\\.txt$")

  db  <- utils::read.delim(res$db_path)
  lca <- utils::read.delim(res$lca_path)

  # DB list: species binomials only, NO bare genus
  expect_setequal(db$Species, c("Lutjanus argentiventris", "Sebastes mystinus"))
  expect_false("Caranx" %in% db$Species)

  # LCA list: same binomials PLUS the retained bare genus
  expect_true("Caranx" %in% lca$Species)
  expect_true(all(c("Lutjanus argentiventris", "Sebastes mystinus") %in% lca$Species))

  # No taxonomizr DB available -> no background taxonomize
  expect_null(res$rds_path)
})

test_that("a taxonomizr DB present triggers the background taxonomize", {
  skip_on_cran()
  sql <- Sys.getenv("REGATTA_TAXONOMIZR_SQL", "accessionTaxa.sql")
  skip_if(!file.exists(sql), "taxonomizr DB not available")

  tmp <- tempfile(); dir.create(tmp)
  local_mocked_bindings(here = function(...) file.path(tmp, ...), .package = "here")

  csv <- file.path(tmp, "local.csv")
  utils::write.csv(data.frame(Genus = "Lutjanus", Species = "argentiventris"),
                   csv, row.names = FALSE)

  res <- suppressMessages(
    build_regional_checklist(region = "r", label = "l",
                             OBIS = FALSE, CSV = csv, sql_path = sql))

  expect_match(res$rds_path, "_for_LCA_taxonomized\\.rds$")
  expect_true(file.exists(res$rds_path))
  tx <- readRDS(res$rds_path)                       # real taxonomized output
  expect_true(all(c("genus", "species") %in% names(tx)))
  expect_equal(tx$genus[tx$input_name == "Lutjanus argentiventris"], "Lutjanus")
})

test_that("pre-made OBIS/GBIF source CSVs are fed in and stacked", {
  tmp <- tempfile(); dir.create(tmp)
  local_mocked_bindings(here = function(...) file.path(tmp, ...), .package = "here")

  obis <- file.path(tmp, "obis.csv")
  gbif <- file.path(tmp, "gbif.csv")
  utils::write.csv(data.frame(Species = "Aaa bbb", Source = "OBIS"), obis, row.names = FALSE)
  utils::write.csv(data.frame(Species = "Ccc ddd", Source = "GBIF"), gbif, row.names = FALSE)

  res <- build_regional_checklist(region = "r", label = "l",
                                  OBIS = obis, GBIF = gbif,
                                  sql_path = file.path(tmp, "absent.sql"))
  lca <- utils::read.delim(res$lca_path)
  expect_true(all(c("Aaa bbb", "Ccc ddd") %in% lca$Species))
})

test_that("region/label are required and at least one source must be on", {
  tmp <- tempfile(); dir.create(tmp)
  local_mocked_bindings(here = function(...) file.path(tmp, ...), .package = "here")

  expect_error(build_regional_checklist(label = "l", OBIS = FALSE, CSV = FALSE), "region")
  expect_error(build_regional_checklist(region = "r", OBIS = FALSE, CSV = FALSE), "label")
  expect_error(
    build_regional_checklist(region = "r", label = "l",
                             OBIS = FALSE, GBIF = FALSE, CSV = FALSE),
    "No species")
})

test_that("a fresh download requires taxa and regional_poly", {
  tmp <- tempfile(); dir.create(tmp)
  local_mocked_bindings(here = function(...) file.path(tmp, ...), .package = "here")
  # OBIS = TRUE (default) but no taxa/poly -> clear error before any network call
  expect_error(build_regional_checklist(region = "r", label = "l"),
               "taxa.*regional_poly|regional_poly")
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
