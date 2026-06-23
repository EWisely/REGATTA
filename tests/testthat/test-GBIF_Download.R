
# Check that searching for Crustacea returns an error
test_that("1+1=2", {
  expect_equal(1+1, 2)
})

test_that("Searching GBIF for Crustacea returns an error", {
  expect_error(GBIF_download(obis_taxa = c("Crustacea"), 
                             regional_poly = "POLYGON ((-117.07 32.55, -117.81 34.27, -121.06 34.56, -120.84 31.69, -117.07 32.55))"))
})

test_that("Searching GBIF for Pancrustacea returns an error", {
  expect_error(GBIF_download(obis_taxa = c("Pancrustacea"), 
                             regional_poly = "POLYGON ((-117.07 32.55, -117.81 34.27, -121.06 34.56, -120.84 31.69, -117.07 32.55))"))
}) 