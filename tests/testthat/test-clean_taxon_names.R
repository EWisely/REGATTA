# Pure, offline test for the internal name-cleaning helper used by
# resolve_names() / taxonomize_checklist().

test_that("clean_taxon_names strips abbreviations and quotes, word-boundary safe", {
  expect_equal(clean_taxon_names("Mugil sp."),          "Mugil")
  expect_equal(clean_taxon_names("Mugil spp."),         "Mugil")
  expect_equal(clean_taxon_names("Mugil cf. curema"),   "Mugil curema")
  expect_equal(clean_taxon_names("Mugil aff. curema"),  "Mugil curema")
  expect_equal(clean_taxon_names('"Mugil" curema'),     "Mugil curema")
  expect_equal(clean_taxon_names("Mugil curema"),       "Mugil curema")  # untouched
})
