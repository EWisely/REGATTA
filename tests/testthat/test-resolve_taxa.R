# Tests for resolve_taxa(). The resolution path queries WoRMS (and optionally
# GBIF), so the network-dependent tests skip on CRAN and when offline. The
# alias table itself is pure data and is tested without a network.

test_that("curated alias table has the expected shorthands", {
  al <- REGATTA:::.regatta_taxon_aliases()
  expect_type(al, "list")
  # one-to-many "fish" shorthand: ray-finned fish + sharks/rays + hagfish + lamprey
  expect_setequal(al[["fish"]],
                  c("Actinopterygii", "Elasmobranchii", "Myxini", "Petromyzonti"))
  expect_setequal(al[["vertebrates"]],
                  c("Actinopterygii", "Elasmobranchii", "Myxini", "Petromyzonti",
                    "Mammalia", "Aves", "Reptilia", "Amphibia"))
  expect_identical(al[["lepidosauria"]], "Reptilia")
})

test_that("kingdom filtering disambiguates Vertebrata to the vertebrate subphylum", {
  skip_on_cran()
  skip_if_offline()
  r <- resolve_taxa("Vertebrata", check_gbif = FALSE)
  expect_equal(nrow(r), 1L)
  expect_equal(r$aphia_id, 146419L)          # the Animalia subphylum, not the red-algae genus
  expect_equal(r$worms_status, "ok")
})

test_that("the AphiaID escape hatch resolves a numeric input directly", {
  skip_on_cran()
  skip_if_offline()
  r <- resolve_taxa("146419", check_gbif = FALSE)
  expect_equal(r$aphia_id, 146419L)
  expect_equal(tolower(r$valid_name), "vertebrata")
})

test_that("one-to-many 'fish' alias expands into four resolved taxa", {
  skip_on_cran()
  skip_if_offline()
  r <- resolve_taxa("fish", check_gbif = FALSE)
  expect_equal(nrow(r), 4L)
  expect_true(all(r$input == "fish"))        # original input preserved on every row
  expect_setequal(r$valid_name,
                  c("Actinopterygii", "Elasmobranchii", "Myxini", "Petromyzonti"))
  expect_true(all(r$worms_status == "ok"))
})

test_that("absent names are aliased (Lepidosauria -> Reptilia)", {
  skip_on_cran()
  skip_if_offline()
  r <- resolve_taxa("Lepidosauria", check_gbif = FALSE)
  expect_equal(r$alias_used, "Reptilia")
  expect_equal(tolower(r$valid_name), "reptilia")
  expect_equal(r$worms_status, "ok")
})

test_that("an unresolvable name errors under the default on_ambiguous = 'error'", {
  skip_on_cran()
  skip_if_offline()
  expect_error(resolve_taxa("Notarealtaxonxyz", check_gbif = FALSE),
               "could not resolve")
})

test_that("on_ambiguous = 'warn' returns the table instead of erroring", {
  skip_on_cran()
  skip_if_offline()
  expect_warning(
    r <- resolve_taxa("Notarealtaxonxyz", check_gbif = FALSE, on_ambiguous = "warn"),
    "could not resolve")
  expect_equal(r$worms_status, "not_found")
  expect_true(is.na(r$aphia_id))
})

test_that("GBIF coverage flags + descent: bony fish descend to keys, sharks use a direct key", {
  skip_on_cran()
  skip_if_offline()
  r <- resolve_taxa(c("Actinopterygii", "Elasmobranchii"), check_gbif = TRUE)
  expect_false(r$gbif_usable[r$valid_name == "Actinopterygii"])  # no usable class node
  expect_true(r$gbif_usable[r$valid_name == "Elasmobranchii"])   # class key 121
  expect_true(all(c("gbif_key", "gbif_usable", "gbif_keys", "note") %in% names(r)))

  # resolve_taxa now returns usable GBIF keys for both: Actinopterygii via
  # order/family descent, Elasmobranchii via its direct key.
  acti <- r$gbif_keys[[which(r$valid_name == "Actinopterygii")]]
  elas <- r$gbif_keys[[which(r$valid_name == "Elasmobranchii")]]
  expect_gt(length(acti), 0L)     # descended keys, not integer(0)
  expect_gt(length(elas), 0L)     # direct key
})
