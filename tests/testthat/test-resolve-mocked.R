# Mocked tests for the network/DB-backed resolvers, so their logic runs in CI
# without WoRMS, GBIF, or a taxonomizr database. worrms/rgbif and
# taxonomizr/name_to_taxid are mocked with local_mocked_bindings().

ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

# ---- resolve_taxa() -------------------------------------------------------

test_that("resolve_taxa disambiguates Vertebrata to the Animalia subphylum by kingdom", {
  local_mocked_bindings(
    wm_records_name = function(name, marine_only = FALSE, ...) {
      data.frame(AphiaID        = c(146419L, 370319L),
                 scientificname = c("Vertebrata", "Vertebrata"),
                 rank           = c("Subphylum", "Genus"),
                 kingdom        = c("Animalia", "Plantae"),
                 status         = c("accepted", "accepted"),
                 stringsAsFactors = FALSE)
    }, .package = "worrms")
  r <- resolve_taxa("Vertebrata", kingdom = "Animalia", check_gbif = FALSE)
  expect_equal(nrow(r), 1L)
  expect_equal(r$aphia_id, 146419L)
  expect_equal(r$worms_status, "ok")
})

test_that("resolve_taxa expands the one-to-many 'fish' alias into four rows", {
  local_mocked_bindings(
    wm_records_name = function(name, marine_only = FALSE, ...) {
      data.frame(AphiaID = 100L, scientificname = name, rank = "Class",
                 kingdom = "Animalia", status = "accepted", stringsAsFactors = FALSE)
    }, .package = "worrms")
  r <- resolve_taxa("fish", check_gbif = FALSE)
  expect_equal(nrow(r), 4L)
  expect_true(all(r$input == "fish"))
  expect_setequal(r$valid_name,
                  c("Actinopterygii", "Elasmobranchii", "Myxini", "Petromyzonti"))
  expect_true(all(r$worms_status == "ok"))
})

test_that("resolve_taxa errors (default) or warns (on_ambiguous='warn') on an absent name", {
  local_mocked_bindings(
    wm_records_name = function(name, marine_only = FALSE, ...) stop("204 No Content"),
    .package = "worrms")
  expect_error(resolve_taxa("Notarealtaxonxyz", check_gbif = FALSE), "could not resolve")
  expect_warning(
    r <- resolve_taxa("Notarealtaxonxyz", check_gbif = FALSE, on_ambiguous = "warn"),
    "could not resolve")
  expect_equal(r$worms_status, "not_found")
  expect_true(is.na(r$aphia_id))
})

test_that("resolve_taxa's GBIF diagnostic flags usable vs unusable backbone keys", {
  local_mocked_bindings(
    wm_records_name = function(name, marine_only = FALSE, ...) {
      data.frame(AphiaID = 100L, scientificname = name, rank = "Class",
                 kingdom = "Animalia", status = "accepted", stringsAsFactors = FALSE)
    }, .package = "worrms")
  local_mocked_bindings(
    name_backbone = function(name, rank = NULL, kingdom = NULL, strict = FALSE, ...) {
      if (name == "Mammalia")
        data.frame(usageKey = 359L, matchType = "EXACT",
                   canonicalName = "Mammalia", scientificName = "Mammalia",
                   stringsAsFactors = FALSE)
      else
        data.frame(usageKey = NA_integer_, matchType = "NONE",
                   canonicalName = NA_character_, scientificName = NA_character_,
                   stringsAsFactors = FALSE)
    }, .package = "rgbif")
  r <- resolve_taxa(c("Mammalia", "Actinopterygii"), check_gbif = TRUE)
  expect_true(r$gbif_usable[r$valid_name == "Mammalia"])
  expect_equal(r$gbif_key[r$valid_name == "Mammalia"], 359L)
  expect_false(r$gbif_usable[r$valid_name == "Actinopterygii"])
})

# ---- resolve_names() ------------------------------------------------------

test_that("clean_taxon_names strips abbreviations/quotes and is word-boundary safe", {
  expect_equal(clean_taxon_names("Sebastes sp."), "Sebastes")
  expect_equal(clean_taxon_names("Sebastes spp."), "Sebastes")
  expect_equal(clean_taxon_names("Sebastes cf. mystinus"), "Sebastes mystinus")
  expect_equal(clean_taxon_names('"Sebastes" mystinus'), "Sebastes mystinus")
  expect_equal(clean_taxon_names("Sebastes"), "Sebastes")   # not corrupted to "ebastes"
})

test_that("resolve_names attaches ranks + match type and leaves unresolved rows NA", {
  tmp <- tempfile(); dir.create(tmp); sql <- file.path(tmp, "db.sql"); file.create(sql)
  local_mocked_bindings(name_to_taxid = function(names, sql_path, accept_types) {
    data.frame(taxID = c(331610L, NA), match_type = c("scientific name", NA),
               stringsAsFactors = FALSE)
  })
  local_mocked_bindings(getTaxonomy = function(ids, sqlFile, desiredTaxa, ...) {
    m <- matrix(NA_character_, nrow = length(ids), ncol = length(desiredTaxa),
                dimnames = list(NULL, desiredTaxa))
    m[1, "genus"] <- "Sebastes"; m[1, "species"] <- "Sebastes mystinus"; m
  }, .package = "taxonomizr")

  out <- suppressMessages(resolve_names(c("Sebastes mystinus", "Notarealtaxon"), sql_path = sql))
  expect_true(all(c("input_name", "name_match_type", "genus", "species") %in% names(out)))
  expect_equal(out$species[1], "Sebastes mystinus")
  expect_true(is.na(out$species[2]))
  expect_equal(out$name_match_type, c("scientific name", NA))
})

# ---- resolve_taxids() -----------------------------------------------------

test_that("resolve_taxids maps a taxID vector to 7 ranks", {
  tmp <- tempfile(); dir.create(tmp); sql <- file.path(tmp, "db.sql"); file.create(sql)
  local_mocked_bindings(getTaxonomy = function(ids, sqlFile, desiredTaxa, ...) {
    m <- matrix(NA_character_, nrow = length(ids), ncol = length(desiredTaxa),
                dimnames = list(NULL, desiredTaxa))
    m[, "genus"]   <- "Sebastes"
    m[, "species"] <- "Sebastes mystinus"
    m
  }, .package = "taxonomizr")

  out <- resolve_taxids(c(331610, 331610), sql_path = sql)
  expect_true(all(c("taxID", ranks) %in% names(out)))
  expect_equal(out$species, c("Sebastes mystinus", "Sebastes mystinus"))
})

test_that("resolve_taxids errors when the SQL DB is absent", {
  expect_error(resolve_taxids(123, sql_path = "/no/such/accessionTaxa.sql"),
               "SQL DB not found")
})
