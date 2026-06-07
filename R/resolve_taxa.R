# resolve_taxa.R
# Eldridge Wisely
#
# Validate and disambiguate the high-level query taxa that GBIF_download()
# and OBIS_download() pull species lists for. This is the front-end fix for
# the WoRMS/GBIF name-resolution traps that otherwise make the downloaders
# crash or, worse, silently fetch the wrong taxon:
#
#   * Ambiguous names. "Vertebrata" is both the vertebrate subphylum
#     (Animalia) and a red-algae genus (Plantae). worrms::wm_name2id()
#     errors (206 Partial Content); robis::checklist("Vertebrata") quietly
#     returns the seaweed. Filtering WoRMS matches by `kingdom` resolves
#     this: within Animalia only the subphylum remains, and we hand the
#     downloaders the AphiaID (not the ambiguous name).
#   * Wrong-kingdom resolution. A name whose only accepted match is outside
#     the target kingdom is flagged rather than fetched.
#   * Absent names. "Lepidosauria" is not in WoRMS (204 No Content); the
#     alias table maps it to "Reptilia", and unknown names error clearly.
#   * GBIF backbone gaps. GBIF's matcher returns no key for broad ray-finned
#     fish names (Actinopterygii / Actinopteri / Teleostei all NONE), so GBIF
#     cannot return bony fish from a class name. resolve_taxa() reports this
#     per taxon (`gbif_usable`) so the gap is explicit, not a silent dropout.

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a[1])) b else a

#' Call a WoRMS lookup with retry, distinguishing absent (204) from transient
#' errors (throttling, network). Returns the result, the string "not_found"
#' for a genuine 204, or a condition object for a transient failure.
#' @keywords internal
#' @noRd
.worms_try <- function(fn, tries = 3L) {
  for (i in seq_len(tries)) {
    res <- tryCatch(fn(), error = function(e) e)
    if (!inherits(res, "condition")) return(res)
    if (grepl("204", conditionMessage(res))) return("not_found")  # genuinely absent
    if (i < tries) Sys.sleep(1)                                    # transient: back off
  }
  res
}

#' Curated WoRMS/GBIF name aliases (internal)
#'
#' Maps lower-cased input names to one OR MORE substitute names used for the
#' WoRMS lookup. A one-to-many alias (e.g. `"fish"`) expands into several
#' resolved taxa. Keeps the known REGATTA traps and convenient shorthands in
#' one place. Users can extend or replace the table via the `aliases` argument
#' of [resolve_taxa()].
#' @keywords internal
#' @noRd
.regatta_taxon_aliases <- function() {
  list(
    "lepidosauria" = "Reptilia",                    # not a WoRMS name
    "osteichthyes" = "Actinopteri",                 # WoRMS class for ray-finned fish
    # ray-finned fish + sharks/rays + hagfishes + lampreys (the fishes a
    # vertebrate metabarcoding primer like MiFish targets):
    "fish"         = c("Actinopterygii", "Elasmobranchii", "Myxini", "Petromyzonti"),
    # All vertebrate classes. Expands like "fish" (rather than mapping to the
    # subphylum "Vertebrata", which has no usable GBIF backbone key) so it works
    # in BOTH OBIS (queried by AphiaID) and GBIF (direct keys for
    # mammals/birds/etc. + order descent for ray-finned fish and reptiles):
    "vertebrates"  = c("Actinopterygii", "Elasmobranchii", "Myxini", "Petromyzonti",
                       "Mammalia", "Aves", "Reptilia", "Amphibia")
  )
}

#' One-name WoRMS resolution (internal)
#' @keywords internal
#' @noRd
.resolve_one_worms <- function(name, kingdom) {
  nf  <- function() list(status = "not_found", aphia = NA_integer_,
                         valid = NA_character_, rank = NA_character_,
                         kingdom = NA_character_, n = 0L, others = character(0))
  err <- function(res) list(status = "lookup_error", aphia = NA_integer_,
                            valid = NA_character_, rank = NA_character_,
                            kingdom = NA_character_, n = 0L,
                            others = conditionMessage(res))

  # Numeric input is treated as an explicit AphiaID escape hatch.
  if (grepl("^[0-9]+$", trimws(name))) {
    rec <- .worms_try(function() worrms::wm_record(as.integer(name)))
    if (identical(rec, "not_found")) return(nf())
    if (inherits(rec, "condition"))  return(err(rec))
    if (is.null(rec) || nrow(rec) == 0) return(nf())
    return(list(status = "ok", aphia = as.integer(rec$AphiaID[1]),
                valid = rec$scientificname[1], rank = rec$rank[1],
                kingdom = rec$kingdom[1], n = 1L, others = character(0)))
  }

  recs <- .worms_try(function() worrms::wm_records_name(name, marine_only = FALSE))
  if (identical(recs, "not_found")) return(nf())
  if (inherits(recs, "condition"))  return(err(recs))
  if (is.null(recs) || nrow(recs) == 0) return(nf())
  # wm_records_name() returns partial matches too (e.g. "Vertebrata" pulls in
  # every "Vertebrata <species>" binomial). Keep only exact name matches.
  recs <- recs[!is.na(recs$scientificname) &
                 tolower(recs$scientificname) == tolower(name), , drop = FALSE]
  if (nrow(recs) == 0) {
    return(list(status = "not_found", aphia = NA_integer_, valid = NA_character_,
                rank = NA_character_, kingdom = NA_character_, n = 0L,
                others = character(0)))
  }

  fmt <- function(df) {
    s <- paste0(df$AphiaID, " (", df$scientificname, ", ", df$rank, ", ", df$kingdom, ")")
    if (length(s) > 8L) c(s[1:8], paste0("... and ", length(s) - 8L, " more")) else s
  }

  acc  <- recs[!is.na(recs$status) & recs$status == "accepted", , drop = FALSE]
  pool <- if (is.null(kingdom)) acc
          else acc[!is.na(acc$kingdom) & acc$kingdom == kingdom, , drop = FALSE]

  if (nrow(pool) == 1L) {
    return(list(status = "ok", aphia = as.integer(pool$AphiaID[1]),
                valid = pool$scientificname[1], rank = pool$rank[1],
                kingdom = pool$kingdom[1], n = 1L, others = character(0)))
  }
  if (nrow(pool) > 1L) {
    return(list(status = "ambiguous", aphia = NA_integer_, valid = NA_character_,
                rank = NA_character_, kingdom = NA_character_,
                n = nrow(pool), others = fmt(pool)))
  }
  # No accepted match in the target kingdom.
  if (nrow(acc) > 0L) {
    return(list(status = "wrong_kingdom", aphia = NA_integer_, valid = NA_character_,
                rank = NA_character_, kingdom = NA_character_,
                n = nrow(acc), others = fmt(acc)))
  }
  list(status = "not_found", aphia = NA_integer_, valid = NA_character_,
       rank = NA_character_, kingdom = NA_character_, n = 0L,
       others = fmt(recs))
}

#' GBIF backbone diagnostic for one resolved taxon (internal)
#' @keywords internal
#' @noRd
.gbif_diag <- function(valid_name, rank, kingdom) {
  if (!requireNamespace("rgbif", quietly = TRUE) || is.na(valid_name)) {
    return(list(key = NA_integer_, usable = NA))
  }
  bb <- tryCatch(
    rgbif::name_backbone(name = valid_name, rank = tolower(rank %||% NA),
                         kingdom = kingdom, strict = FALSE),
    error = function(e) NULL)
  if (is.null(bb) || is.null(bb$usageKey) || identical(bb$matchType, "NONE")) {
    return(list(key = NA_integer_, usable = FALSE))
  }
  canon  <- bb$canonicalName %||% bb$scientificName
  # Usable only if GBIF matched the taxon itself (exact/fuzzy), or returned a
  # higher-rank node that is still the SAME name. A HIGHERRANK bail-up to a
  # different name (e.g. Reptilia -> Chordata, Vertebrata -> Animalia) is not
  # usable: it would over-broaden the occurrence download.
  usable <- bb$matchType %in% c("EXACT", "FUZZY") ||
    (identical(bb$matchType, "HIGHERRANK") &&
       !is.na(canon) && tolower(canon) == tolower(valid_name))
  list(key = if (isTRUE(usable)) as.integer(bb$usageKey) else NA_integer_,
       usable = isTRUE(usable))
}

#' Usable GBIF backbone key(s) for one taxon (internal)
#'
#' Returns the taxon's own backbone key when it is usable; otherwise **walks
#' down** the WoRMS tree (to `descend_to`, with an optional family gap-fill) and
#' returns the GBIF keys those descendants resolve to. This is what makes broad
#' nodes GBIF lacks a usable key for -- ray-finned fish, `Crustacea`,
#' `Vertebrata`, ... -- still return occurrences. OBIS needs no equivalent: a
#' query by the parent AphiaID returns all descendants server-side.
#' @keywords internal
#' @noRd
.gbif_keys_for <- function(aphia_id, direct_key, direct_usable,
                           descend_to = "order", fill_families = TRUE) {
  if (isTRUE(direct_usable)) return(stats::na.omit(as.integer(direct_key)))
  if (!requireNamespace("rgbif", quietly = TRUE) ||
      !requireNamespace("taxize", quietly = TRUE) ||
      is.na(aphia_id)) {
    return(integer(0))
  }

  walk <- function(rank) {
    d <- tryCatch(taxize::worms_downstream(id = aphia_id, downto = rank),
                  error = function(e) NULL)
    if (!is.null(d) && nrow(d) > 0) unique(d$name) else character(0)
  }
  # GBIF backbone match with retry/backoff -- name_backbone_checklist() can
  # return a transient rate-limit error ("Status: 0 - try lower bucket_size or
  # larger sleep"). Returns NULL (treated as "no matches") if it keeps failing.
  nb <- function(nms) {
    for (i in seq_len(3L)) {
      res <- tryCatch(rgbif::name_backbone_checklist(nms), error = function(e) e)
      if (!inherits(res, "condition")) return(res)
      if (i < 3L) Sys.sleep(2L * i)
    }
    NULL
  }

  # Primary descent (default: order -- fast, hits a rank GBIF populates).
  primary_names <- walk(descend_to)
  primary_keys  <- integer(0)
  valid_phyla   <- character(0)
  if (length(primary_names) > 0) {
    bbp <- nb(primary_names)
    if (!is.null(bbp)) {
      keep <- bbp$matchType != "NONE"
      primary_keys <- bbp$usageKey[keep]
      if ("phylum" %in% names(bbp)) valid_phyla <- unique(stats::na.omit(bbp$phylum[keep]))
    }
  }

  # Gap-fill: families whose GBIF parent at the primary rank isn't in the
  # matched set (e.g. families GBIF files with no order). Phylum-guarded so a
  # stray family-name collision can't drag in off-target taxa.
  fill_keys <- integer(0)
  if (isTRUE(fill_families) && descend_to != "family") {
    fam_names <- walk("family")
    bbf <- if (length(fam_names) > 0) nb(fam_names) else NULL
    if (!is.null(bbf)) {
      bbf <- bbf[bbf$matchType != "NONE", , drop = FALSE]
      parent_col <- paste0(descend_to, "Key")
      covered    <- if (parent_col %in% names(bbf)) bbf[[parent_col]] %in% primary_keys else FALSE
      covered[is.na(covered)] <- FALSE
      in_phylum  <- if (length(valid_phyla) && "phylum" %in% names(bbf))
                      bbf$phylum %in% valid_phyla else TRUE
      fill_keys  <- bbf$usageKey[!covered & in_phylum]
    }
  }

  unique(as.integer(c(primary_keys, fill_keys)))
}

#' Validate and disambiguate query taxa for the checklist downloaders
#'
#' Resolves a vector of high-level taxon names to unambiguous WoRMS AphiaIDs,
#' filtering by `kingdom` so homonyms across kingdoms (e.g. *Vertebrata* the
#' vertebrate subphylum vs. *Vertebrata* the red-algae genus) collapse to the
#' intended taxon. Reports, per taxon, whether GBIF's backbone has a usable
#' key, so the well-known GBIF gap for ray-finned fish is explicit rather than
#' a silent dropout. [GBIF_download()] and [OBIS_download()] call this
#' internally; you can also call it directly to sanity-check your group names
#' before launching a long download.
#'
#' @param taxa Character vector of taxon names. A purely numeric element is
#'   treated as an explicit WoRMS AphiaID (the escape hatch for a genuinely
#'   ambiguous name -- look it up once, then pass the id).
#' @param kingdom Target kingdom used to disambiguate WoRMS matches. Default
#'   `"Animalia"`. Pass `NULL` to disable kingdom filtering.
#' @param check_gbif If TRUE (default when `rgbif` is installed), add the GBIF
#'   backbone diagnostic columns.
#' @param on_ambiguous What to do when a name is ambiguous, not found, or has
#'   no accepted match in `kingdom`: `"error"` (default) stops with a report
#'   listing the offending names and their candidate AphiaIDs; `"warn"` keeps
#'   going with those rows left unresolved; `"first"` takes the first candidate
#'   and warns.
#' @param aliases A named list (or named character vector) mapping (lower-cased)
#'   input names to one or more substitute names for the WoRMS lookup. Defaults
#'   to REGATTA's curated table, which includes one-to-many shorthands:
#'   `"fish"` -> ray-finned fish + sharks/rays + hagfishes + lampreys
#'   (`Actinopterygii`, `Elasmobranchii`, `Myxini`, `Petromyzonti`),
#'   `"vertebrates"` -> those plus `Mammalia`, `Aves`, `Reptilia`, `Amphibia`
#'   (the vertebrate classes; works in both OBIS and GBIF, unlike the bare
#'   subphylum `Vertebrata`), `"lepidosauria"` -> `"Reptilia"`,
#'   `"osteichthyes"` -> `"Actinopteri"`. Pass your own to extend or override. A
#'   multi-name alias expands into one result row per taxon.
#' @param gbif_descend_to,gbif_fill_families Control the GBIF descent (used only
#'   when `check_gbif` and a taxon has no usable direct key). `gbif_descend_to`
#'   is the WoRMS rank to walk down to (default `"order"` -- fast, and a rank
#'   GBIF populates); `gbif_fill_families` (default TRUE) additionally adds keys
#'   for families whose GBIF parent at that rank isn't matched, lifting
#'   coverage.
#'
#' @return A data.frame with one row per resolved taxon: `input` (the original
#'   name you passed -- repeated when an alias expanded to several taxa),
#'   `alias_used`, `aphia_id` (the **OBIS** query key -- OBIS descends it
#'   server-side), `valid_name`, `rank`, `kingdom`, `worms_status`
#'   (`ok`/`ambiguous`/`wrong_kingdom`/`not_found`), `n_candidates`, and (when
#'   `check_gbif`) `gbif_key` (the taxon's own backbone key, NA if not usable),
#'   `gbif_usable`, `gbif_keys` (a list-column: the **GBIF** query keys -- the
#'   direct key when usable, otherwise the keys reached by descent), plus a
#'   human-readable `note`.
#'
#' @examples
#' \dontrun{
#' # The well-known trap resolves correctly within Animalia:
#' resolve_taxa("Vertebrata")            # -> AphiaID 146419, gbif_usable FALSE
#'
#' # One-to-many shorthand: "fish" expands to ray-finned fish + sharks/rays.
#' resolve_taxa("fish")                  # -> two rows: Actinopterygii, Elasmobranchii
#'
#' # Pre-check a class list before a long download:
#' resolve_taxa(c("Actinopteri", "Chondrichthyes", "Mammalia",
#'                "Aves", "Lepidosauria"))   # Lepidosauria -> Reptilia
#' }
#'
#' @export
resolve_taxa <- function(taxa,
                         kingdom            = "Animalia",
                         check_gbif         = requireNamespace("rgbif", quietly = TRUE),
                         on_ambiguous       = c("error", "warn", "first"),
                         aliases            = .regatta_taxon_aliases(),
                         gbif_descend_to    = "order",
                         gbif_fill_families = TRUE) {
  if (!requireNamespace("worrms", quietly = TRUE)) {
    stop("Package 'worrms' is required for resolve_taxa().")
  }
  on_ambiguous <- match.arg(on_ambiguous)
  taxa <- as.character(taxa)
  if (!is.null(aliases) && !is.list(aliases)) aliases <- as.list(aliases)

  # Expand aliases first: an input may map to one OR MORE lookup names
  # (e.g. "fish" -> Actinopteri + Elasmobranchii). Each becomes its own row,
  # tagged with the original `input`.
  jobs <- list()
  for (nm in taxa) {
    key <- tolower(trimws(nm))
    if (!is.null(aliases) && key %in% names(aliases)) {
      for (sub in aliases[[key]]) jobs[[length(jobs) + 1L]] <-
        list(input = nm, lookup = sub, alias = sub)
    } else {
      jobs[[length(jobs) + 1L]] <- list(input = nm, lookup = nm, alias = NA_character_)
    }
  }

  rows <- lapply(jobs, function(job) {
    nm    <- job$input
    alias <- job$alias

    w <- .resolve_one_worms(job$lookup, kingdom)

    if (identical(on_ambiguous, "first") && w$status == "ambiguous") {
      # Re-resolve taking the first candidate's AphiaID.
      first_id <- sub(" .*$", "", w$others[1])
      w <- .resolve_one_worms(first_id, kingdom)
      w$status <- "ok"
    }

    g <- if (isTRUE(check_gbif) && w$status == "ok") {
      .gbif_diag(w$valid, w$rank, w$kingdom)
    } else list(key = NA_integer_, usable = NA)

    # Usable GBIF backbone keys for the query: the direct key when usable,
    # otherwise the keys reached by walking WoRMS down to orders/families.
    gbif_keys_vec <- if (isTRUE(check_gbif) && w$status == "ok") {
      .gbif_keys_for(w$aphia, g$key, g$usable, gbif_descend_to, gbif_fill_families)
    } else integer(0)

    note <- switch(
      w$status,
      ok            = if (!isTRUE(check_gbif) || isTRUE(g$usable)) "resolved"
                      else if (length(gbif_keys_vec) > 0)
                        paste0("resolved; GBIF via ", gbif_descend_to, " descent (",
                               length(gbif_keys_vec), " key(s))")
                      else "resolved; no usable GBIF backbone key (use OBIS for this taxon)",
      ambiguous     = paste0(w$n, " accepted matches in kingdom; pass an AphiaID. Candidates: ",
                             paste(w$others, collapse = "; ")),
      wrong_kingdom = paste0("no accepted match in kingdom '", kingdom,
                             "'; found instead: ", paste(w$others, collapse = "; ")),
      not_found     = if (!is.na(alias))
                        paste0("not found in WoRMS (alias '", alias, "' also failed)")
                      else "not found in WoRMS",
      lookup_error  = paste0("WoRMS lookup failed (transient -- retry): ",
                             paste(w$others, collapse = " "))
    )

    data.frame(
      input        = nm,
      alias_used   = alias,
      aphia_id     = w$aphia,
      valid_name   = w$valid,
      rank         = w$rank,
      kingdom      = w$kingdom,
      worms_status = w$status,
      n_candidates = w$n,
      gbif_key     = g$key,
      gbif_usable  = g$usable,
      gbif_keys    = I(list(gbif_keys_vec)),
      note         = note,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  bad <- out[out$worms_status != "ok", , drop = FALSE]
  if (nrow(bad) > 0L) {
    msg <- paste0("resolve_taxa() could not resolve ", nrow(bad), " of ",
                  nrow(out), " taxon name(s):\n",
                  paste0("  - ", bad$input, ": ", bad$note, collapse = "\n"))
    if (identical(on_ambiguous, "error")) stop(msg, call. = FALSE)
    warning(msg, call. = FALSE)
  }
  out
}
