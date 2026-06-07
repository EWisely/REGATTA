#' @keywords internal
"_PACKAGE"

## Imports used across the checklist-building helpers (GBIF_download,
## OBIS_download, build_regional_checklist). The heavy biodiversity packages
## (rgbif, robis, worrms, taxize) are Suggests and are called via `pkg::fn`
## behind requireNamespace() guards, so they are not imported here.
#' @importFrom dplyr %>% filter select rename mutate pull bind_rows
#' @importFrom readr write_delim
#' @importFrom utils read.csv read.delim write.csv head
NULL

## Column names referenced via dplyr non-standard evaluation inside the
## download/build helpers. Declaring them keeps R CMD check from flagging
## "no visible binding for global variable".
utils::globalVariables(c(
  "matchType", "usageKey", "verbatim_name", "verbatim_scientificName",
  "taxonRank", "scientificName", "species", "genus", "family", "order",
  "class", "Species", "Genus", "Source", "taxa",
  "is_marine", "is_brackish", "is_freshwater", "is_terrestrial"
))
