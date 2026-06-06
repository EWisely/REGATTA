# stress_test_california_current_vertebrates.R
# Eldridge Wisely
#
# Full naive-user stress test of REGATTA on a different region + group
# combo: California Current marine vertebrates against the MOSAIC2
# RL2501 MV1 MURI BestTaxon classifier output.
#
# Inputs (the only two REGATTA needs upstream):
#   - WKT polygon (drawn on wktmap.com)
#   - MURI BestTaxon-style classifier output (CSV with BestTaxon column)
#
# Pipeline:
#   1. OBIS_download + GBIF_download for marine vertebrate classes
#   2. build_regional_checklist + taxonomize_checklist
#   3. run_regatta() on the MOSAIC2 file (auto-detects BestTaxon path)
#   4. Inspect outputs

# Source all R/ files (running as a script, not as installed package)
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

# ---- Inputs -----------------------------------------------------------
ca_current_poly <- "POLYGON ((-124.453125 41.967659, -137.373047 41.771312, -133.59375 29.305561, -129.726563 24.527135, -113.818359 28.07198, -124.453125 41.967659))"

# California Current marine vertebrate classes. MURI primer targets
# vertebrates broadly: ray-finned fish + cartilaginous fish + mammals +
# birds (seabirds). All grouped here because the primer intentionally
# amplifies all of them; the off-target check still flags any
# non-vertebrate amplification (invertebrate ASVs, contaminants, etc.).
vert_taxa  <- c("Actinopterygii", "Chondrichthyes", "Mammalia", "Aves")

muri_input <- "/Users/Eldridge/ClaudeCode/research/MOSAIC2/new_lab_protocol_processed_data/results_run2_RL2501/output/run2_RL2501_MV1_taxon_table.csv"

sql_path   <- "/Users/Eldridge/taxonomy_download/accessionTaxa.sql"

# ---- Step 1: OBIS_download --------------------------------------------
# Fast (no auth, smaller queries). Marine subset only — that's what we
# want for a marine vertebrate study.
cat("==== OBIS_download ====\n")
OBIS_download(
  obis_taxa       = vert_taxa,
  regional_poly   = ca_current_poly,
  obis_outputname = "OBIS_CalCurrent_vertebrates",
  marine          = TRUE,
  freshwater      = NA,   # keep species that span habitats
  terrestrial     = FALSE,
  brackish        = NA
)

# ---- Step 2: GBIF_download --------------------------------------------
# Slow (~15 min). Will print "RUNNING" status while polling the
# download API. Credentials must already be in ~/.Renviron.
cat("\n==== GBIF_download ====\n")
GBIF_download(
  obis_taxa       = vert_taxa,
  worms_taxa      = vert_taxa,   # same names; WoRMS will return marine subset
  regional_poly   = ca_current_poly,
  gbif_outputname = "GBIF_CalCurrent_vertebrates"
)

# ---- Step 3: build_regional_checklist + taxonomize_checklist ----------
cat("\n==== build_regional_checklist ====\n")
build_regional_checklist(
  comb_inputnames = c("OBIS_CalCurrent_vertebrates",
                      "GBIF_CalCurrent_vertebrates"),
  comb_outputname = "Comprehensive_CalCurrent_vertebrates_list"
)

cat("\n==== taxonomize_checklist ====\n")
calcurrent_checklist <- taxonomize_checklist(
  input    = "custom_db/Comprehensive_CalCurrent_vertebrates_list.txt",
  sql_path = sql_path
)
saveRDS(calcurrent_checklist,
        "custom_db/Comprehensive_CalCurrent_vertebrates_list_taxonomized.rds")
cat("  Checklist rows:", nrow(calcurrent_checklist), "\n")
cat("  Resolved:",       sum(calcurrent_checklist$resolution_status == "resolved"), "\n")
cat("  Synonym matches:", sum(calcurrent_checklist$name_match_type == "synonym", na.rm = TRUE), "\n")

# ---- Step 4: run_regatta on MOSAIC2 -----------------------------------
cat("\n==== run_regatta on MOSAIC2 MURI BestTaxon ====\n")
res <- run_regatta(
  input     = muri_input,
  checklist = calcurrent_checklist,
  out_dir   = "stress_calcurrent_vertebrates_out",
  sql_path  = sql_path
)

cat("\n==== summary ====\n")
print(res$summary)
