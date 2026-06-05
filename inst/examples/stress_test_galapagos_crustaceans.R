# stress_test_galapagos_crustaceans.R
# Eldridge Wisely
#
# Stress test: run the full REGATTA pipeline on the Galapagos
# crustaceans dataset (BerryCrust primer). This is a different
# taxonomic group, a different primer, and a much larger checklist
# than the MiFish validation case. It exercises:
#
#   - taxonomize_checklist() on a 4535-species crustacean list
#     (broad arthropod taxonomy, not just Chordata)
#   - resolve_taxids() / parse_vsearch_userout() on different ASV scales
#   - reconcile_global_local() / reconcile_checklist() on an off-target
#     -prone primer where %ID-based reconciliation matters more

source("taxonomize_checklist.R")
source("resolve_taxids.R")
source("parse_vsearch_userout.R")
source("reconcile_global_local.R")
source("reconcile_checklist.R")
source("summarize_regatta.R")

sql <- "/Users/Eldridge/taxonomy_download/accessionTaxa.sql"

# ----- Step 1: taxonomize the crustacean checklist (or load cached) -----
cached <- "custom_db/comprehensive_galapagos_crustaceans_list_taxonomized.rds"
if (file.exists(cached)) {
  cat("Loading cached crustacean checklist:", cached, "\n")
  crust_checklist <- readRDS(cached)
} else {
  cat("Taxonomizing crustacean checklist (first run; ~30-60s for 4535 names)\n")
  crust_checklist <- taxonomize_checklist(
    input    = "custom_db/comprehensive_galapagos_crustaceans_list.txt",
    sql_path = sql
  )
  saveRDS(crust_checklist, cached)
}
cat("Checklist:", nrow(crust_checklist), "rows; resolved:",
    sum(crust_checklist$resolution_status == "resolved"),
    sprintf("(%.1f%%)\n",
            100 * sum(crust_checklist$resolution_status == "resolved") / nrow(crust_checklist)))
cat("Synonym matches:",
    sum(crust_checklist$name_match_type == "synonym", na.rm = TRUE), "\n")

# ----- Step 2: prepare BerryCrust classifier inputs -----
cat("\n--- Reading BerryCrust obitools (.tab, ~12 MB) ---\n")
t0 <- Sys.time()
obi_raw <- utils::read.delim("script_03_inputs/BerryCrust_Menu_95_named.tab",
                              sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
cat("  Read in", format(Sys.time() - t0), ";", nrow(obi_raw), "ASVs,",
    ncol(obi_raw), "cols\n")
cat("  Global pct_id (BEST_IDENTITY) range: ",
    paste(range(suppressWarnings(as.numeric(obi_raw$BEST_IDENTITY)), na.rm = TRUE),
          collapse = "-"), "\n")

cat("\n--- Resolving obitools TAXIDs ---\n")
t0 <- Sys.time()
global_in <- resolve_taxids(obi_raw, taxid_col = "TAXID", sql_path = sql)
global_in$pct_id <- suppressWarnings(as.numeric(obi_raw$BEST_IDENTITY))  # 0-1, auto-rescaled
names(global_in)[names(global_in) == "ID"] <- "ASV_id"
cat("  Resolved in", format(Sys.time() - t0), ";",
    sum(!is.na(global_in$domain)), "of", nrow(global_in), "ASVs got taxonomy\n")

cat("\n--- Parsing vsearch BerryCrust userout ---\n")
local_in <- parse_vsearch_userout(
  "script_03_inputs/userout_crustaceans_Galapagos_top5_comprehensive_galapagos_results.txt"
)
cat("  local_in:", nrow(local_in), "ASVs assigned by vsearch+local DB\n")

# ----- Step 3: run the REGATTA pipeline -----
cat("\n--- reconcile_global_local ---\n")
t0 <- Sys.time()
rec <- reconcile_global_local(global_in, local_in,
         output_dir    = "stress_crustaceans_out/reconcile_global_local",
         output_prefix = "berrycrust_reconciled")
cat("  Done in", format(Sys.time() - t0), "\n")
print(rec$stats)

cat("\n--- reconcile_checklist (against crustacean regional checklist) ---\n")
t0 <- Sys.time()
post <- reconcile_checklist(rec$result, crust_checklist, id_col = "ASV_id",
                            output_dir   = "stress_crustaceans_out/reconcile_checklist",
                            output_prefix = "berrycrust_regatta",
                            prior_dir    = "stress_crustaceans_out/reconcile_global_local",
                            prior_prefix = "berrycrust_reconciled")
cat("  Done in", format(Sys.time() - t0), "\n")
print(post$stats)

# ----- Step 4: cross-stage summary -----
cat("\n--- summarize_regatta ---\n")
summ <- summarize_regatta(rec, post, global_in, local_in)
print(summ)
write.csv(summ, "stress_crustaceans_out/berrycrust_regatta_summary.csv", row.names = FALSE)
cat("\nWrote summary to stress_crustaceans_out/berrycrust_regatta_summary.csv\n")
