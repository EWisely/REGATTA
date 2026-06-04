# end_to_end_raw_galapagos.R
# Eldridge Wisely
#
# Worked example: run the full REGATTA pipeline starting from raw
# obitools (.tab) + vsearch (--userout) files for the Galapagos MiFish
# dataset. Uses only public REGATTA functions — no manual data
# munging — and produces the standard 3-CSV-per-stage output triples
# in reconcile_global_local_out/ and reconcile_checklist_out/, plus
# the 21-row Ella-format summary via summarize_regatta().

source("resolve_taxids.R")
source("parse_vsearch_userout.R")
source("reconcile_global_local.R")
source("reconcile_checklist.R")
source("summarize_regatta.R")

sql            <- "/Users/Eldridge/taxonomy_download/accessionTaxa.sql"
fish_checklist <- readRDS("custom_db/comprehensive_galapagos_fish_list_taxonomized.rds")

# ----- Global DB side: obitools .tab file -----
obi_raw   <- utils::read.delim("script_03_inputs/MiFish_Menu_95_named.tab",
                               sep = "\t", stringsAsFactors = FALSE,
                               check.names = FALSE)
# Resolve obitools TAXID → 7-rank taxonomy. pct_id arrives via the
# BEST_IDENTITY column on a 0-1 scale; reconcile_global_local auto-
# detects and rescales it to 0-100.
global_in <- resolve_taxids(obi_raw, taxid_col = "TAXID", sql_path = sql)
global_in$pct_id <- suppressWarnings(as.numeric(obi_raw$BEST_IDENTITY))
names(global_in)[names(global_in) == "ID"] <- "ASV_id"

# ----- Local DB side: vsearch --userout file -----
local_in <- parse_vsearch_userout(
  "script_03_inputs/userout_MiFish5-1db_Galapagos_top5_comprehensive_galapagos_results.txt"
)

cat("Inputs prepared:\n")
cat("  global_in:", nrow(global_in), "ASVs, pct_id range",
    paste(round(range(global_in$pct_id, na.rm = TRUE), 2), collapse = "-"), "\n")
cat("  local_in: ", nrow(local_in),  "ASVs, pct_id range",
    paste(round(range(local_in$pct_id,  na.rm = TRUE), 2), collapse = "-"), "\n\n")

# ----- REGATTA pipeline (3 lines) -----
rec  <- reconcile_global_local(global_in, local_in)             # writes reconcile_global_local_out/
post <- reconcile_checklist(rec$result, fish_checklist)         # writes reconcile_checklist_out/ (augmented)
summ <- summarize_regatta(rec, post, global_in, local_in)

cat("\nSummary table:\n")
print(summ)
