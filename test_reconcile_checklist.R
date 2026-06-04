# test_regatta_checklist_lca.R
# Synthetic 6-row test for regatta_checklist_lca(). Run from project root.

source("reconcile_checklist.R")

# Tiny regional checklist: 3 species across 2 families and 2 orders.
checklist <- data.frame(
  domain  = "Eukaryota",
  phylum  = "Chordata",
  class   = "Actinopterygii",
  order   = c("Scorpaeniformes", "Scorpaeniformes", "Clupeiformes"),
  family  = c("Sebastidae",      "Sebastidae",      "Engraulidae"),
  genus   = c("Sebastes",        "Sebastes",        "Engraulis"),
  species = c("Sebastes mystinus","Sebastes paucispinis","Engraulis mordax"),
  stringsAsFactors = FALSE
)

# 6 test ASVs covering the main cases:
# ASV_1: species match — unchanged
# ASV_2: species NOT on checklist, genus IS — downgrade to genus
# ASV_3: species/genus/family miss, order matches — downgrade to order
# ASV_4: species match — unchanged
# ASV_5: Bacteria, nothing matches — all NA'd
# ASV_6: species/genus/family NA in input, but order matches — match at order
input <- data.frame(
  ASV_id  = paste0("ASV_", 1:6),
  domain  = c("Eukaryota",       "Eukaryota",      "Eukaryota",        "Eukaryota",       "Bacteria",        "Eukaryota"),
  phylum  = c("Chordata",        "Chordata",       "Chordata",         "Chordata",        "Firmicutes",      "Chordata"),
  class   = c("Actinopterygii",  "Actinopterygii", "Actinopterygii",   "Actinopterygii",  "Bacilli",         "Actinopterygii"),
  order   = c("Scorpaeniformes", "Scorpaeniformes","Clupeiformes",     "Clupeiformes",    "Bacillales",      "Clupeiformes"),
  family  = c("Sebastidae",      "Sebastidae",     "Clupeidae",        "Engraulidae",     "Bacillaceae",     NA),
  genus   = c("Sebastes",        "Sebastes",       "Sardinops",        "Engraulis",       "Bacillus",        NA),
  species = c("Sebastes mystinus","Sebastes goodei","Sardinops sagax", "Engraulis mordax","Bacillus subtilis",NA),
  pct_id  = c(99.5,              97.2,              95.0,              100.0,              98.0,             88.0),
  stringsAsFactors = FALSE
)

result <- reconcile_checklist(input, checklist, output_dir = NULL)

# --- $result is exactly 8 columns: ASV_id + 7 ranks. Nothing else. ---
stopifnot(identical(names(result$result),
                    c("ASV_id","domain","phylum","class","order","family","genus","species")))
re <- result$result

# ASV_1: full species match, unchanged
stopifnot(re$species[1] == "Sebastes mystinus")
stopifnot(re$genus[1]   == "Sebastes")

# ASV_2: downgraded to genus
stopifnot(is.na(re$species[2]))
stopifnot(re$genus[2]  == "Sebastes")
stopifnot(re$family[2] == "Sebastidae")

# ASV_3: downgraded to order
stopifnot(is.na(re$species[3]))
stopifnot(is.na(re$genus[3]))
stopifnot(is.na(re$family[3]))
stopifnot(re$order[3] == "Clupeiformes")
stopifnot(re$class[3] == "Actinopterygii")

# ASV_4: full species match
stopifnot(re$species[4] == "Engraulis mordax")

# ASV_5: no match — every rank NA'd
stopifnot(is.na(re$domain[5]))
stopifnot(is.na(re$phylum[5]))
stopifnot(is.na(re$species[5]))

# ASV_6: order match with NAs already in input below order
stopifnot(re$domain[6] == "Eukaryota")
stopifnot(re$phylum[6] == "Chordata")
stopifnot(re$class[6]  == "Actinopterygii")
stopifnot(re$order[6]  == "Clupeiformes")
stopifnot(is.na(re$family[6]))
stopifnot(is.na(re$genus[6]))
stopifnot(is.na(re$species[6]))

# --- $tracking carries the bookkeeping ---
tr <- result$tracking
stopifnot("regatta_match_rank" %in% names(tr))
stopifnot("pct_id" %in% names(tr))   # input metadata preserved
stopifnot("before_species" %in% names(tr))
stopifnot("after_species"  %in% names(tr))
stopifnot("after_scientific_name" %in% names(tr))

# Match-rank assertions go on $tracking
stopifnot(tr$regatta_match_rank[1] == "species")
stopifnot(tr$regatta_match_rank[2] == "genus")
stopifnot(tr$regatta_match_rank[3] == "order")
stopifnot(tr$regatta_match_rank[4] == "species")
stopifnot(is.na(tr$regatta_match_rank[5]))
stopifnot(tr$regatta_match_rank[6] == "order")

# Before/after columns
stopifnot(tr$before_species[2] == "Sebastes goodei")
stopifnot(is.na(tr$after_species[2]))
stopifnot(tr$before_species[5] == "Bacillus subtilis")
stopifnot(is.na(tr$after_species[5]))

# Input metadata preserved
stopifnot(all(tr$pct_id == input$pct_id))
stopifnot(all(tr$ASV_id == input$ASV_id))

# --- $stats has counts ---
st <- result$stats
stopifnot(is.data.frame(st))
stopifnot(all(c("metric", "count") %in% names(st)))
stopifnot(st$count[st$metric == "total ASVs"] == 6)
stopifnot(st$count[st$metric == "matched at species"] == 2)  # ASV_1, ASV_4
stopifnot(st$count[st$metric == "matched at order"]   == 2)  # ASV_3, ASV_6
stopifnot(st$count[st$metric == "not matched (off-target)"] == 1)  # ASV_5

cat("All reconcile_checklist tests passed.\n")
