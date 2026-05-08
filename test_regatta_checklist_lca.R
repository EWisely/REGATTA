# test_regatta_checklist_lca.R
# Synthetic 6-row test for regatta_checklist_lca(). Run from project root.

source("regatta_checklist_lca.R")

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

result <- regatta_checklist_lca(input, checklist)

# --- corrected ---
co <- result$corrected

# ASV_1: full species match, unchanged
stopifnot(co$regatta_match_rank[1] == "species")
stopifnot(co$species[1] == "Sebastes mystinus")
stopifnot(co$genus[1]   == "Sebastes")

# ASV_2: downgraded to genus
stopifnot(co$regatta_match_rank[2] == "genus")
stopifnot(is.na(co$species[2]))
stopifnot(co$genus[2]  == "Sebastes")
stopifnot(co$family[2] == "Sebastidae")

# ASV_3: downgraded to order
stopifnot(co$regatta_match_rank[3] == "order")
stopifnot(is.na(co$species[3]))
stopifnot(is.na(co$genus[3]))
stopifnot(is.na(co$family[3]))
stopifnot(co$order[3] == "Clupeiformes")
stopifnot(co$class[3] == "Actinopterygii")

# ASV_4: full species match
stopifnot(co$regatta_match_rank[4] == "species")
stopifnot(co$species[4] == "Engraulis mordax")

# ASV_5: no match — every rank NA'd, match_rank NA
stopifnot(is.na(co$regatta_match_rank[5]))
stopifnot(is.na(co$domain[5]))
stopifnot(is.na(co$phylum[5]))
stopifnot(is.na(co$species[5]))

# ASV_6: order match with NAs already in input below order
stopifnot(co$regatta_match_rank[6] == "order")
stopifnot(co$domain[6] == "Eukaryota")
stopifnot(co$phylum[6] == "Chordata")
stopifnot(co$class[6]  == "Actinopterygii")
stopifnot(co$order[6]  == "Clupeiformes")
stopifnot(is.na(co$family[6]))
stopifnot(is.na(co$genus[6]))
stopifnot(is.na(co$species[6]))

# Non-rank columns preserved
stopifnot(all(co$pct_id == input$pct_id))
stopifnot(all(co$ASV_id == input$ASV_id))

# --- changes ---
ch <- result$changes

# Should include exactly ASV_2, ASV_3, ASV_5; NOT ASV_1, ASV_4, ASV_6 (ASV_6 had NAs
# below order before AND after, so nothing actually changed for it).
stopifnot(setequal(ch$ASV_id, c("ASV_2", "ASV_3", "ASV_5")))

asv2 <- ch[ch$ASV_id == "ASV_2", ]
stopifnot(asv2$before_species == "Sebastes goodei")
stopifnot(is.na(asv2$after_species))
stopifnot(asv2$before_genus == "Sebastes")
stopifnot(asv2$after_genus  == "Sebastes")
stopifnot(asv2$match_rank   == "genus")

asv5 <- ch[ch$ASV_id == "ASV_5", ]
stopifnot(asv5$before_species == "Bacillus subtilis")
stopifnot(is.na(asv5$after_species))
stopifnot(is.na(asv5$after_domain))
stopifnot(is.na(asv5$match_rank))

# --- before ---
# Should be a faithful snapshot of the input
stopifnot(identical(result$before, input))

cat("All regatta_checklist_lca tests passed.\n")
