# REGATTA

**Reproducible Eco-Geographic Assignment Through Taxonomic Adjustment**

An R package (in development) for correcting off-target species assignments in
eDNA metabarcoding results using a regional species checklist, without manual
curation and while preserving as much taxonomic specificity as possible.

## The problem

Global reference databases (NCBI, EMBL, etc.) confidently assign eDNA reads
to species that don't actually live in your study area, because the closest
match in the database is a relative from somewhere else in the world. The
usual fixes — hand-curating a local reference database, or applying a flat
percent-identity cutoff — are either non-reproducible or sacrifice
specificity unnecessarily.

## The fix

REGATTA pulls a regional species list from public biodiversity sources
(GBIF, OBIS, plus any local CSVs you have), resolves it to NCBI taxonomy,
and uses it as a sanity filter on each ASV's classification. For each
assigned ASV, REGATTA finds the **lowest taxonomic rank shared between
the assignment and the regional checklist** and downgrades only that far.

A *Sebastes mystinus* call in a region that has *S. mystinus* on the
checklist stays at species. A call to *S. goodei* (a Pacific rockfish) in
a region whose checklist has other *Sebastes* but not *S. goodei* is
downgraded to genus *Sebastes* and flagged. A call to a tropical species
in a temperate-region checklist with no *Sebastes* at all gets walked up
the tree until something matches, or flagged as fully off-target if
nothing does.

## Pipeline

```mermaid
flowchart TD
    A[GBIF_download] --> D[build_regional_checklist]
    B[OBIS_download] --> D
    C[Local_csv_download] --> D
    D --> E[taxonomize_checklist]
    E --> F[("Regional checklist<br/>7-rank table<br/>one taxonomic group")]

    G["eDNA classifier output<br/>obitools / vsearch+SINTAX /<br/>Kraken2 / BLAST / etc."] --> H{Input shape?}
    H -->|SINTAX strings| I[parse_sintax]
    H -->|NCBI taxIDs| J[resolve_taxids]
    H -->|7 rank columns| K[no preprocessing needed]
    I --> L[("Taxonomy table<br/>ASV + 7 ranks")]
    J --> L
    K --> L

    F --> M[regatta_checklist_lca]
    L --> M
    M --> N["corrected<br/>(input table with ranks rewritten)"]
