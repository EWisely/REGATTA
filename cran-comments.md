## Submission

New submission. REGATTA creates a regional species checklist from publicly available biodiversity occurrence databases, flags non-local species assignments in eDNA
metabarcoding results, and corrects them by reconciling classifier output taxonomic levels against a regional
species checklist.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
    Maintainer: 'Eldridge Wisely <eldridge.wisely@gmail.com>'
    New submission

(A local "unable to verify current time" note also appears when checking
offline; it is environmental and not expected on CRAN's machines.)

## Test environments

* local: macOS, R 4.4.0
* GitHub Actions (r-lib/actions):
  - macos-latest   (R release)
  - windows-latest (R release)
  - ubuntu-latest  (R devel, release, oldrel-1)

## Notes for the reviewer

* Some examples use \dontrun{} because they require resources unavailable
  during checks: GBIF account credentials, a multi-gigabyte local NCBI taxonomy
  database (taxonomizr), and live WoRMS/GBIF/OBIS web services. Runnable,
  self-contained demonstrations of the core functions are provided in the
  vignettes using small bundled data fixtures.

* No function writes to the user's filespace by default; all output-producing
  functions return their results and write files only when given an output
  directory.

## Downstream dependencies

None (new package).
