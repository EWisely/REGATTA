setup_explainer <- function() {
  print(c("See https://docs.ropensci.org/rgbif/articles/gbif_credentials.html, run usethis::edit_r_environ()
    in the console and input your GBIF login information: 
    GBIF_USER='jwaller'
    GBIF_PWD='safe_fake_password_123'
    GBIF_EMAIL='jwaller@gbif.org'",
    "Make an R project and folders called 'databases' and 'custom_db' in the working directory",
    "Put any existing local species databases you have into 'databases' and make sure 'Genus' and 'Species' columns are named exactly that.",
    "To define your local region, draw it on this webapp: https://wktmap.com, and copy the polygon into your location parameters",
    "devtools::install_github('james-thorson/FishLife') and devtools::install_github('cfree14/freeR')",
    "If you encounter errors with GBIF_Download, try: crul::set_opts(http_version = 2) (from https://github.com/ropensci/crul/issues/174#issuecomment-1599112561)",
    "AND options(timeout=1500) (from https://stackoverflow.com/questions/35282928/how-do-i-set-a-timeout-for-utilsdownload-file-in-r)"
    ))
}

# Test Run
# setup_explainer()