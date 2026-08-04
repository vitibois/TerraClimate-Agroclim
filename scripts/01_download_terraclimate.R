# Download TerraClimate monthly NetCDF files (one file per variable/year).
#
# Originally written by Benjamin Bois (2025-05-14), based on code by
# Sebastien Nicolas (2023). Updated July 2025.
#
# Data source: https://www.climatologylab.org/terraclimate.html
# Catalog: http://thredds.northwestknowledge.net:8080/thredds/catalog/TERRACLIMATE_ALL/data/catalog.html

source("R/config.R")

setwd(file.path(data_root, "Download_TerraClimate"))

terraclimate_base_url <- "https://climate.northwestknowledge.net/TERRACLIMATE-DATA/"
terraclimate_prefix <- "TerraClimate"

year_min <- 2001
year_max <- 2020

# Variables to download
my_vars <- c("tmin", "tmax", "vpd", "pet", "ppt")

terraclimate_url <- function(year, var) {
  paste0(terraclimate_base_url, terraclimate_prefix, "_", var, "_", year, ".nc")
}

for (year in year_min:year_max) {
  t_start <- Sys.time()
  for (var in my_vars) {
    downloader::download(
      url = terraclimate_url(year = year, var = var),
      destfile = paste0(terraclimate_prefix, "_", var, "_", year, ".nc"),
      mode = "wb"
    )
  }
  cat("Done for year", year, "-", date(),
      "\nElapsed time for this loop:", format(Sys.time() - t_start), "\n")
}
