# Average yearly TerraClimate NetCDF files into a multi-year monthly
# climatology (12 monthly rasters per variable, averaged over `year_min`:
# `year_max`).
#
# Originally written by Benjamin Bois (2025-05-14), based on code by
# Sebastien Nicolas (2023). Updated July 2025.

pacman::p_load(terra, sf, tidyverse, tidyterra)

source("R/config.R")

setwd(file.path(data_root, "Download_TerraClimate"))

year_min <- 2001
year_max <- 2020

# Output directory
outdir <- paste0("TerraClimate_obs", substr(year_min, 3, 4), substr(year_max, 3, 4),
                  "_climatologies20y")
ncpath <- file.path("climatologies", outdir)
dir.create(ncpath, recursive = TRUE, showWarnings = FALSE)
prefix <- paste0("TerraClimate", year_min, year_max)

# List of downloaded yearly files
nc_list <- list.files(pattern = "nc$")

my_vars <- c("tmin", "tmax", "vpd", "pet", "ppt")

for (var in my_vars) {
  my_files <- grep(paste0("_", var, "_"), nc_list, value = TRUE)
  r <- rast(my_files)
  crs(r) <- "epsg:4326"

  r_avg <- c()
  for (month in 1:12) {
    month_index <- grep(paste0("_", month, "$"), names(r))
    r_avg <- c(r_avg, mean(subset(r, month_index)))
  }
  r_avg <- rast(r_avg)

  ncfname <- file.path(ncpath, paste0(prefix, "_", var, ".nc"))
  writeCDF(r_avg, ncfname, varname = var, overwrite = TRUE)
}
