# Compile the per-tile annual-series CSVs (script 03) into one summary
# NetCDF raster per agroclimatic index: for each pixel, the value is
# averaged over all extracted years (2000 is discarded, see below).
#
# Originally written by Benjamin Bois, 2023-11-10.

library(data.table)
library(terra)
library(ncdf4)

source("R/config.R")
setwd(data_root)

tiles_df <- read.table("AAA_TilesLists_TerraClimate.txt", header = TRUE)

# Template raster (all cells set to NA), used to place the aggregated values
# back on the grid.
rast_TC <- rast(file.path("Download_TerraClimate", "tmin_2020.nc"))
r_na <- rast_TC[[1]]
values(r_na) <- NA

myfun <- mean
funname <- "Avg"

scenario <- "OBS"
period <- "2001-2020"
# 2000 is only extracted to provide the Southern Hemisphere's Jul-Dec data
# for the 2001 growing season (see script 03); it is not itself a full,
# usable year and must be discarded here.
if (scenario == "OBS") year2discard <- 2000
mypath <- file.path(data_root, "AgroNetCDF")
dir.create(mypath, showWarnings = FALSE)

# TODO: KOPPEN is deliberately excluded here: it is a categorical code, so
# averaging it across years (as done below for the other indices) is not
# meaningful -- it would need a per-pixel mode instead.
my_vars <- c("HI", "CI", "DI", "Q", "WWB", "WI", "ANT", "ANP",
             "GST", "GSP", "GST49", "GSP49", "DSP", "DST",
             "RipeT", "RipeP", "FSP", "FST", "DORT", "SFR", "WCR", "HSI")

t_start <- Sys.time()
# Roughly 5 to 8 minutes per index, depending on the computation involved.
for (var in my_vars) {
  agg <- c()
  for (tile_id in tiles_df$ID) {
    df <- fread(file.path("TC_CSV_DATA",
      paste0("TC_csv_data_TerraClimate_Obs_Tile_", tile_id, ".csv")))
    df <- df[year != year2discard]
    tile_agg <- df[, .(value = myfun(get(var))), by = pixID]
    agg <- rbind(agg, tile_agg)
  }
  r <- r_na
  r[agg$pixID] <- agg$value

  ncname <- paste0(var, "_", funname, "_", scenario, "_", period)
  ncfname <- file.path(mypath, paste0(ncname, ".nc"))
  writeCDF(r, ncfname, var, longname = paste0(var, " ", funname, " - dataset: ",
    scenario, " period ", period), unit = "", split = FALSE, overwrite = TRUE)

  cat("Finished variable", var, "-", date(), "\n")
}
cat("Total time:", format(Sys.time() - t_start), "\n")
