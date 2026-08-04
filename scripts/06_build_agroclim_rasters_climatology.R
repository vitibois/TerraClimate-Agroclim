# Compile the per-tile climatology CSVs (script 04) into one NetCDF raster
# per agroclimatic index. Unlike the annual pipeline, there is nothing to
# aggregate here: each pixel already has a single (climatological) value per
# index, simply copied onto the raster grid.
#
# Originally written by Benjamin Bois, 2023-11-10.
# Reworked 2025-07-31, B. Bois.

library(data.table)
library(terra)
library(ncdf4)

source("R/config.R")
setwd(data_root)

## ---- select which climatology to process ------------------------------------
# Must match the climatology processed by script 04.

# mync_dir <- "TerraClimate_2c_climatologies20y"
# mync_prefix <- "TerraClimate2C_"
# scenario <- "plus2"

# mync_dir <- "TerraClimate_4c_climatologies20y"
# mync_prefix <- "TerraClimate4C_"
# scenario <- "plus4"

# mync_dir <- "TerraClimate_obs8110_climatologies20y"
# mync_prefix <- "TerraClimate19812010_"
# scenario <- "CLIMOBS19812010"

mync_dir <- "TerraClimate_obs0120_climatologies20y"
mync_prefix <- "TerraClimate20012020_"
scenario <- "CLIMOBS20012020"

# mync_dir <- "TerraClimate_obs6190_climatologies30y"
# mync_prefix <- "TerraClimate19611990_"
# scenario <- "CLIMOBS19611990"

# Template raster (all cells set to NA), used to place values back on the grid.
rast_TC <- rast(file.path("Download_TerraClimate", "climatologies", mync_dir,
                           paste0(mync_prefix, "tmin.nc")))
r_na <- rast_TC[[1]]
values(r_na) <- NA

funname <- "Avg" # naming label only, for consistency with the annual outputs

mypath <- file.path(data_root, "AgroNetCDF")
dir.create(file.path(mypath, mync_dir), recursive = TRUE, showWarnings = FALSE)

tiles_df <- read.table(file.path("Download_TerraClimate", "climatologies", mync_dir,
                                  "AAA_TilesLists_TerraClimate.txt"), header = TRUE)

# Load all agroclimatic indices for all pixels
t_start <- Sys.time()
agg <- c()
for (tile_id in tiles_df$ID) {
  tile_agro <- fread(file.path("CLIMATOLOGIES_TC_CSV_DATA", mync_dir,
    paste0(mync_prefix, "TC_csv_data_Tile_", tile_id, ".csv")))
  agg <- rbind(agg, tile_agro)
}
cat("Loaded all tiles in", format(Sys.time() - t_start), "\n")

period <- agg$year[1]

my_vars <- c("HI", "CI", "DI", "Q", "WWB",
             "WI", "ANT", "ANP",
             "GST", "GSP", "GST49", "GSP49",
             "DST", "DSP", "RipeT", "RipeP",
             "FST", "FSP", "DORT",
             "SFR", "WCR", "HSI")

t_start <- Sys.time()
for (var in my_vars) {
  r <- r_na
  r[agg$pixID] <- as.data.frame(agg)[, var]

  ncname <- paste0(var, "_", funname, "_", scenario, "_", period)
  ncfname <- file.path(mypath, mync_dir, paste0(ncname, ".nc"))
  writeCDF(r, ncfname, var, longname = paste0(var, " ", funname, " - dataset: ",
    scenario, " period ", period), unit = "", split = FALSE, overwrite = TRUE)

  cat("Finished variable", var, "-", date(), "\n")
}
cat("Total time:", format(Sys.time() - t_start), "\n")
