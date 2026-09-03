# Extract agroclimatic index rasters (all periods/scenarios found under
# AgroNetCDF/), terrain features (elevation, slope), and raw monthly climate
# variables (one period at a time) at the points of the Winegrape Vineyard
# Geodatabase (VGDB).
#
# Originally written by Benjamin Bois, 2023.
# Reworked 2025-07-31, B. Bois: generalised to extract all periods directly.
# Reworked 2026-09-03, B. Bois: added the monthly climate variables section.

library(terra)
library(sf)
library(elevatr)
library(ncdf4)

source("R/config.R")

setwd(file.path(data_root, "AgroNetCDF"))

# Winegrape Vineyard Geodatabase (version set in R/config.R; scripts 08/09
# read this script's output and must agree on the same version)
vgdbpath <- file.path(sig_root, "world/Wine_Vineyards_GeoDB/v2", paste0("v", vgdb_version),
                       paste0("Points_Winegrape_Vineyard_Geodatabase_V", vgdb_version, ".gpkg"))
mypts <- vect(vgdbpath)
mypts$CN_REG <- paste(mypts$CN, mypts$WINE_REGION, sep = "_")

x <- geom(mypts)[, "x"]
y <- geom(mypts)[, "y"]

## ---- terrain features (elevation, slope) -------------------------------------
# get_elev_raster() with a terra SpatVector errors with "argument inutilise
# (values = FALSE)": elevatr's internal loc_check() calls
# terra::as.points(locations, values = FALSE), and that argument no longer
# exists on this terra version for SpatVector input. Converting to sf first
# routes through a working code path.
#
# Source: NASA Shuttle Radar Topography Mission (SRTM) (2013). Distributed by
# OpenTopography. https://doi.org/10.5069/G9445JDF
# GL1: SRTM GL1 (1 arc-second, ~30 m at the equator).
# If you hit OpenTopography rate limits, set your own API key first with
# elevatr::set_opentopo_key("<your key>") (see the elevatr documentation).

get_terrain_srtm <- function(pt) {
  pt_sf <- pt |> terra::buffer(width = 90) |> sf::st_as_sf()
  r_elev <- get_elev_raster(pt_sf, prj = "epsg:4326", src = "gl1", clip = "location")
  r_slope <- terrain(r_elev, v = "slope", unit = "radians") |> tan()
  elev <- rast(r_elev) |> terra::extract(pt, raw = TRUE, ID = FALSE) |> as.vector()
  slope <- rast(r_slope) |> terra::extract(pt, raw = TRUE, ID = FALSE) |> as.vector()
  data.frame(Elevation = elev, SlopePC = slope * 100)
}

terrain_df <- c()
t_start <- Sys.time()
for (i in 1:length(mypts)) {
  terrain_df <- rbind(terrain_df, cbind(data.frame(ID = i), get_terrain_srtm(mypts[i, ])))
  cat("Done for point", i, "/", length(mypts), "- elapsed:",
      format(Sys.time() - t_start), "\n")
}

## ---- agroclimatic indices, all periods/scenarios -----------------------------

mync_dirs <- sub("^\\./", "", grep("TerraClimate_", list.dirs(), value = TRUE))

# Reference variable names as written by scripts 05/06 (avg over years or
# climatology value); used to detect each scenario directory's file-naming
# suffix ("_<scenario>_<period>.nc").
my_vars <- c("CI_Avg", "DI_Avg", "DORT_Avg", "DSP_Avg", "DST_Avg",
             "FSP_Avg", "FST_Avg", "GSP_Avg", "GSP49_Avg", "GST_Avg", "GST49_Avg",
             "HI_Avg", "HSI_Avg", "KOPPEN_Avg", "Q_Avg", "RipeP_Avg", "SFR_Avg",
             "WCR_Avg", "WWB_Avg", "WI_Avg", "ANT_Avg", "ANP_Avg")

bigout <- c() # combined table of all indices, for all scenarios
for (mync_dir in mync_dirs) {
  myncs <- list.files(pattern = ".nc$", path = mync_dir)
  ref_file <- grep(paste0("^", my_vars[1]), myncs, value = TRUE)
  suffix <- sub(my_vars[1], "", ref_file)
  myncs <- paste0(my_vars, suffix)
  nc_info <- do.call(rbind, strsplit(myncs, "_"))
  myvar <- nc_info[, 1]
  mystat <- nc_info[, 2]
  myscen <- nc_info[, 3]

  out <- NULL
  for (i in seq_along(myncs)) {
    mync <- rast(file.path(mync_dir, myncs[i]))
    if (is.null(out)) {
      cells <- terra::extract(mync, mypts, cells = TRUE)
      out <- data.frame(pixID = cells$cell, x = x, y = y,
                         CN_REG = mypts$CN_REG, scenario = myscen[i])
    }
    out <- cbind(out, terra::extract(mync, mypts)[, 2])
    colnames(out)[ncol(out)] <- paste0(myvar[i], "_", mystat[i])
  }
  # A few VGDB points can fall in the same raster cell; keep one row per
  # unique location.
  out <- out[!duplicated(paste(out$x, out$y)), ]

  bigout <- rbind(bigout, out)
}

dir_out <- file.path(data_root, "Extraction_TerraClimatPoints", paste0("VGDB_v", vgdb_version))
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
name_out <- paste0("AgroIndices_TerraClimate_VGDB_Pts_v", vgdb_version, "AllPeriods.csv")
data.table::fwrite(bigout, file.path(dir_out, name_out), row.names = FALSE)

## ---- monthly climate variables, one climatology period ------------------------
# Extracts raw monthly tmin/tmax/ppt/pet (4 variables x 12 months = 48
# columns) at VGDB points, for a single climatology period at a time -- pick
# mync_dir/mync_prefix below (matching script 04's selection pattern; the
# prefix is the one written by script 02_build_climatologies.R).

# mync_dir <- "TerraClimate_2c_climatologies20y"
# mync_prefix <- "TerraClimate2C_"

# mync_dir <- "TerraClimate_4c_climatologies20y"
# mync_prefix <- "TerraClimate4C_"

# mync_dir <- "TerraClimate_obs8110_climatologies20y"
# mync_prefix <- "TerraClimate19812010_"

mync_dir <- "TerraClimate_obs0120_climatologies20y"
mync_prefix <- "TerraClimate20012020_"

# mync_dir <- "TerraClimate_obs6190_climatologies30y"
# mync_prefix <- "TerraClimate19611990_"

climatology_nc <- function(var) {
  file.path(data_root, "Download_TerraClimate", "climatologies", mync_dir,
            paste0(mync_prefix, var, ".nc"))
}

monthly_out <- data.frame(x = x, y = y, CN_REG = mypts$CN_REG)
for (var in c("tmin", "tmax", "ppt", "pet")) {
  mync <- rast(climatology_nc(var))
  if (!"pixID" %in% names(monthly_out)) {
    cells <- terra::extract(mync, mypts, cells = TRUE)
    monthly_out <- cbind(pixID = cells$cell, monthly_out)
  }
  monthly_out <- cbind(monthly_out, terra::extract(mync, mypts)[, -1])
}
names(monthly_out) <- sub("Z1=", "", names(monthly_out))

# A few VGDB points can fall in the same raster cell; keep one row per
# unique location.
monthly_out <- monthly_out[!duplicated(paste(monthly_out$x, monthly_out$y)), ]

name_out_monthly <- paste0("MonthlyClimate_TerraClimate_VGDB_Pts_v", vgdb_version, "_", mync_dir, ".csv")
data.table::fwrite(monthly_out, file.path(dir_out, name_out_monthly), row.names = FALSE)
