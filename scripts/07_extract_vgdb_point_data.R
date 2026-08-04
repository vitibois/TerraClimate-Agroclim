# Extract agroclimatic index rasters (all periods/scenarios found under
# AgroNetCDF/) and terrain features (elevation, slope) at the points of the
# Winegrape Vineyard Geodatabase (VGDB).
#
# Originally written by Benjamin Bois, 2023.
# Reworked 2025-07-31, B. Bois: generalised to extract all periods directly.

library(terra)
library(sf)
library(elevatr)
library(ncdf4)

source("R/config.R")

setwd(file.path(data_root, "AgroNetCDF"))

# Winegrape Vineyard Geodatabase
vgdbpath <- file.path(sig_root, "world/Wine_Vineyards_GeoDB/v2/v2.3",
                       "Points_Winegrape_Vineyard_Geodatabase_V2.3.gpkg")
vgdb.version <- sub("\\.gpkg", "", sub(".*Geodatabase_v", "", vgdbpath, ignore.case = TRUE))
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

dir_out <- file.path(data_root, "Extraction_TerraClimatPoints", paste0("VGDB_v", vgdb.version))
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
name_out <- paste0("AgroIndices_TerraClimate_VGDB_Pts_v", vgdb.version, "AllPeriods.csv")
data.table::fwrite(bigout, file.path(dir_out, name_out), row.names = FALSE)
