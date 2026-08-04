# Extract every land pixel's monthly climate values from a TerraClimate
# climatology (a single averaged "year") and compute agroclimatic indices.
#
# Pixels are processed in ~100 batches ("tiles"); each tile is written to
# its own CSV file.
#
# Originally written by Benjamin Bois, 2023-11-10.
# Reworked 2025-07-31, B. Bois.

pacman::p_load(tidyverse, sf, terra, tidyterra, geodata)

source("R/config.R")
source("R/koppen.R")
source("R/agroclim_indices.R")
source("R/tiling.R")

setwd(file.path(data_root, "CLIMATOLOGIES_TC_CSV_DATA"))

## ---- select which climatology to process ------------------------------------
# "Fake" year for the +2C / +4C scenarios: the year at which the RCP8.5
# ensemble mean of 39 CMIP5 GCMs reaches +2C / +4C above pre-industrial
# levels (Wang, Jiang & Lang, 2018, https://doi.org/10.1007/s00376-018-7160-4).
# Uncomment the block corresponding to the climatology to process.

# mync_dir <- "TerraClimate_2c_climatologies20y"
# mync_prefix <- "TerraClimate2C_"
# year <- 2042

# mync_dir <- "TerraClimate_4c_climatologies20y"
# mync_prefix <- "TerraClimate4C_"
# year <- 2084

# mync_dir <- "TerraClimate_obs8110_climatologies20y"
# mync_prefix <- "TerraClimate19812010_"
# year <- 1995

mync_dir <- "TerraClimate_obs0120_climatologies20y"
mync_prefix <- "TerraClimate20012020_"
year <- 2010

# mync_dir <- "TerraClimate_obs6190_climatologies30y"
# mync_prefix <- "TerraClimate19611990_"
# year <- 1976

climatology_nc <- function(var) {
  file.path(data_root, "Download_TerraClimate", "climatologies", mync_dir,
            paste0(mync_prefix, var, ".nc"))
}

## ---- split all land pixels into ~100 tiles ---------------------------------

r <- rast(climatology_nc("tmin"))[[1]]
tiles <- make_pixel_tiles(r, n_tiles = 100)
n_tiles <- max(tiles$tile)
write_tile_list(tiles, file.path(data_root, "Download_TerraClimate", "climatologies",
                                  mync_dir, "AAA_TilesLists_TerraClimate.txt"))

## ---- process one tile at a time ---------------------------------------------

t_start <- Sys.time()
for (tile_id in 1:n_tiles) {
  tile_pix <- tiles[tiles$tile == tile_id, c("x", "y", "pixID")]

  dat <- tile_pix
  for (var in c("tmin", "tmax", "ppt", "pet")) {
    rast_TC <- rast(climatology_nc(var))
    dat <- cbind(dat, terra::extract(rast_TC, tile_pix[, 1:2]))
    names(dat) <- sub("Z1=", "", names(dat))
  }
  dat <- cbind(year = year, dat)

  # Some pixels are missing data (e.g. small islands not fully covered by
  # the raster grid) -- drop them.
  col_sums <- Rfast::colsums(as.matrix(dat))
  na_cols <- which(is.na(col_sums))
  if (length(na_cols) > 0) {
    rows_to_remove <- c()
    for (col in na_cols) rows_to_remove <- c(rows_to_remove, which(is.na(dat[, col])))
    dat <- dat[-rows_to_remove, ]
  }

  # Mean monthly temperature
  TNInd <- grep("tmin", names(dat))
  TXInd <- grep("tmax", names(dat))
  tavg <- (dat[, TNInd] + dat[, TXInd]) / 2
  names(tavg) <- sub("tmin", "tavg", names(tavg))
  dat <- cbind(dat, tavg)
  rm(tavg)

  isSH <- which(dat$y < 0)

  # Northern Hemisphere data (calendar year = agricultural year, Jan..Dec)
  sdat <- dat[-isSH, ]
  if (length(isSH) == 0) sdat <- dat
  TNInd <- grep("tmin", names(sdat))
  TXInd <- grep("tmax", names(sdat))
  RRInd <- grep("ppt", names(sdat))
  ETInd <- grep("pet", names(sdat))
  LATInd <- grep("^y$", names(sdat))
  TMInd <- grep("tavg", names(sdat))

  # Southern Hemisphere data: reorder into the local agricultural year
  # (Jul..Jun). A climatology has no "previous year": Jul-Dec is simply
  # taken from the same averaged year.
  if (length(isSH) > 0) {
    sh_dat <- dat[isSH, ]
    sh_dat[, TNInd] <- sh_dat[, TNInd[c(7:12, 1:6)]]
    sh_dat[, TXInd] <- sh_dat[, TXInd[c(7:12, 1:6)]]
    sh_dat[, RRInd] <- sh_dat[, RRInd[c(7:12, 1:6)]]
    sh_dat[, ETInd] <- sh_dat[, ETInd[c(7:12, 1:6)]]
    sh_dat[, TMInd] <- sh_dat[, TMInd[c(7:12, 1:6)]]

    sdat <- rbind(sdat, sh_dat)
  }

  sdat <- as.matrix(sdat)
  isSH <- which(sdat[, "y"] < 0) # positions may have changed after rbind()

  agro_idx <- calc_agroclim_indices(
    sdat, TNInd, TXInd, RRInd, ETInd, TMInd, LATInd, isSH,
    winter_rr = sdat[, c(RRInd[10:12], RRInd[1:3])],
    winter_et = sdat[, c(ETInd[10:12], ETInd[1:3])]
  )

  agro <- cbind(sdat[, c("pixID", "x", "y", "year")], agro_idx)

  dir.create(mync_dir, showWarnings = FALSE)
  data.table::fwrite(agro, file.path(mync_dir,
    paste0(mync_prefix, "TC_csv_data_Tile_", tile_id, ".csv")))

  cat("Finished tile", tile_id, "of", n_tiles, "-", date(), "\n")
}
cat("Total time:", format(Sys.time() - t_start), "\n")
