# Extract every land pixel's monthly climate series from the yearly
# TerraClimate rasters and compute agroclimatic indices, year by year.
#
# Pixels are processed in ~100 batches ("tiles") to keep memory usage and
# per-batch runtime manageable; each tile is written to its own CSV file.
#
# Originally written by Benjamin Bois, 2023-11-10.

pacman::p_load(dplyr, tidyverse, ggplot2, sf, terra, tidyterra, geodata)

source("R/config.R")
source("R/koppen.R")
source("R/agroclim_indices.R")
source("R/tiling.R")

setwd(data_root)
dir.create("TC_CSV_DATA", showWarnings = FALSE)

## ---- split all land pixels into ~100 tiles ---------------------------------

r <- rast(file.path("Download_TerraClimate", "tmin_2020.nc"))[[1]]
tiles <- make_pixel_tiles(r, n_tiles = 100)
n_tiles <- max(tiles$tile)
write_tile_list(tiles, "AAA_TilesLists_TerraClimate.txt")

## ---- process one tile at a time ---------------------------------------------

t_start <- Sys.time()
for (tile_id in 1:n_tiles) {
  tile_pix <- tiles[tiles$tile == tile_id, c("x", "y", "pixID")]

  # One extra year is needed: in the Southern Hemisphere, the growing season
  # spans two calendar years, so year N's indices need Jul-Dec of year N-1.
  dat <- c()
  for (year in 2000:2020) {
    year_dat <- tile_pix
    for (var in c("tmin", "tmax", "ppt", "pet")) {
      rast_TC <- rast(file.path("Download_TerraClimate", paste0(var, "_", year, ".nc")))
      year_dat <- cbind(year_dat, terra::extract(rast_TC, tile_pix[, 1:2]))
    }
    dat <- rbind(dat, cbind(year = year, year_dat))
  }

  # Some pixel/years are missing data (e.g. small islands not fully covered
  # by the raster grid) -- drop them.
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
  # (Jul..Jun), using the actual previous calendar year for Jul-Dec.
  if (length(isSH) > 0) {
    sh_dat <- dat[isSH, ]

    false_years <- sh_dat$year - 1
    false_years[false_years == min(false_years)] <- min(sh_dat$year)
    year_pix <- paste(sh_dat$pixID, sh_dat$year)
    false_year_pix <- paste(sh_dat$pixID, false_years)
    match_year_pix <- match(false_year_pix, year_pix)
    sh_dat[, c(TNInd[7:12], TXInd[7:12], RRInd[7:12], ETInd[7:12])] <-
      sh_dat[match_year_pix, c(TNInd[7:12], TXInd[7:12], RRInd[7:12], ETInd[7:12])]

    sh_dat[, TNInd] <- sh_dat[, TNInd[c(7:12, 1:6)]]
    sh_dat[, TXInd] <- sh_dat[, TXInd[c(7:12, 1:6)]]
    sh_dat[, RRInd] <- sh_dat[, RRInd[c(7:12, 1:6)]]
    sh_dat[, ETInd] <- sh_dat[, ETInd[c(7:12, 1:6)]]
    sh_dat[, TMInd] <- sh_dat[, TMInd[c(7:12, 1:6)]]

    sdat <- rbind(sdat, sh_dat)
  }

  # Previous winter (Oct-Dec of the previous calendar year, Jan-Mar of the
  # current one) precipitation/PET, needed for the winter water balance.
  false_years <- sdat$year - 1
  false_years[false_years == min(false_years)] <- min(sdat$year)
  year_pix <- paste(sdat$pixID, sdat$year)
  false_year_pix <- paste(sdat$pixID, false_years)
  match_year_pix <- match(false_year_pix, year_pix)
  winter_hydro <- sdat[match_year_pix, c(RRInd[10:12], ETInd[10:12])]
  names(winter_hydro) <- paste0("m1_", names(winter_hydro))
  sdat <- cbind(sdat, winter_hydro)
  M1RRInd <- grep("m1_ppt", colnames(sdat))
  M1ETInd <- grep("m1_pet", colnames(sdat))

  sdat <- as.matrix(sdat)
  isSH <- which(sdat[, "y"] < 0) # positions may have changed after rbind()

  agro_idx <- calc_agroclim_indices(
    sdat, TNInd, TXInd, RRInd, ETInd, TMInd, LATInd, isSH,
    winter_rr = sdat[, c(M1RRInd, RRInd[1:3])],
    winter_et = sdat[, c(M1ETInd, ETInd[1:3])]
  )

  agro <- cbind(sdat[, c("pixID", "x", "y", "year")], agro_idx)

  data.table::fwrite(agro, file.path("TC_CSV_DATA",
    paste0("TC_csv_data_TerraClimate_Obs_Tile_", tile_id, ".csv")))

  cat("Finished tile", tile_id, "of", n_tiles, "-", date(), "\n")
}
cat("Total time:", format(Sys.time() - t_start), "\n")
