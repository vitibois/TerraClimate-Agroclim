# Agroclimatic indices for viticulture, computed from monthly climate data.
#
# The core calculation (calc_agroclim_indices()) implements the Tonietto &
# Carbonneau (2004) multicriteria geoviticultural classification (Huglin
# Heliothermal Index, Cool Night Index, Riou Dryness Index), plus a set of
# complementary seasonal indices used throughout this project (Winkler Index,
# growing-season temperature/precipitation, frost/heat extremes...).
#
# It is shared by the two extraction pipelines (annual series and
# climatologies, see scripts/03_... and scripts/04_...), which differ only in
# how they assemble the monthly input matrices (e.g. how the "previous
# winter" precipitation/PET are obtained) but not in the index formulas
# themselves.
#
# Reference:
# Tonietto, J., & Carbonneau, A. (2004). A multicriteria climatic
# classification system for grape-growing regions worldwide. Agricultural
# and Forest Meteorology, 124(1-2), 81-97.
# https://doi.org/10.1016/j.agrformet.2003.06.001

#' Days in each calendar month (non-leap year), January to December.
DAYS_IN_MONTH <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

#' Build an (n x 12) matrix of days-in-month, reordered so that column 1 is
#' the first month of the local "agricultural year" for each row.
#'
#' Northern hemisphere rows keep the calendar order (Jan..Dec, columns 1:12).
#' Southern hemisphere rows are shifted by 6 months (Jul..Jun) to match the
#' reordering already applied to the tmin/tmax/ppt/pet columns upstream, so
#' that column 4 is always "the 4th month of the growing year" (April for the
#' Northern Hemisphere, October for the Southern Hemisphere), etc.
#'
#' @param n number of rows (pixels/pixel-years).
#' @param isSH integer vector, row indices that are in the Southern
#'   Hemisphere. Defaults to none (all Northern Hemisphere).
#'
#' @return an (n x 12) numeric matrix.
days_in_month_matrix <- function(n, isSH = integer(0)) {
  m <- matrix(DAYS_IN_MONTH, nrow = n, ncol = 12, byrow = TRUE)
  if (length(isSH) > 0) m[isSH, ] <- m[isSH, c(7:12, 1:6)]
  m
}

# Huglin's "k" coefficient (day-length correction), Hall & Jones (2010)
# continuous estimate, as used by Tonietto & Carbonneau (2004).
#
# Vendored from Benjamin Bois' personal function library (E:/R_functions/gccm.r).
#
# Reference:
# Hall, A., & Jones, G. V. (2010). Spatial analysis of climate in
# winegrape-growing regions in Australia. Australian Journal of Grape and
# Wine Research, 16(3), 389-404.
#' @param lat latitude in decimal degrees (sign gives the hemisphere).
k.hi.hall <- function(lat, mink = 1, maxk = 1.2, minlat = 34, maxlat = 65) {
  lat <- abs(lat)
  khi <- 5.552e-01 + 2.856e-02 * lat - 6.356e-04 * lat^2 + 5.316e-06 * lat^3
  khi[lat < minlat] <- pmax(mink, khi[lat < minlat])
  khi[lat > maxlat] <- pmin(maxk, khi[lat > maxlat])
  return(khi)
}

#' Compute agroclimatic indices from monthly climate data.
#'
#' All monthly matrices must already be reordered into the local
#' "agricultural year" (column 1 = first month of the growing year: January
#' for the Northern Hemisphere, July for the Southern Hemisphere - see
#' [days_in_month_matrix()]). This reordering, and the handling of the
#' "previous winter" precipitation/PET (`winter_rr`/`winter_et`), are done by
#' the calling script because they differ between pipelines: the annual
#' series pipeline looks up the actual previous calendar year, while the
#' climatology pipeline reuses the same (average) year's Oct-Mar months.
#'
#' @param sdat numeric matrix, one row per pixel (or pixel/year).
#' @param TNInd,TXInd,RRInd,ETInd,TMInd length-12 integer vectors: column
#'   indices in `sdat` of monthly tmin, tmax, precipitation, PET and mean
#'   temperature (in local agricultural-year order, Jan..Dec position).
#' @param LATInd integer, column index of latitude (`y`) in `sdat`.
#' @param isSH integer vector, row indices of `sdat` in the Southern
#'   Hemisphere.
#' @param winter_rr,winter_et (n x 6) matrices: precipitation and PET for the
#'   6 months preceding the growing season (Oct, Nov, Dec, Jan, Feb, Mar, in
#'   local agricultural-year order).
#'
#' @return a data.frame with one row per input row and the following
#'   columns:
#'   - `HI`: Huglin Heliothermal Index (Apr-Sep)
#'   - `CI`: Cool Night Index (mean min. temperature of the ripening month, Sep)
#'   - `DI`: Riou Dryness Index (soil water balance, Apr-Sep)
#'   - `Q`: winter runoff / surplus (Oct-Mar, water balance capped at 200 mm)
#'   - `WWB`: winter water balance (Oct-Mar)
#'   - `WI`: Winkler Index (growing degree-days, base 10C, Apr-Oct)
#'   - `ANT`, `ANP`: mean annual temperature and total annual precipitation
#'   - `KOPPEN`: Koppen-Geiger climate classification (Peel code, see [koppen.fast()])
#'   - `GST`, `GSP`: mean temperature / total precipitation, Apr-Oct
#'   - `GST49`, `GSP49`: mean temperature / total precipitation, Apr-Sep
#'   - `DST`, `DSP`: mean temperature / total precipitation, Apr-Jul (spring)
#'   - `RipeT`, `RipeP`: mean temperature / total precipitation of the
#'     ripening month (Sep)
#'   - `FST`, `FSP`: mean temperature / total precipitation, Jul-Sep (fruit development)
#'   - `DORT`: mean temperature of the dormant season (Oct-Mar)
#'   - `SFR`: minimum temperature, Apr-Jun (spring frost risk)
#'   - `WCR`: minimum temperature over the whole year (winter cold risk)
#'   - `HSI`: maximum temperature over the whole year (heat stress index)
calc_agroclim_indices <- function(sdat, TNInd, TXInd, RRInd, ETInd, TMInd, LATInd,
                                   isSH, winter_rr, winter_et) {
  ndaysmat <- days_in_month_matrix(nrow(sdat), isSH)

  # WI - Winkler Index (growing degree-days, base 10C, Apr-Oct)
  WI <- sdat[, TMInd[4:10]] - 10
  WI[WI < 0] <- 0
  WI <- WI * ndaysmat[, 4:10]
  WI <- rowSums(WI)

  # HI - Huglin Heliothermal Index (base 10C, Apr-Sep, latitude-corrected)
  khi <- k.hi.hall(lat = sdat[, LATInd])
  HI <- sdat[, c(TMInd[4:9], TXInd[4:9])] - 10
  HI[HI < 0] <- 0
  HI <- (HI[, 1:6] + HI[, 7:12]) / 2
  HI <- HI * ndaysmat[, 4:9]
  HI <- rowSums(HI) * khi

  # TODO: Gladstones Biologically Effective Degree Days (BEDD), as revised in
  # Gladstones' "Wine, Terroir and Climate Change". More involved than the
  # indices above (day-length and diurnal-range corrections, monthly cap) -
  # not yet implemented.
  # BEDD <- ...

  # CI - Cool Night Index (mean min. temperature of the ripening month, Sep)
  CI <- sdat[, TNInd[9]]

  # DI - Riou Dryness Index: soil water balance over Apr-Sep, starting from
  # a 200 mm reserve, replenished by rainfall and depleted by soil
  # evaporation (es) and vine transpiration (tv). kdi is the fraction of PET
  # attributed to vine transpiration each month (increasing canopy cover).
  jpm5 <- sdat[, RRInd[4:9]] / 5
  jpm5test <- jpm5 > ndaysmat[, 4:9]
  jpm5[jpm5test] <- ndaysmat[, 4:9][jpm5test]

  kdi <- c(0.1, 0.3, 0.5, 0.5, 0.5, 0.5)
  es <- t(t(sdat[, ETInd[4:9]] / ndaysmat[, 4:9] * jpm5) * (1 - kdi))
  tv <- t(t(sdat[, ETInd[4:9]]) * kdi)
  rr <- sdat[, RRInd[4:9]]
  DI <- NULL
  for (i in 1:ncol(rr)) {
    if (i == 1) {
      DI <- 200 + rr[, i] - es[, i] - tv[, i]
      DI[DI > 200] <- 200
    } else {
      DI <- DI + rr[, i] - es[, i] - tv[, i]
      DI[DI > 200] <- 200
    }
  }

  # WWB - winter water balance (Oct-Mar), continuing from DI, capped at
  # 200 mm; Q accumulates the monthly surplus (runoff) above that cap.
  WWB <- NULL
  Q <- NULL
  for (i in 1:6) {
    if (i == 1) {
      WWB <- DI + winter_rr[, i] - winter_et[, i]
      Q <- WWB - 200
      Q[Q < 0] <- 0
      WWB[WWB > 200] <- 200
      WWB[WWB < 0] <- 0
    } else {
      WWB <- WWB + winter_rr[, i] - winter_et[, i]
      Qm <- WWB - 200
      Qm[Qm < 0] <- 0
      Q <- Q + Qm
      WWB[WWB > 200] <- 200
      WWB[WWB < 0] <- 0
    }
  }

  # Koppen-Geiger climate classification
  KOP <- koppen.fast(teta = sdat[, TMInd], rr = sdat[, RRInd])

  # Annual and seasonal aggregates
  ANT <- rowMeans(sdat[, TMInd])
  ANP <- rowSums(sdat[, RRInd])
  GST <- rowMeans(sdat[, TMInd[4:10]])
  GSP <- rowSums(sdat[, RRInd[4:10]])
  GST49 <- rowMeans(sdat[, TMInd[4:9]])
  GSP49 <- rowSums(sdat[, RRInd[4:9]])
  DORT <- rowMeans(sdat[, TMInd[c(10:12, 1:3)]])
  DST <- rowMeans(sdat[, TMInd[4:7]])
  DSP <- rowSums(sdat[, RRInd[4:7]])
  RipeT <- sdat[, TMInd[9]]
  RipeP <- sdat[, RRInd[9]]
  FST <- rowMeans(sdat[, TMInd[7:9]])
  FSP <- rowSums(sdat[, RRInd[7:9]])
  SFR <- Rfast::rowMins(as.matrix(sdat[, TNInd[4:6]]), value = TRUE)
  WCR <- Rfast::rowMins(as.matrix(sdat[, TNInd]), value = TRUE)
  HSI <- Rfast::rowMaxs(as.matrix(sdat[, TXInd]), value = TRUE)

  data.frame(HI, CI, DI, Q, WWB, WI, ANT, ANP, KOPPEN = KOP,
             GST, GSP, GST49, GSP49, DST, DSP, RipeT, RipeP,
             FST, FSP, DORT, SFR, WCR, HSI)
}
