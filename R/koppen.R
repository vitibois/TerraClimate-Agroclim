# Koppen-Geiger climate classification (Peel, Finlayson & McMahon, 2007)
#
# Vendored from Benjamin Bois' personal function library (E:/R_functions/koppen.r),
# where a larger set of Koppen-related helpers is maintained. Only the fast,
# vectorised classifier actually used by this pipeline is kept here.
#
# Reference:
# Peel, M. C., Finlayson, B. L., & McMahon, T. A. (2007). Updated world map of
# the Koppen-Geiger climate classification. Hydrology and Earth System
# Sciences, 11(5), 1633-1644. https://doi.org/10.5194/hess-11-1633-2007

#' Classify pixels/rows into Koppen-Geiger climate types.
#'
#' Vectorised over rows: each row of `teta`/`rr` is one location (or
#' location/year), with columns 1:12 giving mean monthly temperature (`teta`,
#' degrees C) and monthly precipitation (`rr`, mm) for January to December.
#'
#' @param teta (n x 12) matrix of mean monthly temperatures.
#' @param rr (n x 12) matrix of monthly precipitation totals.
#' @param table.codes.only if `TRUE`, ignore `teta`/`rr` and just return the
#'   lookup table of classification codes (letter code, numeric code, "Peel"
#'   code and a suggested RGB colour) used to decode the output.
#' @param give.peel.code if `TRUE` (default), return the compact 1:30 "Peel"
#'   numeric code; if `FALSE`, return the 3-digit numeric code instead.
#'
#' @return an integer vector of classification codes (one per row of
#'   `teta`/`rr`), or the lookup table if `table.codes.only = TRUE`.
koppen.fast <- function(teta = "matrix.of.12.columns", rr = "matrix.of.12.columns",
                         table.codes.only = FALSE, give.peel.code = TRUE) {
  require(Rfast)
  # A classification code data.frame, that is returned optionally to better
  # read the classification codes
  letcode <- c("Af", "Am", "Aw",
               "BWh", "BWk", "BSh", "BSk",
               "Csa", "Csb", "Csc", "Cwa", "Cwb", "Cwc", "Cfa", "Cfb", "Cfc",
               "Dsa", "Dsb", "Dsc", "Dsd", "Dfa", "Dfb", "Dfc", "Dfd", "Dwa", "Dwb", "Dwc", "Dwd",
               "ET", "EF")
  numcode <- c(110, 120, 130,
               211, 212, 221, 222,
               311, 312, 313, 321, 322, 323, 331, 332, 333,
               411, 412, 413, 414, 421, 422, 423, 424, 431, 432, 433, 434,
               510, 520)
  peelnumcode <- 1:30
  # Colors
  r <- c(0, 0, 70,
         255, 255, 245, 255,
         255, 200, 150, 150, 100, 50, 200, 100, 50,
         255, 200, 150, 150, 170, 90, 75, 50, 0, 55, 0, 0,
         179, 102)
  g <- c(0, 120, 170,
         0, 150, 165, 220,
         255, 200, 150, 255, 200, 150, 255, 255, 200,
         0, 0, 50, 100, 175, 120, 80, 0, 255, 200, 125, 70,
         179, 102)
  b <- c(255, 255, 250,
         0, 150, 0, 100,
         0, 0, 0, 150, 100, 50, 80, 50, 0,
         255, 200, 150, 150, 255, 220, 180, 135, 255, 255, 125, 95,
         179, 102)
  if (table.codes.only == TRUE) {
    return(data.frame(letcode, numcode, peelnumcode, red = r, green = g, blue = b))
  } else {
    if (nrow(teta) != nrow(rr))
      stop("length of lats vector differ from the number of rows in teta or rr)")
    # Mean annual temperature and precipitations
    MAP <- rowSums(rr)
    MAT <- rowMeans(teta)
    # Driest month of the year
    Pdry <- Rfast::rowMins(rr, value = TRUE)
    # Driest months in summer (Psdry) and in winter (Pwdry)
    Psdry <- Rfast::rowMins(rr[, 4:9], value = TRUE)
    Pwdry <- Rfast::rowMins(rr[, c(1:3, 10:12)], value = TRUE)

    # Wettest months in summer (Pswet) and in winter (Pwwet)
    Pswet <- Rfast::rowMaxs(rr[, 4:9], value = TRUE)
    Pwwet <- Rfast::rowMaxs(rr[, c(1:3, 10:12)], value = TRUE)
    # Sum of summer and winter precipitations
    Ps <- rowSums(rr[, 4:9])
    Pw <- rowSums(rr[, c(1:3, 10:12)])

    # Threshold (Pt) for Arid climate classification
    Pt <- rep(NA, nrow(teta))
    Pt[Pw >= (0.7 * MAP)] <- 2 * MAT[Pw >= (0.7 * MAP)]
    Pt[Ps >= (0.7 * MAP)] <- 2 * MAT[Ps >= (0.7 * MAP)] + 28
    Pt[Pw < (0.7 * MAP) & Ps < (0.7 * MAP)] <-
      2 * MAT[Pw < (0.7 * MAP) & Ps < (0.7 * MAP)] + 14
    # Temperatures of the coldest (Tcold) and of the hottest (Thot) months
    Tcold <- Rfast::rowMins(teta, value = TRUE)
    Thot <- Rfast::rowMaxs(teta, value = TRUE)
    Tmon10 <- rowSums(teta > 10)

    # Classification.
    # The priority is given in the following order: B, then A, C
    # (with a priority given to w), D (with a priority given to w either)
    # and then E
    outcode <- vector("numeric", nrow(teta))
    outcode[] <- NA
    # BWh : Arid Desert Hot
    outcode[MAP < (5 * Pt) & MAT >= 18] <- 211
    # BWk : Arid Desert Cold
    outcode[is.na(outcode) & MAP < (5 * Pt) & MAT < 18] <- 212
    # BSh : Arid Steppe Hot
    outcode[is.na(outcode) & MAP >= (5 * Pt) & MAP < (10 * Pt) & MAT >= 18] <- 221
    # BSk : Arid Steppe Cold
    outcode[is.na(outcode) & MAP >= (5 * Pt) & MAP < (10 * Pt) & MAT < 18] <- 222

    # Af : Tropical Rainforest
    outcode[is.na(outcode) & Tcold >= 18 & Pdry >= 60] <- 110
    # Am : Tropical Monsoon
    outcode[is.na(outcode) & Tcold >= 18 & Pdry < 60 & Pdry >= (100 - MAP / 25)] <- 120
    # Aw : Tropical Savannah
    outcode[is.na(outcode) & Tcold >= 18 & Pdry < 60 & Pdry < (100 - MAP / 25)] <- 130

    # Cwa : Temperate with Dry Winter and Hot Summer
    outcode[is.na(outcode) & Thot > 10 & Tcold > 0 & Tcold < 18
            & Pwdry < (Pswet / 10)
            & Thot >= 22] <- 321
    # Cwb : Temperate with Dry Winter and Warm Summer
    outcode[is.na(outcode) & Thot > 10 & Tcold > 0 & Tcold < 18
            & Pwdry < (Pswet / 10)
            & Thot < 22 & Tmon10 >= 4] <- 322
    # Cwc : Temperate with Dry Winter and Cold Summer
    outcode[is.na(outcode) & Thot > 10 & Tcold > 0 & Tcold < 18
            & Pwdry < (Pswet / 10)
            & Thot < 22 & Tmon10 < 4 & Tmon10 >= 1] <- 323
    # Csa : Temperate with Dry and Hot Summer
    outcode[is.na(outcode) & Thot > 10 & Tcold > 0 & Tcold < 18
            & Psdry < 40 & Psdry < (Pwwet / 3)
            & Thot >= 22] <- 311
    # Csb : Temperate with Dry and Warm Summer
    outcode[is.na(outcode) & Thot > 10 & Tcold > 0 & Tcold < 18
            & Psdry < 40 & Psdry < (Pwwet / 3)
            & Thot < 22 & Tmon10 >= 4] <- 312
    # Csc : Temperate with Dry and Cold Summer
    outcode[is.na(outcode) & Thot > 10 & Tcold > 0 & Tcold < 18
            & Psdry < 40 & Psdry < (Pwwet / 3)
            & Thot < 22 & Tmon10 < 4 & Tmon10 >= 1] <- 313
    # Cfa : Temperate without dry season and with Hot Summer
    outcode[is.na(outcode) & Thot > 10 & Tcold > 0 & Tcold < 18
            & Thot >= 22] <- 331
    # Cfb : Temperate without dry season and with Warm Summer
    outcode[is.na(outcode) & Thot > 10 & Tcold > 0 & Tcold < 18
            & Thot < 22 & Tmon10 >= 4] <- 332
    # Cfc : Temperate without dry season and with Cold Summer
    outcode[is.na(outcode) & Thot > 10 & Tcold > 0 & Tcold < 18
            & Thot < 22 & Tmon10 < 4 & Tmon10 >= 1] <- 333

    # Dwa : Cold with Dry Winter and Hot Summer
    outcode[is.na(outcode) & Thot > 10 & Tcold <= 0
            & Pwdry < (Pswet / 10)
            & Thot >= 22] <- 421
    # Dwb : Cold with Dry Winter and Warm Summer
    outcode[is.na(outcode) & Thot > 10 & Tcold <= 0
            & Pwdry < (Pswet / 10)
            & Thot < 22 & Tmon10 >= 4] <- 422
    # Dwc : Cold with Dry Winter and Cold Summer
    outcode[is.na(outcode) & Thot > 10 & Tcold <= 0
            & Pwdry < (Pswet / 10)
            & Thot < 22 & Tmon10 < 4 & Tcold >= -38] <- 423
    # Dwd : Cold with Dry Winter and Very Cold Winter
    outcode[is.na(outcode) & Thot > 10 & Tcold <= 0
            & Pwdry < (Pswet / 10)
            & Thot < 22 & Tmon10 < 4 & Tcold < -38] <- 424

    # Dsa : Cold with Dry and Hot Summer
    outcode[is.na(outcode) & Thot > 10 & Tcold <= 0
            & Psdry < 40 & Psdry < (Pwwet / 3)
            & Thot >= 22] <- 411
    # Dsb : Cold with Dry and Warm Summer
    outcode[is.na(outcode) & Thot > 10 & Tcold <= 0
            & Psdry < 40 & Psdry < (Pwwet / 3)
            & Thot < 22 & Tmon10 >= 4] <- 412
    # Dsc : Cold with Dry and Cold Summer
    outcode[is.na(outcode) & Thot > 10 & Tcold <= 0
            & Psdry < 40 & Psdry < (Pwwet / 3)
            & Thot < 22 & Tmon10 < 4 & Tcold >= -38] <- 413
    # Dsd : Cold with Dry Summer and Very Cold Winter
    outcode[is.na(outcode) & Thot > 10 & Tcold <= 0
            & Psdry < 40 & Psdry < (Pwwet / 3)
            & Thot < 22 & Tmon10 < 4 & Tcold < -38] <- 414
    # Dfa : Cold without dry season and with Hot Summer
    outcode[is.na(outcode) & Thot > 10 & Tcold <= 0
            & Thot >= 22] <- 431
    # Dfb : Cold without dry season and with Warm Summer
    outcode[is.na(outcode) & Thot > 10 & Tcold <= 0
            & Thot < 22 & Tmon10 >= 4] <- 432
    # Dfc : Cold without dry season and with Cold Summer
    outcode[is.na(outcode) & Thot > 10 & Tcold <= 0
            & Thot < 22 & Tmon10 < 4 & Tcold >= -38] <- 433
    # Dsd : Cold without dry season and with Very Cold Winter
    outcode[is.na(outcode) & Thot > 10 & Tcold <= 0
            & Thot < 22 & Tmon10 < 4 & Tcold < -38] <- 434
    # ET Polar Tundra (note a slight change from Peel et al.: Thot < 10 in
    # Peel et al.; here Thot <= 10, to define Polar climates)
    outcode[is.na(outcode) & Thot <= 10 & Thot > 0] <- 510
    # EF Polar Frost (same note as above)
    outcode[is.na(outcode) & Thot <= 10 & Thot <= 0] <- 520
    if (give.peel.code == TRUE)
      outcode <- peelnumcode[match(outcode, numcode)]
    return(outcode)
  }
}
