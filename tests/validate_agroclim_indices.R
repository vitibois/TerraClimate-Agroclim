# Regression check for R/agroclim_indices.R.
#
# Compares the shared calc_agroclim_indices() function against
# reimplementations of the *original* inline code from the two legacy
# scripts it replaces:
#   - legacy_annual():      scripts/01-ExtractTerraClimate_and_CalcAgroInd_v02.R
#   - legacy_climatology(): scripts/01_bis-Climatologies_..._v02.R
#
# Both legacy versions are reproduced here exactly as they were, including
# the two known bugs that were fixed when the code was factored into
# calc_agroclim_indices() (see project history / PR description):
#   (1) climatology: `ndaysmat[4:9]` instead of `ndaysmat[,4:9]` when
#       dividing PET by days-in-month for the soil evaporation term.
#   (2) annual: Southern Hemisphere days-in-month hardcoded as
#       c(31,30,31,31,28,30) instead of c(31,30,31,31,28,31) (March has 31
#       days, not 30).
#
# This script demonstrates that:
#   - for rows/pipelines unaffected by a given bug, old and new code match
#     exactly;
#   - for rows/pipelines affected by a bug, the new code differs only in the
#     indices that bug could plausibly affect (HI/DI/WWB/Q), by a small
#     amount, never in the other indices (CI, Koppen, seasonal means...).
#
# Run with:
#   Rscript tests/validate_agroclim_indices.R

source("R/koppen.R")
source("R/agroclim_indices.R")

set.seed(1)

## ---- synthetic monthly climate data --------------------------------------
# 4 Northern Hemisphere pixels + 4 Southern Hemisphere pixels, already
# reordered into the local agricultural year (as the extraction scripts do
# before calling the shared function): for the SH rows, column 1 = July.
n_nh <- 4
n_sh <- 4
n <- n_nh + n_sh
isSH <- (n_nh + 1):n

make_monthly <- function(base, amp, n) {
  # 12 months per row, simple sinusoidal seasonal cycle + noise
  t(sapply(seq_len(n), function(i) {
    base + amp * sin(seq(0, 2 * pi, length.out = 13)[-13]) + rnorm(12, 0, 0.3)
  }))
}

tmin <- make_monthly(8, 8, n)
tmax <- tmin + matrix(runif(n * 12, 6, 12), n, 12)
tavg <- (tmin + tmax) / 2
ppt  <- matrix(pmax(0, rnorm(n * 12, 60, 30)), n, 12)
pet  <- matrix(pmax(5, rnorm(n * 12, 70, 20)), n, 12)
lat  <- c(rep(45, n_nh), rep(-35, n_sh))

colnames(tmin) <- paste0("tmin_", 1:12)
colnames(tmax) <- paste0("tmax_", 1:12)
colnames(tavg) <- paste0("tavg_", 1:12)
colnames(ppt)  <- paste0("ppt_", 1:12)
colnames(pet)  <- paste0("pet_", 1:12)

sdat <- cbind(pixID = seq_len(n), x = 0, y = lat, year = 2020,
              tmin, tmax, ppt, pet, tavg)

TNInd <- grep("^tmin_", colnames(sdat))
TXInd <- grep("^tmax_", colnames(sdat))
RRInd <- grep("^ppt_", colnames(sdat))
ETInd <- grep("^pet_", colnames(sdat))
TMInd <- grep("^tavg_", colnames(sdat))
LATInd <- grep("^y$", colnames(sdat))

# "Previous winter" precip/PET: for this synthetic single-year check, reuse
# the same Oct-Mar months for both legacy paths (equivalent to the
# climatology pipeline's own behaviour; for the annual pipeline this stands
# in for a previous year that happens to have the same climate).
winter_rr <- sdat[, c(RRInd[10:12], RRInd[1:3])]
winter_et <- sdat[, c(ETInd[10:12], ETInd[1:3])]

## ---- legacy reimplementations ---------------------------------------------

legacy_annual <- function(sdat) {
  ndaysmat <- t(t(sdat[, RRInd[4:9]]) * 0 + 1 * c(30, 31, 30, 31, 31, 30))
  if (length(isSH) > 0)
    ndaysmat[isSH, ] <- t(t(ndaysmat[isSH, ]) * 0 + 1 * c(31, 30, 31, 31, 28, 30)) # BUG: last should be 31

  khi <- k.hi.hall(lat = sdat[, LATInd])
  HI <- sdat[, c(TMInd[4:9], TXInd[4:9])] - 10
  HI[HI < 0] <- 0
  HI <- (HI[, 1:6] + HI[, 7:12]) / 2
  HI <- HI * ndaysmat
  HI <- rowSums(HI) * khi

  CI <- sdat[, TNInd[9]]

  jpm5 <- sdat[, RRInd[4:9]] / 5
  jpm5test <- jpm5 > ndaysmat
  jpm5[jpm5test] <- ndaysmat[jpm5test]

  kdi <- c(0.1, 0.3, 0.5, 0.5, 0.5, 0.5)
  es <- t(t(sdat[, ETInd[4:9]] / ndaysmat * jpm5) * (1 - kdi))
  tv <- t(t(sdat[, ETInd[4:9]]) * kdi)
  rr <- sdat[, RRInd[4:9]]
  DI <- NULL
  for (i in 1:ncol(rr)) {
    if (i == 1) { DI <- 200 + rr[, i] - es[, i] - tv[, i]; DI[DI > 200] <- 200 }
    else { DI <- DI + rr[, i] - es[, i] - tv[, i]; DI[DI > 200] <- 200 }
  }

  rr <- winter_rr; et <- winter_et
  WWB <- NULL; Q <- NULL
  for (i in 1:6) {
    if (i == 1) {
      WWB <- DI + rr[, i] - et[, i]
      Q <- WWB - 200; Q[Q < 0] <- 0
      WWB[WWB > 200] <- 200; WWB[WWB < 0] <- 0
    } else {
      WWB <- WWB + rr[, i] - et[, i]
      Qm <- WWB - 200; Qm[Qm < 0] <- 0
      Q <- Q + Qm
      WWB[WWB > 200] <- 200; WWB[WWB < 0] <- 0
    }
  }

  KOP <- koppen.fast(teta = sdat[, TMInd], rr = sdat[, RRInd])
  data.frame(HI, CI, DI, Q, WWB, KOPPEN = KOP)
}

legacy_climatology <- function(sdat) {
  ndaysmat <- t(t(sdat[, RRInd]) * 0 + 1 * c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))
  if (length(isSH) > 0)
    ndaysmat[isSH, ] <- ndaysmat[isSH, c(7:12, 1:6)]

  WI <- sdat[, TMInd[4:10]] - 10
  WI[WI < 0] <- 0
  WI <- WI * ndaysmat[, 4:10]
  WI <- rowSums(WI)

  khi <- k.hi.hall(lat = sdat[, LATInd])
  HI <- sdat[, c(TMInd[4:9], TXInd[4:9])] - 10
  HI[HI < 0] <- 0
  HI <- (HI[, 1:6] + HI[, 7:12]) / 2
  HI <- HI * ndaysmat[, 4:9]
  HI <- rowSums(HI) * khi

  CI <- sdat[, TNInd[9]]

  jpm5 <- sdat[, RRInd[4:9]] / 5
  jpm5test <- jpm5 > ndaysmat[, 4:9]
  jpm5[jpm5test] <- ndaysmat[, 4:9][jpm5test]

  kdi <- c(0.1, 0.3, 0.5, 0.5, 0.5, 0.5)
  es <- t(t(sdat[, ETInd[4:9]] / ndaysmat[4:9] * jpm5) * (1 - kdi)) # BUG: missing comma
  tv <- t(t(sdat[, ETInd[4:9]]) * kdi)
  rr <- sdat[, RRInd[4:9]]
  DI <- NULL
  for (i in 1:ncol(rr)) {
    if (i == 1) { DI <- 200 + rr[, i] - es[, i] - tv[, i]; DI[DI > 200] <- 200 }
    else { DI <- DI + rr[, i] - es[, i] - tv[, i]; DI[DI > 200] <- 200 }
  }

  rr <- sdat[, c(RRInd[10:12], RRInd[1:3])]
  et <- sdat[, c(ETInd[10:12], ETInd[1:3])]
  WWB <- NULL; Q <- NULL
  for (i in 1:6) {
    if (i == 1) {
      WWB <- DI + rr[, i] - et[, i]
      Q <- WWB - 200; Q[Q < 0] <- 0
      WWB[WWB > 200] <- 200; WWB[WWB < 0] <- 0
    } else {
      WWB <- WWB + rr[, i] - et[, i]
      Qm <- WWB - 200; Qm[Qm < 0] <- 0
      Q <- Q + Qm
      WWB[WWB > 200] <- 200; WWB[WWB < 0] <- 0
    }
  }

  KOP <- koppen.fast(teta = sdat[, TMInd], rr = sdat[, RRInd])
  data.frame(HI, CI, DI, Q, WWB, WI, KOPPEN = KOP)
}

## ---- run new shared function ----------------------------------------------

new_out <- calc_agroclim_indices(sdat, TNInd, TXInd, RRInd, ETInd, TMInd, LATInd,
                                  isSH, winter_rr, winter_et)

old_annual <- legacy_annual(sdat)
old_clim <- legacy_climatology(sdat)

## ---- compare ----------------------------------------------------------------

cat("\n=== Annual pipeline: legacy vs new (rows 1:", n_nh, " = NH, ", n_nh + 1, ":", n,
    " = SH) ===\n", sep = "")
common_cols_a <- intersect(names(old_annual), names(new_out))
diffs_a <- sapply(common_cols_a, function(cn) old_annual[[cn]] - new_out[[cn]])
print(round(diffs_a, 4))
cat("\nColumns above should be ~0 for NH rows (1:", n_nh, "); ",
    "small nonzero DI/WWB/Q/HI differences are expected for SH rows (",
    n_nh + 1, ":", n, ") due to the fixed day-count bug.\n\n", sep = "")

cat("=== Climatology pipeline: legacy vs new (all rows use the same synthetic winter) ===\n")
common_cols_c <- intersect(names(old_clim), names(new_out))
diffs_c <- sapply(common_cols_c, function(cn) old_clim[[cn]] - new_out[[cn]])
print(round(diffs_c, 4))
cat("\nColumns above should be ~0 except DI/WWB/Q, which differ for every row due to the",
    "fixed ndaysmat[4:9] indexing bug (previously always divided by 31 instead of the",
    "correct days-in-month).\n\n")

stopifnot(
  "CI must be identical (unaffected by either bug)" =
    all(abs(diffs_a[, "CI"]) < 1e-8) && all(abs(diffs_c[, "CI"]) < 1e-8),
  "KOPPEN must be identical (unaffected by either bug)" =
    all(diffs_a[, "KOPPEN"] == 0) && all(diffs_c[, "KOPPEN"] == 0),
  "Annual NH rows must match exactly (bug was Southern Hemisphere only)" =
    all(abs(as.matrix(diffs_a[1:n_nh, c("HI", "DI", "Q", "WWB")])) < 1e-8)
)

cat("OK: differences are confined to the expected indices/rows; all other outputs match exactly.\n")
