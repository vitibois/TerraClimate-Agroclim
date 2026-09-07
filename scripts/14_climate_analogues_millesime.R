# 2026-09-03 B. Bois (with Claude Code)
#
# Climate analogues of a vintage (phase 1, site level).
#
# Question: in which of the world's wine regions do we frequently find
# vintages whose vegetative-season climate closely matches that of one local
# site in one given year? Method and full statistical rationale (why
# detrending, why normal scores, why a pooled intra-site covariance metric,
# why three concurrent estimators, mandatory diagnostics) are documented in
# analogues_millesime.md at the project root -- read that file first. This
# script implements it against two concrete datasets:
#
#   WLD -- world vineyards, monthly climate, individual years: the .fst
#          written by scripts/07_extract_vgdb_point_data.R's "individual
#          years" section (one row per VGDB point per year, 48 monthly
#          tmin/tmax/ppt/pet columns).
#   LOC -- one local site's daily or monthly climate time series (its own
#          CSV, e.g. a Davis weather-station export), independent of the
#          rest of this pipeline.
#
# Both hemispheres are handled: for a southern-hemisphere point, vintage N's
# vegetative season is taken as Oct(N-1)-Feb(N) and re-mapped into the same
# slot columns as the northern-hemisphere Apr-Aug window (see
# season_calendar()) -- see analogues_millesime.md section 1.2.

library(data.table)
library(fst)
library(ggplot2)
library(sf)

source("E:/Benjamin/_Recherches_Dijon/_AAA_TerraClimate/Rscripts/R_Project_TerraClimate_Agroclim/R/config.R")

## ---- 0. configuration ----------------------------------------------------

# -- WLD: world vineyard monthly climate, individual years (script 07 output)
wld_path <- file.path(data_root, "Extraction_TerraClimatPoints", paste0("VGDB_v", vgdb_version),
                       paste0("MonthlyClimate_TerraClimate_VGDB_Pts_v", vgdb_version,
                              "_IndividualYears2001_2020.fst"))

# -- LOC: local site climate time series (one file, daily or monthly)
loc_path <- "E:/AAA_bbdocs/Professionnel/Conseil/2023-2025_ClarenceDillon/Data_Clim/data_ready_davis/quotidien/Pessac Haut Lafue Nord_20000101_20260901.csv"
loc_freq <- "daily"   # "daily" (aggregated to monthly below) or "monthly"
loc_cols <- list(date = "DATE", tmin = "TN", tmax = "TX", ppt = "RR", pet = "ET0")
# Hemisphere convention for LOC's vintage year; NA = auto-detect from a LAT
# column in loc_path if present, otherwise assume "N".
loc_hemisphere <- NA

target_year   <- 2026L   # the LOC vintage to characterise
months_season <- 4:8     # vegetative-season window, northern-hemisphere month numbers (Apr-Aug)

# -- optional: a grep pattern to identify LOC's own region among WLD$region,
# for the "control site against itself" calibration diagnostic (8.1). Set to
# NA to skip (LOC here is a private site, not one of the VGDB points).
loc_region_grep <- NA

K       <- 6L    # PCA truncation
lambda  <- NA    # regularised-Gaussian shrinkage; NA = choose by leave-one-vintage-out CV
m_knn   <- 5L    # kNN rank
q_count <- 0.50  # quantile of chisq(K) defining the "close vintage" radius r

run_bootstrap <- TRUE
boot_B <- 30L
boot_N <- 200L

n_top_regions_plot <- 10L

dir_out <- file.path(data_root, "Analogues_Millesime")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
site_label <- tools::file_path_sans_ext(basename(loc_path))
out_prefix <- file.path(dir_out, paste0("AnaloguesMillesime_", site_label, "_", target_year))

set.seed(1)

# ggsave() can hit a transient "cannot open file" error on Windows right
# after a large pdf() device has just closed in the same folder (antivirus/
# indexer briefly locking the newly written file) -- retry a few times.
safe_ggsave <- function(filename, plot, ..., max_tries = 10, wait = 3) {
  for (i in seq_len(max_tries)) {
    ok <- tryCatch({ ggsave(filename, plot, ...); TRUE },
                    error = function(e) { Sys.sleep(wait * i); FALSE })
    if (ok) return(invisible())
  }
  ggsave(filename, plot, ...)
}

# Same retry logic, for a multi-page pdf() built from a list of plots
# (ggplot objects need an explicit print() per page). Larger/heavier PDFs
# (e.g. the multi-page world map) can keep the lock longer, hence the
# growing backoff (wait * i) rather than a fixed delay.
safe_pdf_pages <- function(filename, plot_list, width, height, max_tries = 10, wait = 3) {
  do_pages <- function() { pdf(filename, width = width, height = height)
                            for (pl in plot_list) print(pl)
                            dev.off() }
  for (i in seq_len(max_tries)) {
    ok <- tryCatch({ do_pages(); TRUE },
                    error = function(e) { try(dev.off(), silent = TRUE); Sys.sleep(wait * i); FALSE })
    if (ok) return(invisible())
  }
  do_pages()
}

vars_map <- c(tmin = "Tmin", tmax = "Tmax", ppt = "P", pet = "ET0")
slot_lab <- sprintf("%02d", months_season)
VARS <- as.vector(outer(unname(vars_map), slot_lab, paste, sep = "_"))
pcols <- paste0("PC", 1:K)
wcols <- paste0("W", 1:K)

## ---- 1. hemisphere-aware season alignment ---------------------------------
# For a raw calendar year, which (calendar month, vintage-year offset) feeds
# each requested slot. vintage_year = raw_year - year_delta. Northern
# hemisphere: identity. Southern hemisphere: the window is shifted 6 months
# and split into at most two raw-year groups (see analogues_millesime.md
# section 1.2 for why Apr-Aug <-> Oct-Feb).

season_calendar <- function(months_season, hemisphere) {
  if (hemisphere == "N") {
    data.table(slot = seq_along(months_season), slot_lab = sprintf("%02d", months_season),
               month = months_season, year_delta = 0L)
  } else if (hemisphere == "S") {
    shifted <- ((months_season - 1L + 6L) %% 12L) + 1L
    delta   <- ifelse((months_season + 6L) <= 12L, 1L, 0L)  # raw_year = vintage_year - delta
    data.table(slot = seq_along(months_season), slot_lab = sprintf("%02d", months_season),
               month = shifted, year_delta = -delta)
  } else stop("hemisphere must be \"N\" or \"S\"")
}

## ---- 2. load and reshape WLD ----------------------------------------------

message("loading WLD: ", wld_path)
WLD <- as.data.table(read_fst(wld_path))
WLD[, site := .GRP, by = .(x, y)]

build_wld_season <- function(wld, months_season) {
  make_part <- function(sub, hemi) {
    cal <- season_calendar(months_season, hemi)
    parts <- lapply(split(cal, cal$year_delta), function(rows) {
      yd <- rows$year_delta[1]
      p <- sub[, .(site, region = CN_REG, x, y, year = year - yd)]
      for (i in seq_len(nrow(rows))) {
        for (v in names(vars_map)) {
          set(p, j = paste0(vars_map[[v]], "_", rows$slot_lab[i]),
              value = sub[[paste0(v, "_", rows$month[i])]])
        }
      }
      p
    })
    Reduce(function(a, b) merge(a, b, by = c("site", "region", "x", "y", "year")), parts)
  }
  rbind(make_part(wld[y >= 0], "N"), make_part(wld[y < 0], "S"))
}

D <- build_wld_season(WLD, months_season)
rm(WLD); gc()

# A handful of coastal pixels can carry NA values (nodata neighbour cells);
# drop any site with an incomplete vintage rather than let it propagate NAs
# through the PCA/covariance steps.
na_sites <- unique(D[!complete.cases(D[, ..VARS]), site])
if (length(na_sites) > 0) {
  message("dropping ", length(na_sites), " WLD site(s) with NA in the season window (coastal nodata)")
  D <- D[!site %in% na_sites]
}

stopifnot(all(VARS %in% names(D)), !anyNA(D[, ..VARS]))
nn <- D[, .N, by = site]$N
if (length(unique(nn)) != 1L)
  message("note: number of vintages per site is not constant (range ",
          paste(range(nn), collapse = "-"),
          ") -- expected near the reference period's ends for southern-hemisphere sites")
message("WLD sites: ", uniqueN(D$site), " | regions: ", uniqueN(D$region),
        " | rows: ", format(nrow(D), big.mark = " "))

## ---- 3. load and reshape LOC ----------------------------------------------

message("loading LOC: ", loc_path)
loc_raw <- fread(loc_path)

if (is.na(loc_hemisphere)) {
  lat_col <- grep("^lat$", names(loc_raw), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(lat_col)) {
    loc_hemisphere <- if (loc_raw[[lat_col]][1] < 0) "S" else "N"
    message("LOC hemisphere auto-detected from ", lat_col, ": ", loc_hemisphere)
  } else {
    loc_hemisphere <- "N"
    message("no LAT column found in LOC; assuming hemisphere = \"N\"")
  }
}

loc_dates <- as.IDate(loc_raw[[loc_cols$date]])
loc_raw[, `:=`(year = year(loc_dates), month = month(loc_dates))]

if (loc_freq == "daily") {
  loc_monthly <- loc_raw[, .(
    tmin   = mean(get(loc_cols$tmin), na.rm = TRUE),
    tmax   = mean(get(loc_cols$tmax), na.rm = TRUE),
    ppt    = sum(get(loc_cols$ppt),  na.rm = TRUE),
    pet    = sum(get(loc_cols$pet),  na.rm = TRUE),
    n_days = .N
  ), by = .(year, month)]
  loc_monthly[, days_expected := lubridate::days_in_month(as.Date(paste(year, month, 1, sep = "-")))]
  incomplete <- loc_monthly[n_days < 0.9 * days_expected]
  if (nrow(incomplete) > 0)
    message("note: ", nrow(incomplete), " LOC month(s) have <90% daily coverage (gaps in the raw series)")
} else if (loc_freq == "monthly") {
  loc_monthly <- loc_raw[, .(year, month,
                              tmin = get(loc_cols$tmin), tmax = get(loc_cols$tmax),
                              ppt = get(loc_cols$ppt), pet = get(loc_cols$pet))]
} else stop("loc_freq must be \"daily\" or \"monthly\"")

build_loc_season_vector <- function(loc_monthly, months_season, hemisphere, target_year) {
  cal <- season_calendar(months_season, hemisphere)
  cal[, cal_year := target_year - year_delta]
  m <- merge(cal, loc_monthly, by.x = c("cal_year", "month"), by.y = c("year", "month"), all.x = TRUE)
  setorder(m, slot)
  if (anyNA(m[, .(tmin, tmax, ppt, pet)]))
    stop("LOC is missing one or more months needed for the ", target_year, " season -- check loc_path's date coverage")
  yb <- unlist(lapply(seq_len(nrow(m)), function(i)
    setNames(as.list(m[i, .(tmin, tmax, ppt, pet)]), paste0(vars_map, "_", m$slot_lab[i]))
  ))
  yb[VARS]
}

yb <- build_loc_season_vector(loc_monthly, months_season, loc_hemisphere, target_year)
message("LOC target profile (", site_label, ", ", target_year, "):")
print(yb)

## ---- 4. detrending and rescaling to target_year, normal scores -----------
# See analogues_millesime.md sections 2.1-2.2. Detrending is done before the
# normal-score transform; a physical-units copy is kept for the interpretable
# monthly-profile plot in section 8.

TARGET_YEAR <- target_year

D[, (VARS) := {
  yy  <- year
  ww  <- yy - mean(yy)
  Sww <- sum(ww^2)
  lapply(.SD, function(x) x - (sum(ww * x) / Sww) * (yy - TARGET_YEAR))
}, by = site, .SDcols = VARS]

D_phys  <- copy(D)   # detrended, physical units -- kept for plotting only
yb_phys <- yb

EPS <- 1e-6
for (v in VARS) {
  Fv    <- ecdf(D[[v]])
  yb[v] <- qnorm(pmin(pmax(Fv(yb[v]), EPS), 1 - EPS))
  set(D, j = v, value = qnorm(pmin(pmax(Fv(D[[v]]), EPS), 1 - EPS)))
}

## ---- 5. PCA space and pooled intra-site covariance metric -----------------

X   <- as.matrix(D[, ..VARS])
pca <- prcomp(X, center = TRUE, scale. = TRUE)
message("K = ", K, " PCs | variance retained = ",
        round(100 * sum(pca$sdev[1:K]^2) / sum(pca$sdev^2), 1), " %")

Z  <- pca$x[, 1:K, drop = FALSE]
zb <- drop(predict(pca, matrix(yb, nrow = 1, dimnames = list(NULL, VARS))))[1:K]
rm(X); gc()

DT <- data.table(site = D$site, region = D$region, year = D$year)
DT[, (pcols) := as.data.table(Z)]

ANO <- copy(DT)
ANO[, (pcols) := lapply(.SD, function(x) x - mean(x)), by = site, .SDcols = pcols]
df_pool <- nrow(ANO) - uniqueN(ANO$site)
W <- crossprod(as.matrix(ANO[, ..pcols])) / df_pool
message("pooled degrees of freedom: ", format(df_pool, big.mark = " "),
        " | condition number of W: ", round(kappa(W), 2))
rm(ANO); gc()

Rw <- chol(W)
Zw  <- t(backsolve(Rw, t(Z), transpose = TRUE))
zbw <- backsolve(Rw, zb, transpose = TRUE)

DTw <- data.table(site = D$site, region = D$region, year = D$year)
DTw[, (wcols) := as.data.table(Zw)]
DTw[, d := sqrt(rowSums(sweep(Zw, 2, zbw, "-")^2))]

## ---- 6. three estimators ---------------------------------------------------

# -- 6.1 main: regularised-Gaussian log-likelihood --------------------------
loglik_lovo <- function(M, lambda) {
  n <- nrow(M); s <- 0
  for (i in seq_len(n)) {
    Mi <- M[-i, , drop = FALSE]
    S  <- lambda * W + (1 - lambda) * cov(Mi)
    R  <- chol(S)
    x  <- backsolve(R, M[i, ] - colMeans(Mi), transpose = TRUE)
    s  <- s - sum(log(diag(R))) - 0.5 * sum(x^2)
  }
  s
}

if (is.na(lambda)) {
  ech  <- sample(unique(DT$site), min(400L, uniqueN(DT$site)))
  grid <- seq(0.1, 1, by = 0.1)
  LL <- vapply(ech, function(s) {
    M <- as.matrix(DT[site == s, ..pcols])
    vapply(grid, function(l) loglik_lovo(M, l), numeric(1))
  }, numeric(length(grid)))
  lambda <- grid[which.max(rowSums(LL))]
  message("lambda chosen by leave-one-vintage-out CV: ", lambda)
} else {
  message("lambda fixed by config: ", lambda)
}

sc_gau <- DT[, {
  M  <- as.matrix(.SD)
  S  <- lambda * W + (1 - lambda) * cov(M)
  R  <- chol(S)
  x  <- backsolve(R, zb - colMeans(M), transpose = TRUE)
  .(score_gau = -sum(log(diag(R))) - 0.5 * sum(x^2) - 0.5 * K * log(2 * pi))
}, by = .(site, region), .SDcols = pcols]

# -- 6.2 control: Gaussian-kernel density -------------------------------------
h <- sqrt(qchisq(0.50, df = K)) / 2
sc_kde <- DTw[, .(score_kde = mean(exp(-0.5 * (d / h)^2))), by = .(site, region)]

# -- 6.3 control: adaptive kNN density ----------------------------------------
sc_knn <- DTw[, {
  dm <- sort(d)[m_knn]
  .(score_knn = log(m_knn) - log(.N) - K * log(dm), d_m = dm)
}, by = .(site, region)]

# -- 6.4 communication statistic: analogue count -----------------------------
r <- sqrt(qchisq(q_count, df = K))
sc_cnt <- DTw[, .(n_proches = sum(d <= r), n_tot = .N,
                   d_min = min(d), an_proche = year[which.min(d)]),
              by = .(site, region)]

res <- Reduce(function(a, b) merge(a, b, by = c("site", "region")),
              list(sc_gau, sc_kde, sc_knn, sc_cnt))
res[, `:=`(rg_gau = frank(-score_gau), rg_kde = frank(-score_kde), rg_knn = frank(-score_knn))]
setorder(res, rg_gau)

message("\ntop 30 analogue sites (ranked by regularised-Gaussian log-likelihood):")
print(res[1:30, .(site, region, rg_gau, rg_kde, rg_knn, n_proches, n_tot, d_min, an_proche)])

## ---- 7. monthly profile: the month that disqualifies ----------------------

profil_mensuel <- function(mo) {
  v   <- paste0(unname(vars_map), "_", sprintf("%02d", mo))
  Xm  <- scale(as.matrix(D[, ..v]))
  ctr <- attr(Xm, "scaled:center"); scl <- attr(Xm, "scaled:scale")

  Am <- data.table(site = D$site)
  Am[, (v) := as.data.table(Xm)]
  Am[, (v) := lapply(.SD, function(x) x - mean(x)), by = site, .SDcols = v]
  Wm <- crossprod(as.matrix(Am[, ..v])) / (nrow(Am) - uniqueN(Am$site))
  Rm <- chol(Wm)

  xb <- (yb[v] - ctr) / scl
  Xw <- t(backsolve(Rm, t(Xm), transpose = TRUE))
  xw <- backsolve(Rm, xb, transpose = TRUE)
  sqrt(rowSums(sweep(Xw, 2, xw, "-")^2))
}

Dm <- vapply(months_season, profil_mensuel, numeric(nrow(D)))
colnames(Dm) <- slot_lab

prof <- data.table(site = D$site, region = D$region, Dm)
prof_site <- prof[, lapply(.SD, median), by = .(site, region), .SDcols = slot_lab]
prof_site[, mois_pire := month.abb[months_season][apply(.SD, 1, which.max)], .SDcols = slot_lab]

res <- merge(res, prof_site, by = c("site", "region"))
setorder(res, rg_gau)

fwrite(res, paste0(out_prefix, "_table.csv"))
message("wrote ", paste0(out_prefix, "_table.csv"))

## ---- 8. diagnostics ---------------------------------------------------------

pdf(paste0(out_prefix, "_diagnostics.pdf"), width = 9, height = 6)

# -- 8.1 internal control: LOC's own region against itself, if identifiable --
if (!is.na(loc_region_grep)) {
  own <- res[grepl(loc_region_grep, region, ignore.case = TRUE)]
  if (nrow(own) > 0) {
    message("\ncontrol -- LOC's own region (\"", loc_region_grep, "\") n_proches summary:")
    print(summary(own$n_proches))
  } else {
    message("loc_region_grep \"", loc_region_grep, "\" matched no WLD region")
  }
}

# -- 8.2 bootstrap over vintages: is the ranking reproducible? ---------------
if (run_bootstrap) {
  boot_top <- function(B = boot_B, N = boot_N) {
    ref <- res[order(rg_gau)][1:N, site]
    ov  <- numeric(B)
    for (b in seq_len(B)) {
      idx <- DT[, .I[sample(.N, .N, replace = TRUE)], by = site]$V1
      Db  <- DT[idx]
      sb  <- Db[, {
        M <- as.matrix(.SD)
        S <- lambda * W + (1 - lambda) * cov(M)
        R <- chol(S)
        x <- backsolve(R, zb - colMeans(M), transpose = TRUE)
        .(sc = -sum(log(diag(R))) - 0.5 * sum(x^2))
      }, by = site, .SDcols = pcols]
      ov[b] <- length(intersect(ref, sb[order(-sc)][1:N, site])) / N
    }
    ov
  }
  ov <- boot_top()
  message("\nbootstrap overlap of top-", boot_N, " (B=", boot_B, "): ",
          round(100 * mean(ov)), "% [", round(100 * min(ov)), "-", round(100 * max(ov)), "]")
  if (mean(ov) < 0.5)
    message("WARNING: bootstrap overlap below 50% -- the site-level ranking is not reproducible; ",
            "treat results as exploratory only (see analogues_millesime.md section 9, phase 2)")
}

# -- 8.3 concordance between estimators --------------------------------------
topN <- function(col, N = 200L) res[order(res[[col]])][1:N, site]
message("\nconcordance of top-200 between estimators:")
for (p in list(c("rg_gau", "rg_knn"), c("rg_gau", "rg_kde"), c("rg_kde", "rg_knn")))
  message(p[1], " / ", p[2], " : ", length(intersect(topN(p[1]), topN(p[2]))), "/200")

# -- 8.5 spatial coherence within region -------------------------------------
coh <- res[, .(n = .N, rg_med = median(rg_gau), rg_iqr = IQR(rg_gau),
               part_top5pct = mean(rg_gau <= 0.05 * nrow(res))), by = region][order(rg_med)]
# log = "x" only: single-site regions have rg_iqr = 0, which a log y-axis
# cannot show.
print(plot(coh$rg_med, coh$rg_iqr, log = "x",
           xlab = "median rank", ylab = "rank IQR",
           main = "Spatial coherence within region (low-left & high part_top5pct = credible)"))

# -- 8.6 Gaussian calibration --------------------------------------------------
qqplot(qchisq(ppoints(1e5), df = K), sample(rowSums(Zw^2), min(1e5, nrow(Zw))),
       xlab = "theoretical chi-squared", ylab = "observed squared distances",
       main = "Gaussian calibration")
abline(0, 1, col = "red")

dev.off()
message("wrote ", paste0(out_prefix, "_diagnostics.pdf"))

## ---- 9. monthly profile plot: top regions vs LOC ---------------------------

top_regions <- coh[order(rg_med)][1:n_top_regions_plot, region]
top_regions5 <- top_regions[1:5]
# Fixed palette/shapes, reused as-is for the map in section 10 so all plots agree.
region_colors <- setNames(scales::hue_pal()(length(top_regions)), top_regions)
region_shapes <- setNames(0:(length(top_regions) - 1), top_regions)

phys_long <- melt(D_phys[region %in% top_regions, c("site", "region", "year", VARS), with = FALSE],
                   id.vars = c("site", "region", "year"), variable.name = "col")
phys_long[, `:=`(variable = sub("_.*$", "", col), slot = sub("^.*_", "", col))]
phys_long[, month_lab := factor(slot, levels = slot_lab,
                                 labels = month.abb[months_season])]
phys_long[, region := factor(region, levels = top_regions)]
prof_region <- phys_long[, .(value = median(value)), by = .(region, variable, month_lab)]

yb_long <- data.table(variable = sub("_.*$", "", VARS), slot = sub("^.*_", "", VARS), value = yb_phys[VARS])
yb_long[, month_lab := factor(slot, levels = slot_lab, labels = month.abb[months_season])]

# regions_n: which regions this page shows (top 10 or top 5); colour and
# shape both vary by region, reusing the fixed palette/shapes above so the
# same region always looks the same across pages and across the map.
build_monthly_plot <- function(regions_n) {
  pr <- prof_region[region %in% regions_n]
  pr[, region := factor(region, levels = regions_n)]

  ggplot(pr, aes(month_lab, value, group = region, color = region, shape = region)) +
    geom_line(alpha = 0.6) +
    geom_point(alpha = 0.6, size = 2) +
    geom_line(data = yb_long, aes(month_lab, value, group = 1), color = "black",
              linetype = "dashed", linewidth = 0.6, inherit.aes = FALSE) +
    geom_point(data = yb_long, aes(month_lab, value, group = 1), color = "black",
               size = 1.5, inherit.aes = FALSE) +
    facet_wrap(~variable, scales = "free_y") +
    scale_color_manual(values = region_colors[regions_n]) +
    scale_shape_manual(values = region_shapes[regions_n]) +
    labs(x = NULL, y = NULL, color = "Region", shape = "Region",
         title = paste0("Monthly profile: top ", length(regions_n), " analogue regions vs ",
                         site_label, " ", target_year),
         subtitle = "Black dashed = target vintage; coloured = median of each region's detrended vintages") +
    theme_minimal()
}

profile_pages <- list(build_monthly_plot(top_regions), build_monthly_plot(top_regions5))
safe_pdf_pages(paste0(out_prefix, "_monthly_profile.pdf"), profile_pages, width = 10, height = 7)
message("wrote ", paste0(out_prefix, "_monthly_profile.pdf"),
        " (2 pages: top", length(top_regions), ", top5)")

## ---- 10. world map of WLD points and analogue regions ---------------------
# Grey world boundaries; every WLD site as a small transparent black dot
# (site is already keyed on unique (x,y), so this is spatial-duplicate-free);
# the target site as a large white/black circle; the analogue regions (same
# palette as the monthly-profile plot above) as bigger coloured circles at
# the mean coordinates of their member sites. Four pages: world / zoomed-in,
# each for the top 10 and top 5 analogue regions.

world_bounds <- st_read(file.path(sig_root, "world/World_AdmBoundaries_Countries/world-administrative-boundaries.shp"),
                         quiet = TRUE)

site_coords <- unique(D[, .(site, region, x, y)])

region_coords <- site_coords[region %in% top_regions, .(x = mean(x), y = mean(y)), by = region]

loc_lon_col <- grep("^lon$", names(loc_raw), ignore.case = TRUE, value = TRUE)[1]
loc_lat_col <- grep("^lat$", names(loc_raw), ignore.case = TRUE, value = TRUE)[1]
if (is.na(loc_lon_col) || is.na(loc_lat_col))
  stop("no LON/LAT column found in LOC -- cannot place the target site on the map")
loc_coord <- data.table(x = loc_raw[[loc_lon_col]][1], y = loc_raw[[loc_lat_col]][1])

# regions_n: which regions to show: bbox = NULL for a world view, or a
# data.table(x, y) point set (mean region coordinates + LOC) to crop to.
build_analogue_map <- function(regions_n, bbox_pts = NULL) {
  rc_n <- region_coords[region %in% regions_n]
  rc_n[, region := factor(region, levels = regions_n)]

  p <- ggplot() +
    geom_sf(data = world_bounds, fill = "grey93", color = "grey65", linewidth = 0.1) +
    geom_point(data = site_coords, aes(x, y), color = "black", alpha = 0.2, shape = 20, size = 1) +
    geom_point(data = rc_n, aes(x, y, fill = region), shape = 21, color = "black",
               size = 4, stroke = 0.3) +
    geom_point(data = loc_coord, aes(x, y), shape = 21, fill = "white", color = "black",
               size = 6, stroke = 1.2) +
    scale_fill_manual(values = region_colors[regions_n], name = "Analogue region") +
    labs(x = NULL, y = NULL,
         title = paste0("WLD vineyard points and top ", length(regions_n), " analogue regions"),
         subtitle = paste0("White circle = ", site_label, " (target site)")) +
    theme_minimal()

  if (is.null(bbox_pts)) {
    p + coord_sf(ylim = c(-60, 85), expand = FALSE)
  } else {
    xr <- range(bbox_pts$x); yr <- range(bbox_pts$y)
    xpad <- max(diff(xr) * 0.15, 2); ypad <- max(diff(yr) * 0.15, 2)
    p + coord_sf(xlim = xr + c(-xpad, xpad), ylim = yr + c(-ypad, ypad), expand = FALSE)
  }
}

bbox10 <- rbind(region_coords[region %in% top_regions, .(x, y)], loc_coord)
bbox5  <- rbind(region_coords[region %in% top_regions5, .(x, y)], loc_coord)

map_pages <- list(
  build_analogue_map(top_regions),
  build_analogue_map(top_regions, bbox10),
  build_analogue_map(top_regions5),
  build_analogue_map(top_regions5, bbox5)
)

safe_pdf_pages(paste0(out_prefix, "_map.pdf"), map_pages, width = 12, height = 7)
message("wrote ", paste0(out_prefix, "_map.pdf"), " (4 pages: world/top", length(top_regions),
        ", zoom/top", length(top_regions), ", world/top5, zoom/top5)")

## ---- 11. parameter sensitivity (manual, see analogues_millesime.md 8.4) ---
# Not automated here (re-running sections 5-6 several times is costly on
# ~925k rows). To check robustness, manually rerun with K in c(4,6,8),
# lambda in c(0.3, <CV value>, 0.9), m_knn in c(3,5,8), and q_count in
# c(0.35, 0.50, 0.65), and compare topN() overlaps as in section 8.3.
