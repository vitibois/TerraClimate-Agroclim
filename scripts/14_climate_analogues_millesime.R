# 2026-09-03 B. Bois (with Claude Code)
# Reworked 2026-09-04, B. Bois (with Claude, then Claude Code): migrated the
# main estimator from a regularised-Gaussian score to an empirical-Bayes NIW/
# Student predictive, coordinate-based station anchoring, source-delta
# calibration, sequential regional cohorts, and the map/report outputs below
# (see 00_note_migration.md for the full list of what changed and why).
#
# Climate analogues of a vintage (phase 1, site level).
#
# Question: in which of the world's wine regions do we frequently find
# vintages whose vegetative-season climate closely matches that of one local
# site in one given year?
#
# Method and full statistical rationale are documented in:
#   analogues_millesime.md        (methode, phase 1 parcellaire)
#   note_station_delta_carte.md   (ancrage station, delta, carto, rapport)
#   00_note_migration.md          (ce qui a change par rapport a la version precedente)
# Read those first.
#
# Datasets:
#   WLD -- world vineyards, monthly climate, individual years: the .fst
#          written by scripts/07_extract_vgdb_point_data.R ("individual
#          years" section). One row per VGDB point per year, 48 monthly
#          tmin/tmax/ppt/pet columns. ~47 000 points, ~700 wine regions,
#          density ~1 point per km2 -- these are POINTS with their own
#          coordinates, not a raster to extract.
#   LOC -- one local weather station's daily or monthly series (its own CSV,
#          e.g. a Davis export), independent of the rest of this pipeline.
#
# Both hemispheres are handled: for a southern-hemisphere point, vintage N's
# vegetative season is Oct(N-1)-Feb(N), re-mapped into the same slot columns
# as the northern-hemisphere Apr-Aug window (see season_calendar()).
#
# ---------------------------------------------------------------------------
# ESTIMATEUR PRINCIPAL : predictif bayesien a a priori empirique (NIW /
# Student). Terminologie imposee : "empirical Bayes" ou "predictif bayesien a
# a priori empirique", JAMAIS "analyse bayesienne" tout court. L'a priori W
# etant estime sur les memes donnees, les incertitudes a posteriori sont
# sous-estimees : le score est une STATISTIQUE DE CLASSEMENT ORDINAL, jamais
# une probabilite. Voir analogues_millesime.md sections 5.1.0 et 5.1.0bis.
# ---------------------------------------------------------------------------
#
# AVERTISSEMENT DE FOND (analogues_millesime.md 0.4) : 20 millesimes en
# dimension 6 est un regime statistique defavorable. Le top-1 parcellaire n'a
# aucune signification. Les diagnostics 8.2 (bootstrap) et 8.5 (coherence
# spatiale) conditionnent toute publication au niveau parcellaire.
#
# RAPPEL (note_station_delta_carte.md J) : fenetre avril-aout => analogues de
# SAISON VEGETATIVE uniquement. Aucune contrainte hivernale (besoins en
# froid, gel hivernal, dormance) n'entre dans le calcul.

library(data.table)
library(fst)
library(ggplot2)
library(sf)
library(terra)
library(tidyterra)
library(maptiles)
library(ggnewscale)

source("R/config.R")

## ---- 0. configuration ------------------------------------------------------

# -- WLD: world vineyard monthly climate, individual years (script 07 output)
wld_path <- file.path(data_root, "Extraction_TerraClimatPoints", paste0("VGDB_v", vgdb_version),
                      paste0("MonthlyClimate_TerraClimate_VGDB_Pts_v", vgdb_version,
                             "_IndividualYears2001_2025.fst"))

# -- colonnes reelles de WLD -- A VERIFIER avec names(WLD) avant execution
#    (note_station_delta_carte.md section B ; 00_note_migration.md section 5)
COL_REGION  <- "CN_REG"    # spec: "WINE_REGION"
COL_COUNTRY <- "CNT"       # optionnel : mis a NA si absent
COL_LON     <- "x"         # spec: "lon"
COL_LAT     <- "y"         # spec: "lat"

# -- LOC: local station climate time series (one file, daily or monthly)
loc_path <- "E:/AAA_bbdocs/Professionnel/Conseil/2023-2025_ClarenceDillon/Data_Clim/data_ready_davis/quotidien/Pessac Haut Lafue Nord_20000101_20260901.csv"
loc_freq <- "daily"   # "daily" (aggregated to monthly below) or "monthly"
loc_cols <- list(date = "DATE", tmin = "TN", tmax = "TX", ppt = "RR", pet = "ET0")

# -- ANCRAGE PAR COORDONNEES (remplace l'ancrage par nom de region) ----------
#    note_station_delta_carte.md section B. La region est DEDUITE du point
#    apparie ; aucun nom de region a saisir.
# NA = auto-detectees depuis les colonnes LON/LAT de loc_path (voir plus bas,
# section 0.1) -- fixer une valeur ici seulement pour forcer un point d'ancrage
# different de celui inscrit dans le fichier LOC.
LOC_LON  <- NA_real_   # longitude de la station, degres decimaux
LOC_LAT  <- NA_real_   # latitude
D_MAX_KM <- 15        # distance d'appariement maximale toleree (garde-fou
                      # contre une inversion lon/lat ou une erreur de signe)
VOIS_KM  <- 10        # rayon du voisinage de controle (B.2)

# Hemisphere convention for LOC's vintage year; NA = deduced from LOC_LAT.
loc_hemisphere <- NA

# -- SOURCE DU VECTEUR CIBLE (note_station_delta_carte.md section E) --------
#    "station" : millesime issu de LOC, corrige du delta de source (run B)
#    "pixel"   : millesime issu du point WLD apparie, sans delta      (run A)
#    Le protocole de validation a deux sources consiste a lancer ce script
#    deux fois, tout identique sauf cible_source, sur une annee presente dans
#    LES DEUX sources ET ORDINAIRE -- jamais sur le millesime "ovni" qu'on
#    cherche a diagnostiquer (E.3). Comparer ensuite les deux tableaux avec
#    comparer_runs() (fin de script).
cible_source <- "station"

target_year   <- 2026L   # the LOC vintage to characterise
months_season <- 4:8     # vegetative-season window, northern-hemisphere month numbers

# -- estimateurs -------------------------------------------------------------
K       <- 6L    # troncature ACP (regularisation, pas une commodite)
M0      <- 15L   # poids d'a priori en "millesimes fictifs" -- DECLARE, PAS
                 # OPTIMISE. Ne pas ecrire de validation croisee pour le
                 # choisir (analogues_millesime.md 5.1.3).
KAPPA0  <- 0.01  # a priori faible sur la position : on regularise la FORME
m_knn   <- 5L    # rang kNN
q_count <- 0.50  # quantile de chisq(K) definissant le rayon r de "millesime proche"

# -- calibration station / point apparie (sections C et D) -------------------
VARS4    <- c("Tmin", "Tmax", "P", "ET0")
ADD      <- c("Tmin", "Tmax")   # delta additif
MUL      <- c("P", "ET0")       # delta multiplicatif
P_MIN    <- 5                   # mm, plancher du denominateur du ratio
COMPL_MIN <- 0.95               # completude minimale d'un mois LOC journalier
N_AN_MIN <- 5L                  # annees communes minimales pour calibrer
delta_mode <- "moyenne"         # "moyenne" ou "median" : statistique du delta

# -- diagnostics et sorties --------------------------------------------------
run_brut_corr   <- TRUE   # rejouer le pipeline sur yb_brut et comparer (C/E.2)
run_bootstrap   <- TRUE   # 8.2 -- le diagnostic central au niveau parcellaire
run_sens_M0     <- TRUE   # 5.1.3 -- sensibilite a M0
save_continuite <- TRUE   # ecrit l'entree du script de controle 14b
boot_B <- 30L
boot_N <- 200L
M0_grid <- c(5L, 15L, 30L, 60L)

N_REG          <- 10L   # nombre de regions dans la synthese
n_zoom_regions <- 5L    # nombre de zooms cartographiques par region

# -- fond de carte (cartes du top-N et du rapport) ---------------------------
use_basemap_tiles <- TRUE     # tuiles OpenTopoMap ; repli automatique sur les
                              # frontieres vectorielles si indisponible (pas
                              # de reseau, serveur de tuiles en echec, etc.)
basemap_alpha     <- 0.22     # luminosite/transparence du fond de tuiles
vine_pts_alpha    <- 0.45     # opacite de la couche rasterisee des points WLD
                              # (carte mondiale uniquement -- pas sur les zooms)
vine_pts_res_deg  <- 0.1      # resolution (degres) de cette couche rasterisee

site_label <- tools::file_path_sans_ext(basename(loc_path))
# nom de station = debut du nom de fichier LOC, sans le suffixe de dates
# (ex. "Pessac Haut Lafue Nord_20000101_20260901" -> "Pessac Haut Lafue Nord")
station_name <- sub("_[0-9]{8}(_[0-9]{8})?$", "", site_label)

dir_out <- file.path(data_root, "Analogues_Millesime", paste(station_name, target_year))
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
out_prefix <- file.path(dir_out, paste0("AnaloguesMillesime_", site_label, "_", target_year,
                                        if (cible_source == "pixel") "_pixel" else ""))
dir_sorties <- file.path(dir_out, "sorties")
dir.create(dir_sorties, recursive = TRUE, showWarnings = FALSE)

set.seed(1)

# pdf() can hit a transient "cannot open file" error on Windows right after a
# large pdf() device has just closed in the same folder -- retry.
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

`%||%` <- function(a, b) if (is.null(a)) b else a

vars_map <- c(tmin = "Tmin", tmax = "Tmax", ppt = "P", pet = "ET0")
slot_lab <- sprintf("%02d", months_season)
VARS  <- as.vector(outer(unname(vars_map), slot_lab, paste, sep = "_"))
pcols <- paste0("PC", 1:K)
wcols <- paste0("W", 1:K)

## ---- 1. hemisphere-aware season alignment ----------------------------------
# For a raw calendar year, which (calendar month, vintage-year offset) feeds
# each requested slot. Relation: vintage_year = raw_year - year_delta, hence
# raw_year = vintage_year + year_delta. Northern hemisphere: identity.
# Southern hemisphere: window shifted 6 months, split over two raw years.
# analogues_millesime.md section 1.2. NE PAS MODIFIER (00_note_migration.md
# section 4) : une erreur ici decale silencieusement la moitie du corpus.

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

# distance grand-cercle (haversine), km -- note_station_delta_carte.md B.1
dist_km <- function(lon1, lat1, lon2, lat2) {
  R <- 6371
  p1 <- lat1 * pi/180; p2 <- lat2 * pi/180
  dp <- p2 - p1; dl <- (lon2 - lon1) * pi/180
  a  <- sin(dp/2)^2 + cos(p1) * cos(p2) * sin(dl/2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}

## ---- 1.1 resolve LOC_LON / LOC_LAT before anchoring ------------------------
# note_station_delta_carte.md checklist I: "LOC_LON / LOC_LAT corrects, ordre
# et signe verifies". Auto-detected from loc_path's own LON/LAT columns
# rather than hand-typed, so this can never silently drift from the station's
# actual recorded position.
if (is.na(LOC_LON) || is.na(LOC_LAT)) {
  peek     <- fread(loc_path, nrows = 1)
  lon_peek <- grep("^lon$", names(peek), ignore.case = TRUE, value = TRUE)[1]
  lat_peek <- grep("^lat$", names(peek), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(lon_peek) || is.na(lat_peek))
    stop("LOC_LON/LOC_LAT non renseignes et aucune colonne LON/LAT trouvee dans loc_path -- ",
         "les fixer explicitement en section 0.")
  LOC_LON <- peek[[lon_peek]][1]; LOC_LAT <- peek[[lat_peek]][1]
  message("LOC_LON / LOC_LAT auto-detectes depuis ", basename(loc_path), " : ",
          LOC_LON, " / ", LOC_LAT)
}

## ---- 2. load WLD, check columns, anchor by coordinates ---------------------

message("loading WLD: ", wld_path)
WLD <- as.data.table(read_fst(wld_path))

req <- c(COL_REGION, COL_LON, COL_LAT, "year", "pixID")
if (!all(req %in% names(WLD)))
  stop("colonnes WLD manquantes : ", paste(setdiff(req, names(WLD)), collapse = ", "),
       "\nnames(WLD) = ", paste(head(names(WLD), 40), collapse = ", "))
has_country <- COL_COUNTRY %in% names(WLD)
if (!has_country) {
  # scripts/07_extract_vgdb_point_data.R builds CN_REG as paste(CN, WINE_REGION,
  # sep = "_"): the country code is recoverable as the prefix before the first
  # "_", so we don't have to fall back to NA everywhere.
  COL_COUNTRY <- ".country_derived"
  WLD[, (COL_COUNTRY) := sub("^([^_]+)_.*$", "\\1", get(COL_REGION))]
  has_country <- TRUE
  message("note: colonne pays absente de WLD, derivee du prefixe de \"", COL_REGION, "\"")
}

setnames(WLD, c(COL_LON, COL_LAT), c("lon", "lat"), skip_absent = TRUE)

# -- 2.0 deduplication par pixel TerraClimate --------------------------------
# Le semis VGDB (~1 point/km2) est bien plus dense que la grille TerraClimate
# (~4 km) : plusieurs points VGDB tombent alors dans le meme pixel (meme
# pixID) et portent des donnees climatiques STRICTEMENT IDENTIQUES. Les garder
# tous gonflerait artificiellement l'effectif de ce pixel (poids indu dans la
# covariance intra-site poolee, le KDE, le kNN et tous les comptages de
# points) sans aucune information climatique supplementaire pour le CALCUL.
# Un seul point est retenu par pixID ; si plusieurs CN_REG se partagent un
# meme pixel, on garde celui qui y a le plus gros effectif (nombre de points
# VGDB de ce pixel). Le nombre de points VGDB originaux derriere chaque pixel
# retenu (n_pts_pixel) et le total original par region (region_totals_pts)
# sont conserves : ils ne servent JAMAIS au calcul d'homologie, seulement au
# rapport (section 9), pour dire a combien de points reels un pixel analogue
# correspond.
pts_uniq <- unique(WLD[, .(lon, lat, region = get(COL_REGION), pixID)])
region_totals_pts <- pts_uniq[, .(n_pts_region = .N), by = region]

eff <- pts_uniq[, .N, by = .(pixID, region)]
setorder(eff, pixID, -N)
maj <- eff[!duplicated(pixID), .(pixID, region_maj = region, n_pts_pixel = N)]
pts_keep <- merge(pts_uniq, maj, by = "pixID")[region == region_maj]
pts_keep <- pts_keep[!duplicated(pixID), .(lon, lat, n_pts_pixel)]   # ties within the majority region

n_before <- nrow(pts_uniq)
n_after  <- nrow(pts_keep)
if (n_after < n_before) {
  message("deduplication par pixel TerraClimate (pixID) : ", n_before, " -> ", n_after,
          " points geographiques retenus pour le calcul (", n_before - n_after,
          " doublons de pixel retires -- climat identique bit a bit ; en cas de CN_REG ",
          "multiples sur un meme pixel, la region la plus representee dans ce pixel est ",
          "conservee). Les effectifs de points originaux sont gardes pour le rapport.")
} else {
  message("aucun doublon de pixel TerraClimate detecte dans WLD")
}
WLD <- merge(WLD, pts_keep, by = c("lon", "lat"))

WLD[, site := .GRP, by = .(lon, lat)]
SITE_PTS <- unique(WLD[, .(site, n_pts_pixel)])   # section 9 : n_cohorte_pts, jamais le calcul

# -- 2.1 appariement geographique station -> point WLD le plus proche --------
# note_station_delta_carte.md section B. Rapprochement POINT A POINT (WLD est
# un semis de points a ~1/km2), pas une extraction raster. Aucun centroide
# regional n'intervient nulle part dans ce script.

apparier_station <- function(WLD) {
  PTS <- unique(WLD[, .(site, region = get(COL_REGION),
                        country = if (has_country) get(COL_COUNTRY) else NA_character_,
                        lon, lat)])
  PTS[, d_km := dist_km(LOC_LON, LOC_LAT, lon, lat)]
  setorder(PTS, d_km)

  if (PTS$d_km[1] > D_MAX_KM)
    stop("Point WLD le plus proche a ", round(PTS$d_km[1], 1), " km (> D_MAX_KM = ",
         D_MAX_KM, "). Verifier LOC_LON / LOC_LAT : ordre lon/lat inverse ? ",
         "signe de la longitude a l'ouest de Greenwich ?")

  anc <- as.list(PTS[1])
  message("station appariee au point ", anc$site, " | ", anc$region,
          " (", anc$country, ") | ", round(anc$d_km, 2), " km")
  message("  5 points les plus proches : ",
          paste0(PTS$site[1:5], " (", round(PTS$d_km[1:5], 1), " km)", collapse = " | "))
  anc$PTS <- PTS
  anc
}

ANC  <- apparier_station(WLD)
VOIS <- ANC$PTS[d_km <= VOIS_KM, site]
message("voisinage de controle : ", length(VOIS), " point(s) dans un rayon de ", VOIS_KM, " km")

if (is.na(loc_hemisphere)) {
  loc_hemisphere <- if (LOC_LAT < 0) "S" else "N"
  message("hemisphere LOC deduit de LOC_LAT : ", loc_hemisphere)
}
cal_loc  <- season_calendar(months_season, loc_hemisphere)
mois_cal <- cal_loc$month   # mois calendaires reellement utilises par LOC

## ---- 2.2 serie mensuelle du point apparie (pour le delta, section C) -------
# Extraite AVANT le remodelage saisonnier et avant rm(WLD). Comparaison au
# POINT APPARIE, jamais a une moyenne regionale.

extraire_mensuel <- function(WLD, sites) {
  sub <- WLD[site %in% sites]
  rbindlist(lapply(1:12, function(mo) {
    cols <- paste0(names(vars_map), "_", mo)
    if (!all(cols %in% names(WLD))) return(NULL)
    dt <- data.table(site = sub$site, year = sub$year, mois = mo)
    for (v in names(vars_map)) set(dt, j = vars_map[[v]], value = sub[[paste0(v, "_", mo)]])
    dt
  }))
}

PX_M   <- extraire_mensuel(WLD, ANC$site)          # point apparie
PX_VOIS <- extraire_mensuel(WLD, VOIS)[, lapply(.SD, mean, na.rm = TRUE),
                                       by = .(year, mois), .SDcols = VARS4]  # variante B.2
if (!nrow(PX_M)) stop("point apparie ", ANC$site, " absent de WLD apres extraction mensuelle")

## ---- 2.3 reshape WLD into the season window -------------------------------

build_wld_season <- function(wld, months_season) {
  make_part <- function(sub, hemi) {
    cal <- season_calendar(months_season, hemi)
    parts <- lapply(split(cal, cal$year_delta), function(rows) {
      yd <- rows$year_delta[1]
      p <- sub[, .(site, region = get(COL_REGION), lon, lat, year = year - yd)]
      for (i in seq_len(nrow(rows))) {
        for (v in names(vars_map)) {
          set(p, j = paste0(vars_map[[v]], "_", rows$slot_lab[i]),
              value = sub[[paste0(v, "_", rows$month[i])]])
        }
      }
      p
    })
    Reduce(function(a, b) merge(a, b, by = c("site", "region", "lon", "lat", "year")), parts)
  }
  rbind(make_part(wld[lat >= 0], "N"), make_part(wld[lat < 0], "S"))
}

REG_INFO <- unique(WLD[, .(region = get(COL_REGION),
                           country = if (has_country) get(COL_COUNTRY) else NA_character_)])
REG_INFO <- REG_INFO[!duplicated(region)]

D <- build_wld_season(WLD, months_season)
rm(WLD); gc()

# A handful of coastal pixels can carry NA values (nodata neighbour cells);
# drop any site with an incomplete vintage rather than propagate NAs.
na_sites <- unique(D[!complete.cases(D[, ..VARS]), site])
if (length(na_sites) > 0) {
  message("dropping ", length(na_sites), " WLD site(s) with NA in the season window (coastal nodata)")
  if (ANC$site %in% na_sites) stop("le point d'ancrage est parmi les sites a NA -- ancrage impossible")
  D <- D[!site %in% na_sites]
}

stopifnot(all(VARS %in% names(D)), !anyNA(D[, ..VARS]))
nn <- D[, .N, by = site]$N
n_millesimes_typ <- as.integer(round(median(nn)))   # for report text (section 14); NOT hardcoded
if (length(unique(nn)) != 1L)
  message("note: number of vintages per site is not constant (range ",
          paste(range(nn), collapse = "-"),
          ") -- expected near the reference period's ends for southern-hemisphere sites")
message("WLD sites: ", uniqueN(D$site), " | regions: ", uniqueN(D$region),
        " | rows: ", format(nrow(D), big.mark = " "))
message("NOTE d'inference : la taille d'echantillon effective est de l'ordre de ",
        uniqueN(D$region), " regions, PAS de ", uniqueN(D$site), " points independants.")

## ---- 3. load and aggregate LOC --------------------------------------------

message("loading LOC: ", loc_path)
loc_raw   <- fread(loc_path)
loc_dates <- as.IDate(loc_raw[[loc_cols$date]])
loc_raw[, `:=`(year = year(loc_dates), mois = month(loc_dates))]

if (loc_freq == "daily") {
  LOC_m <- loc_raw[, .(Tmin = mean(get(loc_cols$tmin), na.rm = TRUE),
                       Tmax = mean(get(loc_cols$tmax), na.rm = TRUE),
                       P    = sum(get(loc_cols$ppt),  na.rm = TRUE),
                       ET0  = sum(get(loc_cols$pet),  na.rm = TRUE),
                       n_j  = .N), by = .(year, mois)]
  jours_att <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  LOC_m[, compl := n_j / jours_att[mois]]
  if (LOC_m[compl < COMPL_MIN, .N] > 0) {
    warning("mois LOC incomplets (<", 100 * COMPL_MIN, " % de jours) -- ecartes :")
    print(LOC_m[compl < COMPL_MIN, .(year, mois, n_j, compl = round(compl, 2))])
  }
  LOC_m <- LOC_m[compl >= COMPL_MIN]
} else if (loc_freq == "monthly") {
  LOC_m <- loc_raw[, .(year, mois,
                       Tmin = get(loc_cols$tmin), Tmax = get(loc_cols$tmax),
                       P    = get(loc_cols$ppt),  ET0  = get(loc_cols$pet))]
} else stop("loc_freq must be \"daily\" or \"monthly\"")

## ---- 4. delta de source : calibration et qualification --------------------
# note_station_delta_carte.md sections C et D. Un pixel de ~4 km et une
# station ne mesurent pas la meme chose ; ce biais entre dans toutes les
# distances, y compris celle du site a lui-meme. Non corrige, il deplace le
# classement entier de facon systematique et invisible.

calibrer_delta <- function(LOC_m, px_m, mois_sel = mois_cal) {
  comm <- intersect(LOC_m$year, px_m$year)
  if (length(comm) < N_AN_MIN)
    stop("seulement ", length(comm), " annee(s) commune(s) LOC / point apparie : ",
         "calibrage impossible (minimum ", N_AN_MIN, ")")
  message("calibrage du delta sur ", length(comm), " annees communes : ",
          paste(range(comm), collapse = "-"))

  M <- merge(LOC_m[year %in% comm & mois %in% mois_sel],
             px_m [year %in% comm & mois %in% mois_sel],
             by = c("year", "mois"), suffixes = c("_loc", "_px"))
  for (v in ADD) M[[paste0("d_", v)]] <- M[[paste0(v, "_loc")]] - M[[paste0(v, "_px")]]
  for (v in MUL) M[[paste0("d_", v)]] <- M[[paste0(v, "_loc")]] /
                                          pmax(M[[paste0(v, "_px")]], P_MIN)
  M[]
}

qualifier_delta <- function(DLT) {
  out <- rbindlist(lapply(VARS4, function(v) {
    col <- paste0("d_", v)
    DLT[, {
      ok <- is.finite(get(col)); x <- get(col)[ok]; yy <- year[ok]; n <- length(x)
      mu <- mean(x); md <- median(x); sdv <- sd(x)
      tt <- if (n >= 6L) {
              ct <- suppressWarnings(cor.test(yy, x, method = "spearman"))
              c(rho = unname(ct$estimate), p = ct$p.value)
            } else c(rho = NA_real_, p = NA_real_)
      .(variable = v, n = n,
        delta = round(mu, 3), delta_med = round(md, 3), sd = round(sdv, 3),
        snr = round(abs(mu) / sdv, 2),
        rho_an = round(unname(tt["rho"]), 2), p_tend = signif(unname(tt["p"]), 2))
    }, by = mois]
  }))
  out[, fiabilite := fifelse(!is.na(p_tend) & p_tend < 0.05 & abs(rho_an) > 0.5,
                             "DERIVE - ne pas extrapoler",
                     fifelse(snr >= 2, "fiable",
                     fifelse(snr >= 1, "reserve", "ECARTER - non systematique")))]
  setorder(out, variable, mois)
  out[]
}

DLT <- calibrer_delta(LOC_m, PX_M)
QD  <- qualifier_delta(DLT)
print(QD)
message("--- synthese fiabilite du delta ---")
print(QD[, .N, by = fiabilite])
if (QD[grepl("ECARTER|DERIVE", fiabilite), .N] > 0)
  warning("Certains couples variable x mois ne sont pas corrigibles de facon fiable ",
          "(colonne 'fiabilite'). Les resultats sur ces mois doivent etre presentes ",
          "avec reserve explicite.")

fwrite(QD, paste0(out_prefix, "_delta_qualification.csv"))

# variante de robustesse : delta calibre sur le voisinage 10 km (B.2)
QD_VOIS <- tryCatch(qualifier_delta(calibrer_delta(LOC_m, PX_VOIS)), error = function(e) NULL)
if (!is.null(QD_VOIS))
  message("controle voisinage ", VOIS_KM, " km : ecart median des deltas au point apparie = ",
          round(median(abs(QD$delta - QD_VOIS$delta)), 3))

## ---- 4.1 vecteur cible : yb_brut et yb_corr -------------------------------
# Sens de la correction : on ramene la station VERS l'echelle du point WLD,
# jamais l'inverse. La base de reference n'est jamais touchee.

vecteur_cible <- function(LOC_m, QD, hemisphere, an, corriger = TRUE, mode = delta_mode) {
  cal <- season_calendar(months_season, hemisphere)
  # raw_year = vintage_year + year_delta  (cf. section 1)
  cal[, cal_year := an + year_delta]
  m <- merge(cal, LOC_m, by.x = c("cal_year", "month"), by.y = c("year", "mois"), all.x = TRUE)
  setorder(m, slot)
  if (nrow(m) != nrow(cal) || anyNA(m[, ..VARS4]))
    stop("millesime ", an, " incomplet dans LOC : mois manquants ou ecartes pour completude. ",
         "Verifier la couverture de loc_path.")
  out <- numeric(0)
  for (i in seq_len(nrow(m))) {
    for (v in VARS4) {
      val <- m[[v]][i]
      if (corriger) {
        d <- QD[variable == v & mois == m$month[i],
                if (mode == "median") delta_med else delta]
        if (length(d) != 1L || !is.finite(d))
          stop("delta indisponible pour ", v, " mois ", m$month[i])
        val <- if (v %in% ADD) val - d else val / d
      }
      out[paste0(v, "_", m$slot_lab[i])] <- val
    }
  }
  out[VARS]
}

# vecteur cible tire du point WLD apparie : aucun delta a appliquer, la cible
# est deja a l'echelle de la base de reference (run A du protocole E.1).
vecteur_cible_pixel <- function(px_m, hemisphere, an) {
  cal <- season_calendar(months_season, hemisphere)
  cal[, cal_year := an + year_delta]
  m <- merge(cal, px_m, by.x = c("cal_year", "month"), by.y = c("year", "mois"), all.x = TRUE)
  setorder(m, slot)
  if (nrow(m) != nrow(cal) || anyNA(m[, ..VARS4]))
    stop("millesime ", an, " absent ou incomplet dans WLD pour le point apparie ", ANC$site,
         " -- le run \"pixel\" n'est possible que sur une annee couverte par WLD")
  out <- numeric(0)
  for (i in seq_len(nrow(m)))
    for (v in VARS4) out[paste0(v, "_", m$slot_lab[i])] <- m[[v]][i]
  out[VARS]
}

if (cible_source == "station") {
  yb_corr <- vecteur_cible(LOC_m, QD, loc_hemisphere, target_year, corriger = TRUE)
  yb_brut <- vecteur_cible(LOC_m, QD, loc_hemisphere, target_year, corriger = FALSE)
} else if (cible_source == "pixel") {
  yb_corr <- vecteur_cible_pixel(PX_M, loc_hemisphere, target_year)
  yb_brut <- yb_corr
  run_brut_corr <- FALSE
  message("cible tiree du point WLD apparie : aucun delta applique")
} else stop("cible_source doit valoir \"station\" ou \"pixel\"")

message("profil cible (", site_label, ", ", target_year, ", source ", cible_source, ") :")
print(round(yb_corr, 2))
if (cible_source == "station") {
  message("ecart brut - corrige :")
  print(round(yb_brut - yb_corr, 2))
}

yb <- yb_corr   # vecteur de reference du pipeline

## ---- 5. detrending, rescaling to target_year, normal scores ---------------
# analogues_millesime.md sections 2.1-2.2. Retrait de tendance AVANT la
# transformation en scores normaux. Une copie en unites physiques est
# conservee pour les graphiques interpretables.

TARGET_YEAR <- target_year

D[, (VARS) := {
  yy  <- year
  ww  <- yy - mean(yy)
  Sww <- sum(ww^2)
  lapply(.SD, function(x) x - (sum(ww * x) / Sww) * (yy - TARGET_YEAR))
}, by = site, .SDcols = VARS]

D_phys       <- copy(D)     # detrended, physical units -- plots only
yb_phys      <- yb_corr
yb_phys_brut <- yb_brut

# Scores normaux (copule gaussienne) sur TOUTES les variables. Remplace toute
# transformation ad hoc de type log1p sur les precipitations.
# ORDRE IMPORTANT : transformer les cibles AVANT d'ecraser la colonne de D.
EPS <- 1e-6
for (v in VARS) {
  Fv         <- ecdf(D[[v]])
  yb_corr[v] <- qnorm(pmin(pmax(Fv(yb_corr[v]), EPS), 1 - EPS))
  yb_brut[v] <- qnorm(pmin(pmax(Fv(yb_brut[v]), EPS), 1 - EPS))
  set(D, j = v, value = qnorm(pmin(pmax(Fv(D[[v]]), EPS), 1 - EPS)))
}
yb <- yb_corr

## ---- 6. PCA space, diagnostics, pooled intra-site covariance metric -------
# ---------------------------------------------------------------------------
#  METRIQUE DE FOND : distance de Mahalanobis sur covariance intra-site poolee
#    Mahalanobis (1936), Proc. Natl. Inst. Sci. India 2(1), 49-55
#                       reed. Sankhya A doi:10.1007/s13171-019-00164-5
#  APPLICATION AUX ANALOGUES CLIMATIQUES
#    Mahony et al. (2017)  doi:10.1111/gcb.13645  [sigma dissimilarity]
#    Grenier et al. (2013) doi:10.1175/JAMC-D-12-0170.1  [comparaison de six
#                          metriques de dissimilarite climatique]
#  COVARIANCE INTRA-GROUPE POOLEE COMME METRIQUE
#    Mardia, Kent & Bibby (1979), Multivariate Analysis [metrique de l'ADL]
# ---------------------------------------------------------------------------

X   <- as.matrix(D[, ..VARS])
pca <- prcomp(X, center = TRUE, scale. = TRUE)

# -- 6.1 bloc diagnostique ACP (analogues_millesime.md section 4) -----------
# AVERTISSEMENT : une seule ACP pour tout le jeu mondial. Ces pourcentages
# sont GLOBAUX, ils ne decrivent aucun site en particulier.
vp  <- pca$sdev^2
acp <- data.table(CP      = paste0("PC", seq_along(vp)),
                  var_pct = round(100 * vp / sum(vp), 2),
                  cum_pct = round(100 * cumsum(vp) / sum(vp), 2))
print(head(acp, 12))
message("K = ", K, " -> variance cumulee retenue : ", acp$cum_pct[K], " %")

charg <- round(pca$rotation[, 1:K] * rep(pca$sdev[1:K], each = nrow(pca$rotation)), 2)
message("--- correlations variables d'origine / composantes ---")
print(charg)
for (j in 1:K) {
  o <- order(abs(charg[, j]), decreasing = TRUE)[1:5]
  message("PC", j, " (", acp$var_pct[j], " %) : ",
          paste0(rownames(charg)[o], " ", charg[o, j], collapse = " | "))
}
fwrite(acp, paste0(out_prefix, "_acp_variance.csv"))
fwrite(data.table(variable = rownames(charg), as.data.table(charg)),
       paste0(out_prefix, "_acp_correlations.csv"))

Z       <- pca$x[, 1:K, drop = FALSE]
zb_corr <- drop(predict(pca, matrix(yb_corr, nrow = 1, dimnames = list(NULL, VARS))))[1:K]
zb_brut <- drop(predict(pca, matrix(yb_brut, nrow = 1, dimnames = list(NULL, VARS))))[1:K]
zb      <- zb_corr
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

Rw  <- chol(W)                                   # W = Rw'Rw
Zw  <- t(backsolve(Rw, t(Z), transpose = TRUE))  # euclidien dans Zw == Mahalanobis dans Z
DTw <- data.table(site = D$site, region = D$region, year = D$year)
DTw[, (wcols) := as.data.table(Zw)]

## ---- 7. les trois estimateurs ---------------------------------------------

# ===========================================================================
# 7.1  ESTIMATEUR PRINCIPAL -- predictif bayesien a a priori empirique
#      (Normal-inverse-Wishart, predictif Student). Forme fermee : une
#      decomposition de Cholesky K x K par site, aucun MCMC.
#
#  SOURCE DE LA METHODE
#    Raiffa & Schlaifer (1961)  doi:10.1002/9781118625125   [cadre conjugue]
#    Gelman et al. (2013), Bayesian Data Analysis 3e, ch. 3
#                               doi:10.1201/b16018          [predictif Student, MAJ NIW]
#    Schafer & Strimmer (2005)  doi:10.2202/1544-6115.1175  [shrinkage de covariance]
#    Robbins (1956) ; Efron (2010) doi:10.1017/CBO9780511761362  [empirical Bayes]
#  EXEMPLE D'APPLICATION EN ENVIRONNEMENT
#    Mahony et al. (2017), Global Change Biology
#                               doi:10.1111/gcb.13645       [Mahalanobis regularise,
#                                                            analogues climatiques]
#
#  STATUT : a priori estime sur les memes donnees => EMPIRICAL BAYES, pas
#  bayesien strict. Les incertitudes a posteriori sont sous-estimees. Le
#  score est utilise comme STATISTIQUE DE CLASSEMENT ORDINAL, jamais comme
#  une probabilite (analogues_millesime.md 5.1.0bis).
# ===========================================================================

nu0  <- K + 2 + M0
Psi0 <- (nu0 - K - 1) * W          # => E[Sigma] = W a priori

# log-densite d'une Student multivariee, parametree par matrice d'echelle
ldstudent <- function(x, mu, Scale, df) {
  R <- chol(Scale)
  z <- backsolve(R, x - mu, transpose = TRUE)
  lgamma((df + K) / 2) - lgamma(df / 2) -
    (K / 2) * log(df * pi) - sum(log(diag(R))) -
    ((df + K) / 2) * log1p(sum(z^2) / df)
}

scorer_bayesien <- function(zb, m0 = M0) {
  nu0_  <- K + 2 + m0
  Psi0_ <- (nu0_ - K - 1) * W
  DT[, {
    M    <- as.matrix(.SD)
    n    <- nrow(M)
    xbar <- colMeans(M)

    kn   <- KAPPA0 + n
    nu_n <- nu0_   + n
    mu0 <- xbar                                # a priori non informatif sur la position ;
    mn  <- (KAPPA0 * mu0 + n * xbar) / kn      # forme generale conservee pour permettre
    d0  <- xbar - mu0                          # de tester un mu0 regional (=> d0 non nul)
    Sn  <- if (n >= 2L) (n - 1) * cov(M) else matrix(0, K, K)
    Pn  <- Psi0_ + Sn + (KAPPA0 * n / kn) * tcrossprod(d0)

    df    <- nu_n - K + 1
    Scale <- Pn * (kn + 1) / (kn * df)
    .(score_bay = ldstudent(zb, mn, Scale, df), df_post = df)
  }, by = .(site, region), .SDcols = pcols]
}

# ---------------------------------------------------------------------------
#  7.2 CONTROLE 1 -- noyau gaussien (KDE)
#    Rosenblatt (1956)  doi:10.1214/aoms/1177728190
#    Parzen (1962)      doi:10.1214/aoms/1177704472
#  EXEMPLES D'APPLICATION EN ENVIRONNEMENT
#    Broennimann et al. (2012), Global Ecol. Biogeogr.
#                       doi:10.1111/j.1466-8238.2011.00698.x
#                       [noyaux en espace climatique issu d'une ACP -- cas le
#                        plus proche du notre ; package ecospat]
#    Qiao et al. (2017), Global Ecol. Biogeogr.  doi:10.1111/geb.12492
#                       [GARDE-FOU : biais du KDE selon dimension et effectif]
#
#  7.3 CONTROLE 2 -- densite kNN adaptative
#    Loftsgaarden & Quesenberry (1965)  doi:10.1214/aoms/1177700079
#  EXEMPLE D'APPLICATION EN ENVIRONNEMENT
#    Lall & Sharma (1996), Water Resources Research  doi:10.1029/95WR02966
#                       [kNN sur series climatiques multivariees]
# ---------------------------------------------------------------------------

h <- sqrt(qchisq(0.50, df = K)) / 2
r <- sqrt(qchisq(q_count, df = K))

# mediane de 3 sans apply() : 47 000 lignes x 3 colonnes
med3  <- function(a, b, c) pmax(pmin(a, b), pmin(pmax(a, b), c))
span3 <- function(a, b, c) pmax(a, b, c) - pmin(a, b, c)

scorer_cible <- function(zb, etiquette = "") {
  zbw <- backsolve(Rw, zb, transpose = TRUE)
  d   <- sqrt(rowSums(sweep(Zw, 2, zbw, "-")^2))
  DTd <- data.table(site = DT$site, region = DT$region, year = DT$year, d = d)

  sc_bay <- scorer_bayesien(zb)
  sc_kde <- DTd[, .(score_kde = mean(exp(-0.5 * (d / h)^2))), by = .(site, region)]
  sc_knn <- DTd[, { dm <- sort(d)[m_knn]
                    .(score_knn = log(m_knn) - log(.N) - K * log(dm), d_m = dm) },
                by = .(site, region)]
  # statistique de COMMUNICATION, jamais de classement
  sc_cnt <- DTd[, .(n_proches = sum(d <= r), n_tot = .N,
                    d_min = min(d), an_proche = year[which.min(d)]), by = .(site, region)]

  res <- Reduce(function(a, b) merge(a, b, by = c("site", "region")),
                list(sc_bay, sc_kde, sc_knn, sc_cnt))

  res[, `:=`(rg_bay = frank(-score_bay, ties.method = "min"),
             rg_kde = frank(-score_kde, ties.method = "min"),
             rg_knn = frank(-score_knn, ties.method = "min"))]

  # consensus par RANG MEDIAN des trois methodes. Ce n'est PAS du model
  # averaging : aucun poids de vraisemblance de modele n'intervient. C'est une
  # agregation ordinale robuste (analogues_millesime.md 5.1.0bis, point 3).
  res[, `:=`(rg_med  = med3(rg_bay, rg_kde, rg_knn),
             rg_span = span3(rg_bay, rg_kde, rg_knn))]
  res[, rg_cons := frank(rg_med, ties.method = "min")]

  N1 <- ceiling(0.01 * nrow(res))
  res[, accord := (rg_bay <= N1) + (rg_kde <= N1) + (rg_knn <= N1)]
  res[, solidite := fifelse(accord == 3L, "3/3",
                    fifelse(accord == 2L, "2/3",
                    fifelse(accord == 1L, "1/3", "0/3")))]
  setorder(res, rg_cons)
  if (nzchar(etiquette)) message("scores calcules (", etiquette, ") : ", nrow(res), " sites")
  res[]
}

res <- scorer_cible(zb_corr, "cible corrigee")
res <- merge(res, SITE_PTS, by = "site")   # n_pts_pixel : reporting seulement, cf. section 2.0
message("degres de liberte a posteriori : ", res$df_post[1],
        "  (plus eleve = plus proche de la gaussienne)")

## ---- 7.4 tableau principal : 4 colonnes de rang ---------------------------

TAB <- res[, .(site, region,
               bayesien = rg_bay, kde = rg_kde, knn = rg_knn,
               consensus = rg_cons, span = rg_span, solidite,
               n_proches, n_tot, d_min, an_proche)]
message("\ntop 30 points (rang consensus) -- RAPPEL : le top-1 parcellaire n'a aucune signification")
print(head(TAB, 30))

top_par_methode <- function(dat, N = 10L, unite = "region") {
  f <- function(col) {
    d <- dat[order(dat[[col]])]
    if (unite == "region") d <- d[!duplicated(region)]
    d[1:N, get(unite)]
  }
  data.table(rang = 1:N, bayesien = f("rg_bay"), kde = f("rg_kde"),
             knn = f("rg_knn"), consensus = f("rg_cons"))[]
}
cote_a_cote <- top_par_methode(res, N = N_REG, unite = "region")
print(cote_a_cote)
freq <- sort(table(unlist(cote_a_cote[, -1])), decreasing = TRUE)
message("--- presence dans les top-", N_REG, " (max 4) ---")
print(freq)

## ---- 7.5 effet de la correction de source : yb_brut vs yb_corr ------------
# note_station_delta_carte.md C (conserver yb_brut) et E.2 (metriques).

if (run_brut_corr) {
  res_brut <- scorer_cible(zb_brut, "cible brute")
  M <- merge(res[, .(site, region, rgC = rg_cons)],
             res_brut[, .(site, rgB = rg_cons)], by = "site")
  topC <- res[!duplicated(region)][1:50, region]
  topB <- res_brut[!duplicated(region)][1:50, region]
  EFF <- list(recouvrement_top50 = length(intersect(topC, topB)) / 50,
              spearman = cor(M$rgC, M$rgB, method = "spearman"),
              only_corr = setdiff(topC, topB), only_brut = setdiff(topB, topC))
  message("\neffet de la correction de source (corrige vs brut) : recouvrement top-50 = ",
          round(100 * EFF$recouvrement_top50), " % | Spearman = ", round(EFF$spearman, 3))
  if (EFF$recouvrement_top50 < 0.5)
    warning("recouvrement < 50 % : le biais de source domine le classement. ",
            "Revoir le delta (section D) avant toute interpretation.")
}

## ---- 8. profil mensuel : le mois qui disqualifie ---------------------------
# analogues_millesime.md section 6.

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
prof      <- data.table(site = D$site, region = D$region, Dm)
prof_site <- prof[, lapply(.SD, median), by = .(site, region), .SDcols = slot_lab]
prof_site[, mois_pire := month.abb[months_season][apply(.SD, 1, which.max)], .SDcols = slot_lab]

res <- merge(res, prof_site, by = c("site", "region"))
setorder(res, rg_cons)
fwrite(res, paste0(out_prefix, "_table.csv"))
message("wrote ", paste0(out_prefix, "_table.csv"))

## ---- 9. cohortes regionales sequentielles ---------------------------------
# note_station_delta_carte.md G.2. Aucune moyenne, aucun centroide : la
# diversite climatique interne d'une region est une propriete reelle de cette
# region, pas du bruit. Un analogue est UN LIEU.

construire_cohortes <- function(res, N_REG = 10L) {
  # n_cohorte_pix compte les PIXELS retenus (l'unite du calcul, dedupliquee en
  # 2.0) ; n_cohorte_pts compte les points VGDB originaux qu'ils representent
  # (region_totals_pts / n_pts_pixel, cf. section 2.0) -- utilises ICI
  # seulement pour le rapport, jamais pour reponderer le classement.
  setorder(res, rg_cons)
  R <- copy(res)[, rang_pt := .I]

  R[, bloc := rleid(region)]
  COH <- R[, .(region        = region[1],
               rang_deb      = min(rang_pt),
               rang_fin      = max(rang_pt),
               n_cohorte_pix = .N,
               n_cohorte_pts = sum(n_pts_pixel),
               rg_med        = median(as.double(rg_cons)),
               solidite      = solidite[1]), by = bloc]
  # une region peut reapparaitre plus loin : on garde son premier bloc
  COH <- COH[!duplicated(region)]
  setorder(COH, rang_deb)
  COH[, rang_region := .I]
  COH <- merge(COH, region_totals_pts, by = "region", all.x = TRUE)
  setnames(COH, "n_pts_region", "n_region_pts")
  COH[, pct_region := round(100 * n_cohorte_pts / n_region_pts, 1)]
  setorder(COH, rang_region)
  head(COH, N_REG)[]
}

mesurer_dispersion <- function(COH, res, DTw) {
  rgs <- res[region %in% COH$region,
             .(rg_med_reg = median(as.double(rg_cons)),
               rg_q25     = quantile(rg_cons, .25),
               rg_q75     = quantile(rg_cons, .75),
               n_top1pct  = sum(rg_cons <= ceiling(0.01 * nrow(res)))), by = region]

  # etalement climatique : dispersion des centres de points, espace blanchi
  ctr <- DTw[region %in% COH$region, lapply(.SD, mean), by = .(site, region), .SDcols = wcols]
  spr <- ctr[, { M <- as.matrix(.SD); mu <- colMeans(M)
                 dd <- sqrt(rowSums(sweep(M, 2, mu, "-")^2))
                 .(etalement = round(median(dd), 2), etal_p90 = round(quantile(dd, .9), 2)) },
             by = region, .SDcols = wcols]

  out <- merge(merge(COH, rgs, by = "region"), spr, by = "region")
  out[, profil := fifelse(pct_region >= 50, "region homogene et proche",
                  fifelse(pct_region >= 15, "part substantielle de la region",
                  fifelse(n_cohorte_pix >= 5, "sous-secteur localise",
                                              "point isole - a signaler")))]
  setorder(out, rang_region)
  out[]
}

COH <- construire_cohortes(res, N_REG)
SYN <- mesurer_dispersion(COH, res, DTw)
SYN <- merge(SYN, REG_INFO, by = "region", all.x = TRUE)
setorder(SYN, rang_region)

n_proches_reg <- res[region %in% SYN$region & rg_cons <= max(SYN$rang_fin),
                     .(n_proches = round(median(as.double(n_proches)), 1)), by = region]
SYN <- merge(SYN, n_proches_reg, by = "region", all.x = TRUE)
setorder(SYN, rang_region)

TAB_SYN <- SYN[, .(rang = rang_region, region, country,
                   rang_deb, rang_fin, n_cohorte_pts, n_region_pts, pct_region,
                   n_cohorte_pix, rg_med_reg, rg_q25, rg_q75, etalement, etal_p90,
                   n_proches, solidite, profil)]
print(TAB_SYN)

fwrite(TAB_SYN, file.path(dir_sorties, "synthese_top10_regions.csv"))
fwrite(res[rg_cons <= max(SYN$rang_fin)], file.path(dir_sorties, "points_retenus.csv"))

## ---- 10. diagnostics obligatoires -----------------------------------------

pdf(paste0(out_prefix, "_diagnostics.pdf"), width = 9, height = 6)

# -- 10.1 etalon interne : la region d'ancrage contre elle-meme -------------
# analogues_millesime.md 8.1 ; note_station_delta_carte.md F.3.
# n_proches ne veut rien dire dans l'absolu : il se lit relativement au
# n_proches de la region d'ancrage.
anc_pt  <- res[site == ANC$site]
anc_reg <- res[region == ANC$region]
message("\n--- etalon interne : point et region d'ancrage ---")
if (nrow(anc_pt)) {
  message("point apparie ", ANC$site, " : rang consensus ", anc_pt$rg_cons,
          " / ", nrow(res), " | n_proches = ", anc_pt$n_proches, " / ", anc_pt$n_tot)
  n_proches_anc <- anc_pt$n_proches
} else { n_proches_anc <- NA_integer_ }
if (nrow(anc_reg)) {
  message("region d'ancrage (", ANC$region, ", ", nrow(anc_reg), " points) : n_proches")
  print(summary(anc_reg$n_proches))
  print(anc_reg[order(rg_cons)][1:min(20, .N),
        .(site, rg_bay, rg_kde, rg_knn, rg_cons, n_proches, n_tot)])
}
# NUANCE (note_station_delta_carte.md E.3) : un millesime "ovni" DOIT mal
# classer sa propre region. Le controle d'ancrage n'est interpretable que sur
# un millesime ordinaire.

# -- 10.2 score d'ancrage en leave-one-out (H.1, adapte au predictif NIW) ---
# Si le millesime cible est present dans WLD, il participe au nuage de sa
# propre region et l'avantage mecaniquement. Ici la cible vient de la station
# et l'annee cible n'est en general pas dans WLD : le controle est alors sans
# objet, et le message le dit.
score_ancrage_loo <- function(DT, region_anc, an, zb) {
  M <- as.matrix(DT[region == region_anc & year != an, ..pcols])
  n <- nrow(M); xbar <- colMeans(M)
  kn <- KAPPA0 + n; nu_n <- nu0 + n
  Pn <- Psi0 + (n - 1) * cov(M)
  df <- nu_n - K + 1
  ldstudent(zb, xbar, Pn * (kn + 1) / (kn * df), df)
}
if (target_year %in% unique(DT$year)) {
  message("score d'ancrage (region ", ANC$region, ") en leave-one-out de ", target_year,
          " : ", round(score_ancrage_loo(DT, ANC$region, target_year, zb), 3))
} else {
  message("millesime ", target_year, " absent de WLD : aucune exclusion leave-one-out necessaire")
}

# -- 10.3 bootstrap sur millesimes -- diagnostic central (8.2) --------------
if (run_bootstrap) {
  boot_top <- function(B = boot_B, N = boot_N) {
    ref <- res[order(rg_bay)][1:N, site]
    ov  <- numeric(B)
    for (b in seq_len(B)) {
      idx <- DT[, .I[sample(.N, .N, replace = TRUE)], by = site]$V1
      sb  <- DT[idx][, {
        M  <- as.matrix(.SD); n <- nrow(M)
        kn <- KAPPA0 + n; nu_n <- nu0 + n
        Pn <- Psi0 + (n - 1) * cov(M)
        df <- nu_n - K + 1
        .(sc = ldstudent(zb, colMeans(M), Pn * (kn + 1) / (kn * df), df))
      }, by = site, .SDcols = pcols]
      ov[b] <- length(intersect(ref, sb[order(-sc)][1:N, site])) / N
    }
    ov
  }
  ov <- boot_top()
  message("\nrecouvrement bootstrap du top-", boot_N, " (B=", boot_B, ") : ",
          round(100 * mean(ov)), " % [", round(100 * min(ov)), "-", round(100 * max(ov)), "]")
  if (mean(ov) < 0.5)
    warning("recouvrement bootstrap < 50 % : NE RIEN PUBLIER au niveau parcellaire. ",
            "Passer a la phase 2 regionale (analogues_millesime.md section 9).")
}

# -- 10.4 sensibilite a M0 (5.1.3) -- M0 est DECLARE, pas optimise ----------
if (run_sens_M0) {
  ref200 <- res[order(rg_bay)][1:200, site]
  SENS <- rbindlist(lapply(M0_grid, function(m) {
    s <- scorer_bayesien(zb, m0 = m)
    tp <- s[order(-score_bay)][1:200, site]
    data.table(M0_val = m, recouvrement_top200 = length(intersect(ref200, tp)) / 200)
  }))
  print(SENS)
  if (SENS[M0_val %in% c(5L, 30L), min(recouvrement_top200)] < 0.85)
    warning("recouvrement du top-200 < 0,85 entre M0 = 5 et M0 = 30 : les nuages locaux ",
            "sont trop peu informatifs et SEUL l'a priori parle -- a signaler explicitement.")
  fwrite(SENS, paste0(out_prefix, "_sensibilite_M0.csv"))
}

# -- 10.5 concordance entre estimateurs (8.3) -------------------------------
topN <- function(col, N = 200L) res[order(res[[col]])][1:N, site]
message("\nconcordance des top-200 entre estimateurs :")
for (p in list(c("rg_bay", "rg_knn"), c("rg_bay", "rg_kde"), c("rg_kde", "rg_knn")))
  message(p[1], " / ", p[2], " : ", length(intersect(topN(p[1]), topN(p[2]))), "/200")

# -- 10.6 coherence spatiale intra-region -- le juge de paix (8.5) ----------
# La specification 8.5 raisonne sur rg_bay, la section G.8 sur rg_cons : les
# deux sont rapportes, le graphique utilise le consensus.
coh <- res[, .(n = .N,
               rg_med_bay = median(as.double(rg_bay)), rg_iqr_bay = IQR(rg_bay),
               rg_med = median(as.double(rg_cons)), rg_iqr = IQR(rg_cons),
               part_top5pct = mean(rg_cons <= 0.05 * nrow(res))), by = region][order(rg_med)]
print(head(coh, 25))
# log = "x" seulement : les regions a un seul point ont rg_iqr = 0.
plot(coh$rg_med, coh$rg_iqr, log = "x", pch = 16, cex = .5,
     xlab = "rang median", ylab = "IQR des rangs",
     main = "Coherence spatiale intra-region (bas-gauche = credible)")

# -- 10.7 curseur du seuil de comptage (F.1) --------------------------------
zbw_main <- backsolve(Rw, zb, transpose = TRUE)
d_main   <- sqrt(rowSums(sweep(Zw, 2, zbw_main, "-")^2))
DTd_main <- data.table(site = DT$site, region = DT$region, d = d_main)
for (q in c(0.25, 0.35, 0.50, 0.65)) {
  rq  <- sqrt(qchisq(q, df = K))
  cnt <- DTd_main[, .(n = sum(d <= rq)), by = .(site, region)]
  message("q = ", q, " (r = ", round(rq, 2), ") | n_proches median = ", median(cnt$n),
          " | ancrage = ", cnt[region == ANC$region, round(median(n), 1)])
}

# -- 10.8 calibration gaussienne (8.6) --------------------------------------
qqplot(qchisq(ppoints(1e5), df = K), sample(rowSums(Zw^2), min(1e5, nrow(Zw))),
       xlab = "chi2 theorique", ylab = "distances au carre observees",
       main = "Calibration gaussienne")
abline(0, 1, col = "red")

dev.off()
message("wrote ", paste0(out_prefix, "_diagnostics.pdf"))

## ---- 11. palette commune a toutes les sorties (G.5) -----------------------

PAL <- c(setNames(hcl.colors(nrow(SYN), "Dark 3"), SYN$region), ancrage = "#111111")
SYN[, couleur := PAL[region]]
# base-R plotting characters (pch 0-9: square, circle, triangle, plus, cross,
# diamond, ... -- one shape per region, reused by graph_profils())
PCH <- setNames(0:(nrow(SYN) - 1), SYN$region)

## ---- 12. carte -- points reels uniquement, aucun centroide (G.6) ----------
# Fond de carte : tuiles OpenTopoMap (attenuees) si disponibles, sinon repli
# sur les frontieres vectorielles. Les ~46 000 points WLD (contexte mondial,
# sans poids analytique propre) sont RASTERISES en une seule fois pour
# alleger le PDF ; les points RETENUS et de CONTEXTE REGIONAL (le resultat)
# restent des points vectoriels reels, jamais de centroide.

world_bounds <- st_read(file.path(sig_root, "world/World_AdmBoundaries_Countries/world-administrative-boundaries.shp"),
                        quiet = TRUE)

site_coords <- unique(D[, .(site, region, lon, lat)])
PT_TOP <- merge(res[rg_cons <= max(SYN$rang_fin), .(site, region, rg_cons)],
                site_coords[, .(site, lon, lat)], by = "site")
PT_TOP[, couleur := PAL[region]]
fwrite(PT_TOP, file.path(dir_sorties, "points_carte.csv"))

loc_coord <- data.table(lon = LOC_LON, lat = LOC_LAT)
anc_coord <- data.table(lon = ANC$lon, lat = ANC$lat)

# -- 12.1 couche rasterisee des ~46 000 points WLD (une seule fois, reutilisee
#    par toutes les pages : coord_sf() ne fait que cadrer la vue) -----------
wld_pts_v  <- vect(site_coords, geom = c("lon", "lat"), crs = "EPSG:4326")
wld_grid_r <- rast(ext(-180, 180, -60, 85), resolution = vine_pts_res_deg, crs = "EPSG:4326")
wld_pres_r <- terra::rasterize(wld_pts_v, wld_grid_r, field = 1, background = NA)

# -- 12.2 fond de tuiles OpenTopoMap, avec repli et nouvelles tentatives ----
# Un echec au tout premier appel (latence/DNS) est courant et transitoire.
choose_tile_zoom <- function(width_deg, target_tiles = 6, zmin = 2, zmax = 11)
  max(zmin, min(zmax, round(log2(target_tiles * 360 / max(width_deg, 0.01)))))

safe_get_tiles <- function(xlim, ylim, provider = "OpenTopoMap", max_tries = 3, wait = 3) {
  zoom <- choose_tile_zoom(diff(xlim))
  bbox_v <- vect(ext(xlim[1], xlim[2], ylim[1], ylim[2]), crs = "EPSG:4326")
  for (i in seq_len(max_tries)) {
    r <- tryCatch(get_tiles(bbox_v, provider = provider, zoom = zoom, crop = TRUE),
                  error = function(e) { message("basemap tile fetch attempt ", i, " failed: ",
                                                conditionMessage(e)); NULL })
    if (!is.null(r)) return(r)
    Sys.sleep(wait * i)
  }
  message("basemap tiles unavailable after ", max_tries, " tries -- falling back to vector boundaries")
  NULL
}

build_carte <- function(regions_n, bbox_pts = NULL, titre = NULL) {
  ctx <- site_coords[region %in% regions_n]           # contexte : tous les points des memes regions
  sel <- PT_TOP[region %in% regions_n]                # cohortes retenues
  ctx[, region := factor(region, levels = regions_n)]
  sel[, region := factor(region, levels = regions_n)]

  if (is.null(bbox_pts)) {
    xlim <- c(-180, 180); ylim <- c(-60, 85)
  } else {
    xr <- range(bbox_pts$lon); yr <- range(bbox_pts$lat)
    xpad <- max(diff(xr) * 0.15, 1); ypad <- max(diff(yr) * 0.15, 1)
    xlim <- xr + c(-xpad, xpad); ylim <- yr + c(-ypad, ypad)
  }

  is_world  <- is.null(bbox_pts)
  basemap_r <- if (use_basemap_tiles) safe_get_tiles(xlim, ylim) else NULL

  p <- ggplot()
  p <- if (!is.null(basemap_r)) {
    p + geom_spatraster_rgb(data = basemap_r, alpha = basemap_alpha, interpolate = TRUE)
  } else {
    p + geom_sf(data = world_bounds, fill = "grey93", color = "grey65", linewidth = 0.1)
  }

  # Rasterised WLD context points: world page only. On a zoom page the
  # extent is already small enough that plotting ctx/sel as real vector
  # points (below) is cheap, and the raster haze would only obscure detail.
  if (is_world) {
    p <- p +
      # geom_spatraster(fill = <fixed color>) ignores the raster's own NA
      # mask and paints the whole extent -- go through a constant-colour
      # gradient scale instead so empty cells (na.value) stay transparent.
      geom_spatraster(data = wld_pres_r, na.rm = TRUE, alpha = vine_pts_alpha) +
      scale_fill_gradient(low = "grey15", high = "grey15", na.value = NA, guide = "none") +
      new_scale_fill()
  }

  p +
    geom_point(data = ctx, aes(lon, lat, color = region), shape = 1, size = 1.4, stroke = 0.3) +
    geom_point(data = sel, aes(lon, lat, color = region), shape = 16, size = 1.8) +
    geom_point(data = anc_coord, aes(lon, lat), shape = 8, color = "black", size = 3.5, stroke = 0.9) +
    geom_point(data = loc_coord, aes(lon, lat), shape = 4, color = "black", size = 4, stroke = 1.2) +
    scale_color_manual(values = PAL[regions_n], name = "Region analogue", drop = FALSE) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(x = NULL, y = NULL,
         title = titre %||% paste0("Points analogues -- top ", length(regions_n), " regions"),
         subtitle = paste0("croix = station ", site_label, " | etoile = point WLD apparie (",
                           round(ANC$d_km, 2), " km) | cercles pleins = points retenus,\n",
                           "cercles vides = autres points des memes regions",
                           if (is_world) " | fond gris = points WLD rasterises" else "",
                           if (!is.null(basemap_r)) " | fond de carte : (c) OpenStreetMap contributors, style OpenTopoMap (CC-BY-SA)" else "")) +
    theme_minimal()
}

reg_top   <- SYN$region
map_pages <- list(build_carte(reg_top),
                  build_carte(reg_top, rbind(PT_TOP[, .(lon, lat)], loc_coord),
                             paste0("Points analogues -- top ", length(reg_top), " regions")))
for (rg in head(reg_top, n_zoom_regions)) {
  pts <- site_coords[region == rg, .(lon, lat)]
  map_pages <- c(map_pages, list(build_carte(reg_top, pts, paste0("Zoom -- ", rg))))
}
safe_pdf_pages(paste0(out_prefix, "_carte.pdf"), map_pages, width = 12, height = 7)
message("wrote ", paste0(out_prefix, "_carte.pdf"), " (", length(map_pages), " pages)")

## ---- 13. graphiques de synthese (G.7 - G.9) -------------------------------

# -- plan factoriel ACP, axes 1 et 2 ---------------------------------------
# Fond gris allege : echantillon de 10 % des sites du monde (pas de
# rasterisation ici, graphique base R) -- les points colores (le resultat)
# restent a pleine resolution.
plan_acp <- function(res, DT, zb, SYN, PAL, bg_frac = 0.10) {
  op <- par(mar = c(7, 4, 4, 2) + 0.1)
  on.exit(par(op))
  ctr    <- DT[, lapply(.SD, mean), by = .(site, region), .SDcols = pcols]
  ctr_bg <- ctr[sample(.N, ceiling(bg_frac * .N))]
  plot(ctr_bg$PC1, ctr_bg$PC2, pch = 16, cex = .25, col = "grey85",
       xlab = paste0("PC1 (", acp$var_pct[1], " %)"),
       ylab = paste0("PC2 (", acp$var_pct[2], " %)"),
       main = "Plan factoriel -- regions analogues")
  for (rg in SYN$region) {
    sub <- ctr[region == rg]
    points(sub$PC1, sub$PC2, pch = 16, cex = .8, col = PAL[rg])
    if (nrow(sub) >= 5 && requireNamespace("ellipse", quietly = TRUE)) {
      Mm <- as.matrix(sub[, .(PC1, PC2)])
      lines(ellipse::ellipse(cov(Mm), centre = colMeans(Mm), level = 0.80),
            col = PAL[rg], lwd = 1.5)
    }
  }
  points(zb[1], zb[2], pch = 4, cex = 2.2, lwd = 3, col = "black")
  legend("topright", legend = c("cible (station)", SYN$region),
         col = c("black", PAL[SYN$region]), pch = c(4, rep(16, nrow(SYN))),
         bty = "n", cex = .7)
  mtext(paste0("axes 1-2 : ", acp$cum_pct[2], " % de la variance seulement -- ",
               "plan ILLUSTRATIF, le classement vient des ", K, " dimensions"),
        side = 1, line = 4, cex = .7)
  mtext("grey points: 10% random subsample of world wine-growing regions",
        side = 1, line = 5, cex = .65, font = 3)
}

# -- distribution des rangs par region (G.8) -------------------------------
graph_dispersion <- function(res, SYN, PAL) {
  dat <- res[region %in% SYN$region]
  dat[, region := factor(region, levels = rev(SYN$region))]
  op <- par(mar = c(4, 12, 3, 1))
  boxplot(rg_cons ~ region, data = dat, horizontal = TRUE, log = "x",
          las = 1, outline = FALSE, border = PAL[levels(dat$region)],
          xlab = "rang consensus (echelle log)", ylab = "",
          main = "Distribution des rangs, tous points de la region")
  tete <- res[region %in% SYN$region, .(mn = min(rg_cons)), by = region]
  tete[, y := match(region, levels(dat$region))]
  points(tete$mn, tete$y, pch = 18, cex = 1.6, col = "black")
  abline(v = ceiling(0.01 * nrow(res)), lty = 2, col = "red")
  legend("bottomright", c("point de tete", "seuil 1 %"), pch = c(18, NA),
         lty = c(NA, 2), col = c("black", "red"), bty = "n", cex = .7)
  par(op)
}

# -- profils mensuels du top-N (G.9) ---------------------------------------
# Unites physiques (D_phys, yb_phys) : le graphique doit rester lisible.
# Medianes des POINTS RETENUS, pas de la region entiere : on decrit le
# secteur analogue, pas un climat moyen inexistant.
# ATTENTION (F.2) : une mediane peut coller alors qu'aucun millesime
# individuel ne ressemble a la cible. Controle visuel, PAS le resultat.
UNITS <- c(Tmin = "°C", Tmax = "°C", P = "mm", ET0 = "mm")

graph_profils <- function(SYN, res, Dp, ybp, PAL, PCH, mois_lab = slot_lab) {
  # x-axis carries two rows of month labels: the northern-hemisphere calendar
  # month used to name the slots, and its southern-hemisphere equivalent
  # (6-month shift, see season_calendar()) -- the window means Apr-Aug for a
  # site in the north and Oct-Feb for a site in the south.
  cal_sud <- season_calendar(months_season, "S")

  op <- par(mfrow = c(2, 2), mar = c(6.5, 4.3, 2, 1))
  on.exit(par(op))
  for (v in VARS4) {
    cols <- paste0(v, "_", mois_lab)
    dat_reg <- lapply(SYN$region, function(rg) {
      st <- res[region == rg & rg_cons <= max(SYN$rang_fin), site]
      if (!length(st)) st <- res[region == rg][order(rg_cons)][1, site]
      unlist(Dp[site %in% st, lapply(.SD, median), .SDcols = cols])
    })
    yl <- range(unlist(dat_reg), ybp[cols], na.rm = TRUE)
    plot(seq_along(cols), ybp[cols], type = "n", ylim = yl, xaxt = "n", xlab = NA,
         ylab = paste0(v, " [", UNITS[v], "]"), main = v)
    axis(1, at = seq_along(cols), labels = month.abb[months_season], line = 0)
    axis(1, at = seq_along(cols), labels = month.abb[cal_sud$month], line = 1.6,
         tick = FALSE, font = 3, cex.axis = 0.85)
    mtext("Month North./South. hemisphere", side = 1, line = 3.4, cex = 0.75)

    for (i in seq_along(SYN$region)) {
      lines(seq_along(cols), dat_reg[[i]], col = PAL[SYN$region[i]], lwd = 1.8)
      points(seq_along(cols), dat_reg[[i]], col = PAL[SYN$region[i]],
             pch = PCH[SYN$region[i]], cex = 1)
    }
    # target vintage: dashed, so it reads as the reference curve rather than
    # one more region among the others
    lines(seq_along(cols), ybp[cols], col = "black", lwd = 3, lty = 2)
    points(seq_along(cols), ybp[cols], pch = 4, cex = 1.3, lwd = 2)

    legend("topleft", legend = c(paste0("target (", site_label, ")"), SYN$region),
           col = c("black", PAL[SYN$region]), pch = c(4, PCH[SYN$region]),
           lty = c(2, rep(1, nrow(SYN))), bty = "n", cex = 0.55, seg.len = 1.6)
  }
}

## ---- 14. rapport PDF (G.10) -----------------------------------------------

rapport_pdf <- function(chemin = file.path(dir_sorties, "rapport_analogues.pdf")) {
  pdf(chemin, width = 11, height = 8.5)

  # p.1 en-tete
  plot.new()
  text(0, 1, "Analogues climatiques de millesime", adj = 0, cex = 1.8, font = 2)
  txt <- c(
    paste0("Station : ", site_label, "  |  ", LOC_LON, " / ", LOC_LAT,
           "  |  millesime ", TARGET_YEAR),
    paste0("Point WLD apparie : ", ANC$site, " (", round(ANC$d_km, 2), " km) -- ",
           "appariement par coordonnees, point a point"),
    paste0("Region d'ancrage : ", ANC$region, " (", ANC$country, ")"),
    paste0("K = ", K, " CP, variance cumulee ", acp$cum_pct[K],
           " %  |  M0 = ", M0, " millesimes fictifs (declare, non optimise)"),
    paste0("n_proches du point d'ancrage (etalon interne) : ", n_proches_anc),
    paste0("Fenetre : mois ", paste(range(months_season), collapse = "-"),
           " -- analogues de SAISON VEGETATIVE uniquement"),
    "",
    "Estimateur principal : predictif bayesien a a priori empirique (NIW / Student).",
    "Score = statistique de classement ordinal, JAMAIS une probabilite.",
    paste0("Le top-1 parcellaire n'a aucune signification (", n_millesimes_typ,
           " millesimes en dimension ", K, ")."))
  text(0, seq(0.90, 0.35, length.out = length(txt)), txt, adj = 0, cex = .85)

  # p.2 tableau de synthese
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    grid::grid.newpage()
    grid::grid.text("Synthese -- top 10 des regions (cohortes sequentielles)",
                    y = 0.95, gp = grid::gpar(fontsize = 14, fontface = "bold"))
    gridExtra::grid.table(TAB_SYN[, .(rang, region, rang_deb, rang_fin,
                                      n_cohorte_pts, n_region_pts, pct_region,
                                      n_cohorte_pix, solidite, profil)],
                          rows = NULL)
  } else message("package gridExtra absent : tableau de synthese non insere dans le PDF")

  # p.3+ graphiques
  graph_dispersion(res, SYN, PAL)
  plan_acp(res, DT, zb, SYN, PAL)
  graph_profils(SYN, res, D_phys, yb_phys, PAL, PCH)

  # deltas de source : Tmin/Tmax partagent une echelle Y en degC (delta
  # additif station - pixel) ; P/ET0 partagent une echelle Y en % (delta
  # MULTIPLICATIF station/pixel, reexprime en biais relatif -- P et ET0 sont
  # des flux bornes a 0, un ratio est plus stable qu'une difference, voir
  # note_station_delta_carte.md C.2).
  DLT <- copy(DLT)
  DLT[, `:=`(d_P_pct = 100 * (d_P - 1), d_ET0_pct = 100 * (d_ET0 - 1))]
  yr_add <- range(DLT$d_Tmin, DLT$d_Tmax, na.rm = TRUE)
  yr_mul <- range(DLT$d_P_pct, DLT$d_ET0_pct, na.rm = TRUE)
  op <- par(mfrow = c(2, 2), mar = c(6.5, 4, 2, 1))
  for (v in ADD) {
    boxplot(as.formula(paste0("d_", v, " ~ mois")), data = DLT, ylim = yr_add,
            main = paste0("delta station - pixel : ", v), xlab = "mois",
            ylab = "delta station - pixel (°C)")
    abline(h = 0, col = "red", lty = 2)
  }
  for (v in MUL) {
    boxplot(as.formula(paste0("d_", v, "_pct ~ mois")), data = DLT, ylim = yr_mul,
            main = paste0("delta station - pixel : ", v), xlab = "mois",
            ylab = "delta station / pixel (%)")
    abline(h = 0, col = "red", lty = 2)
    mtext("delta multiplicatif (station/pixel - 1) x 100 ; P, ET0 sont des flux",
          side = 1, line = 4.2, cex = .6, font = 3)
    mtext("bornes a 0, un ratio evite l'explosion d'une difference proche de zero",
          side = 1, line = 5.0, cex = .6, font = 3)
  }
  par(op)

  # eboulis ACP
  plot(acp$var_pct[1:15], type = "b", pch = 19, xlab = "composante",
       ylab = "variance (%)", main = "Eboulis ACP (variance GLOBALE, pas locale)")
  abline(v = K, col = "red", lty = 2)

  # cartes
  for (pl in map_pages) print(pl)

  dev.off()
  message("rapport ecrit : ", chemin)
}

rapport_pdf()

## ---- 15. entrees du script de controle de continuite ----------------------
# 00_note_migration.md 3.1 : l'ancien estimateur gaussien regularise ne fait
# plus partie du pipeline. Il est execute UNE SEULE FOIS, dans le script
# separe 14b_controle_continuite_gaussien.R, qui lit ce fichier.

if (save_continuite) {
  f <- paste0(out_prefix, "_continuite_inputs.rds")
  saveRDS(list(DT_pc  = DT[, c("site", pcols), with = FALSE],
               W = W, zb = zb, K = K, pcols = pcols, M0 = M0,
               score_bay = res[, .(site, score_bay, rg_bay)]), f)
  message("wrote ", f, " -- lancer ensuite 14b_controle_continuite_gaussien.R")
}

## ---- 16. protocole de validation a deux sources (E.2) ---------------------
# A lancer APRES avoir produit deux tableaux avec cible_source = "pixel" puis
# "station", sur une annee ORDINAIRE presente dans les deux sources. Un
# millesime atypique DOIT mal classer sa propre region : le controle
# d'ancrage n'est interpretable que sur un millesime median (E.3).
#
#   comparer_runs(<..._pixel_table.csv>, <..._table.csv>)
#
# Criteres : recouvrement top-50 > 0,70 = station validee ; 0,50-0,70 =
# utilisable avec reserve, ne publier que les regions communes ; < 0,50 =
# biais de source dominant, revoir le delta ou renoncer au vecteur station.

comparer_runs <- function(path_A, path_B, N = 50L) {
  A <- fread(path_A)[, .(site, region, rgA = rg_cons)]
  B <- fread(path_B)[, .(site, rgB = rg_cons)]
  M <- merge(A, B, by = "site")
  topA <- A[order(rgA)][!duplicated(region)][1:N, region]
  topB <- merge(B, unique(A[, .(site, region)]), by = "site")[
            order(rgB)][!duplicated(region)][1:N, region]
  list(recouvrement_topN = length(intersect(topA, topB)) / N,
       spearman          = cor(M$rgA, M$rgB, method = "spearman"),
       rang_ancrage_A    = A[site == ANC$site, rgA],
       rang_ancrage_B    = B[site == ANC$site, rgB],
       communes = intersect(topA, topB), only_A = setdiff(topA, topB),
       only_B = setdiff(topB, topA))
}

## ---- 17. sensibilite aux autres parametres (8.4, manuel) ------------------
# Non automatise (rejouer les sections 6-7 est couteux sur ~940 000 lignes).
# Rejouer manuellement avec K dans c(4, 6, 8), m_knn dans c(3, 5, 8) et
# q_count dans c(0.35, 0.50, 0.65), puis comparer les recouvrements de top-200
# comme en 10.5. Conclusion a publier : une REGION d'analogues, pas une
# parcelle unique (Grenier et al. 2013, doi:10.1175/JAMC-D-12-0170.1).
