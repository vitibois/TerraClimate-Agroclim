# Note annexe — station locale, delta de source, validation et cartographie

Complement a `analogues_millesime.md`. Couvre l'usage d'une **station meteo** comme millesime cible
face a une base mondiale **TerraClimate** (pixels ~4 km), le calibrage du biais de source, le
protocole de validation, et les sorties cartographiques.

**A implementer dans cet ordre.** Chaque section depend de la precedente.

---

## A. Contexte et probleme pose

| | source | resolution | disponibilite |
|---|---|---|---|
| base de reference `WLD` | TerraClimate | pixel ~4 km | differee (millesime N-1 au mieux) |
| millesime cible `LOC` | station meteo | ponctuelle | quasi temps reel |

La station permet de diagnostiquer un millesime en cours (2026) que TerraClimate ne couvre pas
encore. Mais **un pixel de 4 km et une station ne mesurent pas la meme chose** : la station est
ponctuelle, souvent plus ouverte ou plus urbanisee, et le pixel est une moyenne spatiale. Le biais
attendu porte surtout sur **Tmax** et sur **ET0**.

Ce biais entre dans toutes les distances, y compris celle du site a lui-meme. Non corrige, il
deplace le classement entier de facon systematique et invisible.

> **Regle.** Ne jamais injecter un vecteur station dans le pipeline sans avoir mesure le delta sur
> les annees communes (section C) et verifie sa systematicite (section D).

---

## A.1 Precision de vocabulaire

Ce qui est calcule n'est **pas une classification**. Il n'y a pas de partition, pas de groupes a
trouver, pas de CAH. C'est un **score de plausibilite** calcule site par site : pour chaque
vignoble, le millesime cible tombe-t-il dans la zone ou ce vignoble produit habituellement ses
millesimes ?

Le mot « gaussien » ne qualifie que la facon de decrire la **forme** du nuage de 20 millesimes
(centre + etirement directionnel). L'asymetrie des variables est traitee en amont par les scores
normaux (§2.2 du document principal), et le kNN sert de garde-fou sans hypothese de forme.

---

## B. Parametres d'ancrage — coordonnees, pas nom de region

La station est reperee par ses **coordonnees**, et on retient le **point de `WLD` le plus proche
geographiquement**. Aucun nom de region a saisir : la region est deduite du point apparie.

`WLD` contient des points avec leurs coordonnees propres (le climat vient d'un raster, mais les
lignes sont des points). Il s'agit donc d'un appariement **point a point**, pas d'une extraction
raster.

Avec ~1 point par kilometre carre, le point apparie est tres proche, et la comparaison
station / point est bien plus stricte qu'une moyenne regionale.

```r
# ===========================================================================
# PARAMETRES  (a renseigner)
# ===========================================================================
LOC_LON     <- -0.5800     # longitude de la station, degres decimaux
LOC_LAT     <- 44.8300     # latitude
TARGET_YEAR <- 2026L
D_MAX_KM    <- 15          # distance d'appariement maximale toleree

# colonnes reelles de WLD  -- A VERIFIER avec names(WLD)
COL_SITE   <- "site"
COL_REGION <- "WINE_REGION"
COL_COUNTRY<- "CNT"
COL_LON    <- "lon"
COL_LAT    <- "lat"
```

### B.1 Appariement geographique

```r
# distance grand-cercle (haversine), en km
dist_km <- function(lon1, lat1, lon2, lat2) {
  R <- 6371
  p1 <- lat1 * pi/180; p2 <- lat2 * pi/180
  dp <- p2 - p1; dl <- (lon2 - lon1) * pi/180
  a  <- sin(dp/2)^2 + cos(p1) * cos(p2) * sin(dl/2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}

apparier_station <- function(WLD) {
  PTS <- unique(WLD[, .(site    = get(COL_SITE),
                        region  = get(COL_REGION),
                        country = get(COL_COUNTRY),
                        lon     = get(COL_LON),
                        lat     = get(COL_LAT))])

  PTS[, d_km := dist_km(LOC_LON, LOC_LAT, lon, lat)]
  setorder(PTS, d_km)

  if (PTS$d_km[1] > D_MAX_KM)
    stop("Point WLD le plus proche a ", round(PTS$d_km[1], 1),
         " km (> D_MAX_KM = ", D_MAX_KM, "). Verifier lon/lat et leur ordre.")

  anc <- as.list(PTS[1])
  message("Station appariee au point ", anc$site, " | ", anc$region,
          " (", anc$country, ") | ", round(anc$d_km, 2), " km")
  message("  5 points les plus proches : ",
          paste0(PTS$site[1:5], " (", round(PTS$d_km[1:5], 1), " km)",
                 collapse = " | "))

  anc$PTS <- PTS
  anc
}

ANC <- apparier_station(WLD)
```

> **Piege classique.** Une inversion longitude / latitude produit un point plausible mais faux. Le
> garde-fou `D_MAX_KM` l'attrape. Verifier aussi le signe de la longitude a l'ouest de Greenwich.

### B.2 Voisinage de controle

Le point apparie peut etre atypique par accident (relief local). Les voisins immediats servent de
controle : leur climat doit etre proche de celui du point retenu.

```r
VOIS <- ANC$PTS[d_km <= 10, site]
message("voisinage de controle : ", length(VOIS), " points dans un rayon de 10 km")
```
## C. Delta de source, calibre sur les annees communes

### C.1 Principe

`LOC` est fourni comme un jeu de donnees **pluriannuel** (journalier ou mensuel). On ne corrige pas
a l'aveugle : on mesure l'ecart station / pixel **sur l'intersection des annees**, mois par mois,
puis on applique ce delta au millesime cible — que celui-ci soit dans `WLD` ou non.

Il suffit de quelques annees communes pour calibrer. Cinq est un minimum tres faible, dix est
confortable.

### C.2 Forme du delta selon la variable

| variable | forme | justification |
|---|---|---|
| `Tmin`, `Tmax` | **additive** (station − pixel) | echelle d'intervalle, pas de borne |
| `P` | **multiplicative** (station / pixel) | flux borne a 0, ratio plus stable que la difference |
| `ET0` | **multiplicative** | idem, flux borne a 0 |

Un ratio sur un mois quasi sec explose. Deux garde-fous : plancher sur le denominateur, et bascule
sur un ratio calcule au niveau du **cumul saisonnier** quand les valeurs mensuelles sont trop
faibles.

### C.3 Agregation mensuelle de LOC

```r
# LOC journalier -> mensuel. Si LOC est deja mensuel, sauter cette etape.
agreger_loc <- function(LOC_daily) {
  L <- as.data.table(LOC_daily)
  L[, `:=`(year = year(date), mois = month(date))]
  L[, .(Tmin = mean(tmin, na.rm = TRUE),
        Tmax = mean(tmax, na.rm = TRUE),
        P    = sum(prec,  na.rm = TRUE),
        ET0  = sum(et0,   na.rm = TRUE),
        n_j  = .N),
    by = .(year, mois)]
}

LOC_m <- agreger_loc(LOC)

# controle de completude : un mois lacunaire biaise cumul P et ET0
jours_att <- c(31,28,31,30,31,30,31,31,30,31,30,31)
LOC_m[, compl := n_j / jours_att[mois]]
if (LOC_m[compl < 0.95, .N] > 0) {
  warning("mois incomplets (<95 % de jours) :")
  print(LOC_m[compl < 0.95, .(year, mois, n_j, compl = round(compl, 2))])
}
LOC_m <- LOC_m[compl >= 0.95]
```

### C.4 Calcul du delta

```r
VARS4 <- c("Tmin", "Tmax", "P", "ET0")
ADD   <- c("Tmin", "Tmax")          # delta additif
MUL   <- c("P", "ET0")              # delta multiplicatif
P_MIN <- 5                          # mm, plancher du denominateur

calibrer_delta <- function(LOC_m, WLD, ANC, mois = 4:8) {

  # POINT APPARIE (jamais la moyenne regionale) : comparaison stricte
  px_m <- WLD[get(COL_SITE) == ANC$site,
              c(list(year = year, mois = mois), .SD), .SDcols = VARS4]
  if (!nrow(px_m)) stop("Point apparie ", ANC$site, " absent de WLD.")

  # variante de robustesse : moyenne du voisinage 10 km (B.2), a comparer
  # px_m <- WLD[get(COL_SITE) %in% VOIS,
  #             lapply(.SD, mean), by = .(year, mois), .SDcols = VARS4]

  comm <- intersect(LOC_m$year, px_m$year)
  if (length(comm) < 5L)
    stop("Seulement ", length(comm), " annee(s) commune(s) : calibrage impossible.")
  message("Calibrage sur ", length(comm), " annees communes : ",
          paste(range(comm), collapse = "-"))

  M <- merge(LOC_m[year %in% comm & mois %in% mois],
             px_m[year %in% comm & mois %in% mois],
             by = c("year", "mois"), suffixes = c("_loc", "_px"))

  # deltas annuels, par mois
  for (v in ADD) M[[paste0("d_", v)]] <- M[[paste0(v, "_loc")]] - M[[paste0(v, "_px")]]
  for (v in MUL) M[[paste0("d_", v)]] <- M[[paste0(v, "_loc")]] /
                                          pmax(M[[paste0(v, "_px")]], P_MIN)
  M[]
}

DLT <- calibrer_delta(LOC_m, WLD, ANC)
```

---

## D. Systematicite du delta — corriger ou ecarter

Un biais **stable** se corrige sans hesiter. Un ecart **erratique** signale que la station n'est pas
representative de son pixel, et la correction ajoute alors du bruit. Il faut trancher explicitement.

### D.1 Trois criteres

1. **Rapport signal / bruit** : `|delta moyen| / ecart-type du delta`. Eleve = biais net et stable.
2. **Coefficient de variation du delta** : dispersion relative. Faible = reproductible.
3. **Test de tendance** sur le delta en fonction de l'annee : une pente significative revele une
   derive (changement de capteur, deplacement de station, urbanisation), ce qui interdit
   l'extrapolation naive au millesime cible.

### D.2 Code

```r
qualifier_delta <- function(DLT) {
  out <- rbindlist(lapply(VARS4, function(v) {
    col <- paste0("d_", v)
    DLT[, {
      x  <- get(col); n <- sum(is.finite(x)); x <- x[is.finite(x)]
      mu <- mean(x);  sdv <- sd(x)
      # tendance temporelle du delta
      tt <- if (n >= 6L) {
              ct <- suppressWarnings(cor.test(year[is.finite(get(col))], x,
                                              method = "spearman"))
              c(rho = unname(ct$estimate), p = ct$p.value)
            } else c(rho = NA_real_, p = NA_real_)
      .(variable = v, n = n,
        delta    = round(mu, 3),
        sd       = round(sdv, 3),
        snr      = round(abs(mu) / sdv, 2),
        rho_an   = round(tt["rho"], 2),
        p_tend   = signif(tt["p"], 2))
    }, by = mois]
  }))

  out[, fiabilite := fifelse(!is.na(p_tend) & p_tend < 0.05 & abs(rho_an) > 0.5,
                             "DERIVE - ne pas extrapoler",
                     fifelse(snr >= 2, "fiable",
                     fifelse(snr >= 1, "reserve", "ECARTER - non systematique")))]
  setorder(out, variable, mois)
  out[]
}

QD <- qualifier_delta(DLT)
print(QD)

message("--- synthese ---")
print(QD[, .N, by = fiabilite])
if (QD[grepl("ECARTER|DERIVE", fiabilite), .N] > 0)
  warning("Certains couples variable x mois ne sont pas corrigibles de facon fiable. ",
          "Voir colonne 'fiabilite'. Les resultats du millesime cible sur ces mois ",
          "doivent etre presentes avec reserve explicite.")

# visualisation : delta par mois, dispersion interannuelle
par(mfrow = c(2, 2))
for (v in VARS4) {
  boxplot(as.formula(paste0("d_", v, " ~ mois")), data = DLT,
          main = v, xlab = "mois", ylab = "delta station - pixel")
  abline(h = if (v %in% ADD) 0 else 1, col = "red", lty = 2)
}
par(mfrow = c(1, 1))
```

### D.3 Application au millesime cible

```r
appliquer_delta <- function(LOC_m, QD, an = TARGET_YEAR, mois = 4:8,
                            mode = c("median", "moyenne")) {
  mode <- match.arg(mode)
  cible <- LOC_m[year == an & mois %in% mois]
  if (nrow(cible) != length(mois))
    stop("Millesime ", an, " incomplet dans LOC : ", nrow(cible), "/", length(mois), " mois.")

  D_ <- copy(cible)
  for (v in VARS4) {
    for (m in mois) {
      d <- QD[variable == v & mois == m, delta]
      if (v %in% ADD) D_[mois == m, (v) := get(v) - d]     # station -> echelle pixel
      else            D_[mois == m, (v) := get(v) / d]
    }
  }
  # mise au format large attendu par le pipeline : Tmin_04 ... ET0_08
  yb <- unlist(lapply(VARS4, function(v)
          setNames(D_[order(mois)][[v]], paste0(v, "_", sprintf("%02d", mois)))))
  yb[VARS]
}

yb_corr <- appliquer_delta(LOC_m, QD)
yb_brut <- appliquer_delta(LOC_m, QD[, .(variable, mois,
                                         delta = fifelse(variable %in% ADD, 0, 1))])
```

> **Sens de la correction.** On ramene la station **vers l'echelle pixel**, jamais l'inverse : la
> base de reference de 940 000 lignes ne doit pas etre touchee.

> **Toujours conserver `yb_brut`.** Faire tourner le pipeline sur les deux vecteurs et comparer les
> classements : c'est la mesure directe de l'effet de la correction (section E.2).

---

## E. Protocole de validation a deux sources

Le test decisif. Il valide — ou invalide — l'usage de la station pour les millesimes que
TerraClimate ne couvre pas encore.

### E.1 Principe

Prendre une annee presente dans **les deux** sources (p. ex. 2025, apres extension de `WLD` a
2001-2025). Faire tourner le pipeline deux fois, tout identique sauf le vecteur cible :

- **Run A** : millesime issu du pixel TerraClimate de la region.
- **Run B** : millesime issu de la station, corrige du delta.

Tout le reste est constant, donc l'ecart observe mesure exactement l'effet de source.

### E.2 Metriques a rapporter

```r
comparer_runs <- function(resA, resB, N = 50L) {
  key <- c("site", "region")
  A <- resA[, .(site, region, rgA = rg_cons)]
  B <- resB[, .(site, region, rgB = rg_cons)]
  M <- merge(A, B, by = key)

  topA <- A[order(rgA)][!duplicated(region)][1:N, region]
  topB <- B[order(rgB)][!duplicated(region)][1:N, region]

  list(
    recouvrement_topN = length(intersect(topA, topB)) / N,
    spearman          = cor(M$rgA, M$rgB, method = "spearman"),
    rang_ancrage_A    = A[region == ANC$region, min(rgA)],
    rang_ancrage_B    = B[region == ANC$region, min(rgB)],
    communes          = intersect(topA, topB),
    only_A            = setdiff(topA, topB),
    only_B            = setdiff(topB, topA)
  )
}

VAL <- comparer_runs(res_pixel, res_station)
str(VAL[1:4])
```

### E.3 Criteres de decision

| recouvrement top-50 | lecture |
|---|---|
| > 0,70 | station validee, utilisable pour les millesimes hors `WLD` |
| 0,50 – 0,70 | utilisable avec reserve ; ne publier que les regions communes aux deux runs |
| < 0,50 | biais de source dominant ; revoir le delta ou renoncer au vecteur station |

**Second critere, independant du premier : le rang de la region d'ancrage.** Si le vecteur cible est
coherent avec `WLD`, le millesime doit ressortir comme relativement typique de sa propre region, ou
au moins de regions voisines. Un rang de 30 000 sur 47 000 signale un probleme d'unites, de source
ou d'alignement — pas un resultat climatique.

> **Nuance importante pour un millesime atypique.** Un millesime « ovni » (2026 a Bordeaux) *doit*
> mal classer sa propre region : c'est le resultat attendu, pas une erreur. Le controle d'ancrage
> n'est interpretable que sur un millesime **ordinaire**. Faire la validation E.1 sur une annee
> mediane, jamais sur le millesime extreme qu'on cherche a diagnostiquer.

---

## F. Interpretation : « a combien de millesimes ressemble-t-il ? »

C'est la colonne `n_proches` du tableau principal, et elle repond directement a la question.

`n_proches = 7` sur `n_tot = 20` se lit : **dans ce vignoble, 7 millesimes sur 20 — un sur trois —
ressemblent au millesime cible**. C'est la statistique de communication, celle qui s'ecrit dans un
article.

### F.1 Le seuil est un curseur, a expliciter

```r
# r = sqrt(qchisq(q, df = K)) : q fixe la severite
for (q in c(0.25, 0.35, 0.50, 0.65)) {
  rq <- sqrt(qchisq(q, df = K))
  cnt <- DTw[, .(n = sum(d <= rq)), by = .(site, region)]
  message("q = ", q, " (r = ", round(rq, 2), ") | n_proches median = ",
          median(cnt$n), " | ancrage = ",
          cnt[region == ANC$region, round(median(n), 1)])
}
```

Serrer `q` repond a « combien de millesimes ressemblent *de tres pres* », desserrer repond a
« combien sont dans la meme famille ». Toujours publier la valeur retenue.

### F.2 Ne pas confondre avec le graphique des profils mensuels

Le graphique qui superpose le millesime cible a la **mediane** des 20 ans du site est un controle
visuel, **pas** le resultat. Une mediane peut coller parfaitement alors qu'aucun millesime
individuel ne ressemble a la cible : il suffit que les millesimes se repartissent de part et d'autre.
C'est le piege de la moyenne, et la raison pour laquelle tout le pipeline travaille sur les
millesimes individuels.

Sortie recommandee : mediane **et** enveloppe interquartile des 20 millesimes, cible en surimpression.

```r
profil_regional <- function(reg, res, D, yb, mois = 4:8) {
  st <- res[region == reg, site]
  S  <- D[site %in% st]
  rbindlist(lapply(VARS4, function(v) {
    cols <- paste0(v, "_", sprintf("%02d", mois))
    rbindlist(lapply(seq_along(cols), function(i) {
      x <- S[[cols[i]]]
      .(variable = v, mois = mois[i],
        q25 = quantile(x, .25), med = median(x), q75 = quantile(x, .75),
        cible = yb[cols[i]])
    }))
  }))
}
```

### F.3 Etalon interne obligatoire

`n_proches` ne veut rien dire dans l'absolu. Il se lit **relativement au `n_proches` de la region
d'ancrage** (§8.1 du document principal). Si le millesime cible ne ressemble qu'a 2 millesimes de sa
propre region, une region etrangere a 7 est un resultat fort. S'il ressemble a 12 des siens, il n'a
rien de singulier et le cadrage change.

---

## G. Agregation par region et sorties de synthese

### G.1 Le probleme du nombre de points par region

Les points sont a ~1 par km2. Une grande region peut en compter des centaines, une petite quelques
unites. Un top-10 brut de points peut donc etre entierement occupe par une seule region.

**Principe retenu.** La diversite climatique interne d'une region est une **propriete reelle** de
cette region, pas du bruit. On ne moyenne donc jamais, et on ne resume pas une region a un
centroide. On construit un classement de regions qui **respecte l'ordre des points**.

### G.2 Cohortes sequentielles

On parcourt les points par rang croissant. Chaque fois qu'une **nouvelle** region apparait, elle
prend le rang regional suivant, et sa **cohorte** est l'ensemble de ses points contigus dans ce
classement.

Exemple : les rangs 1 a 50 sont tous dans la region X → X est la region n°1, cohorte de 50 points,
rang de fin 50. Les rangs 51 a 53 sont dans Y → Y est la region n°2, cohorte de 3 points, rang de
fin 53. Etc.

Ce qui donne, pour chaque region du top-10 : ou commence sa cohorte, ou elle finit, combien de
points elle contient, et **quelle part de la region** cela represente.

```r
construire_cohortes <- function(res, WLD, N_REG = 10L) {
  setorder(res, rg_cons)
  R <- copy(res)[, rang_pt := .I]

  # effectif total de points par region dans WLD
  eff <- unique(WLD[, .(site = get(COL_SITE), region = get(COL_REGION))]
               )[, .(n_region = .N), by = region]

  # cohortes contigues : nouveau bloc a chaque changement de region
  R[, bloc := rleid(region)]

  COH <- R[, .(region     = region[1],
               rang_deb   = min(rang_pt),
               rang_fin   = max(rang_pt),
               n_cohorte  = .N,
               rg_med     = median(rg_cons),
               solidite   = solidite[1]),
           by = bloc]

  # une region peut reapparaitre plus loin : on garde son premier bloc,
  # mais on totalise ses points sur l'ensemble du top considere
  COH <- COH[!duplicated(region)]
  setorder(COH, rang_deb)
  COH[, rang_region := .I]

  COH <- merge(COH, eff, by = "region", all.x = TRUE)
  COH[, pct_region := round(100 * n_cohorte / n_region, 1)]
  setorder(COH, rang_region)
  head(COH, N_REG)[]
}

COH <- construire_cohortes(res, WLD)
```

### G.3 Dispersion des points retenus au sein de la region

Deux mesures complementaires, toutes deux en **ecarts de millesime** (unite de la metrique poolee),
donc comparables entre regions.

1. **Representativite du point de tete** : son rang par rapport a la distribution des rangs de tous
   les points de sa region. Le point de tete est-il typique de sa region, ou une singularite isolee ?
2. **Etalement climatique de la region** : dispersion des points de la region dans l'espace blanchi.

```r
mesurer_dispersion <- function(COH, res, DTw, WLD) {

  # 1. distribution des rangs, tous points de la region
  rgs <- res[region %in% COH$region,
             .(rg_med_reg  = median(rg_cons),
               rg_q25      = quantile(rg_cons, .25),
               rg_q75      = quantile(rg_cons, .75),
               n_top1pct   = sum(rg_cons <= ceiling(0.01 * nrow(res)))),
             by = region]

  # 2. etalement climatique : dispersion des centres de points, espace blanchi
  ctr <- DTw[, lapply(.SD, mean), by = .(site, region), .SDcols = wcols]
  spr <- ctr[region %in% COH$region, {
      M  <- as.matrix(.SD)
      mu <- colMeans(M)
      dd <- sqrt(rowSums(sweep(M, 2, mu, "-")^2))
      .(etalement = round(median(dd), 2), etal_p90 = round(quantile(dd, .9), 2))
    }, by = region, .SDcols = wcols]

  out <- merge(merge(COH, rgs, by = "region"), spr, by = "region")

  # flag de lecture
  out[, profil := fifelse(pct_region >= 50, "region homogene et proche",
                  fifelse(pct_region >= 15, "part substantielle de la region",
                  fifelse(n_cohorte >= 5,   "sous-secteur localise",
                                            "point isole - a signaler")))]
  setorder(out, rang_region)
  out[]
}

SYN <- mesurer_dispersion(COH, res, DTw, WLD)
```

### G.4 Tableau de synthese — sortie principale

```r
TAB_SYN <- SYN[, .(rang        = rang_region,
                   region, country = NA_character_,
                   rang_deb, rang_fin,
                   n_cohorte, n_region, pct_region,
                   rg_med_reg, rg_q25, rg_q75,
                   etalement, etal_p90,
                   n_proches   = NA_integer_,
                   solidite, profil)]
print(TAB_SYN)

fwrite(TAB_SYN, "sorties/synthese_top10_regions.csv")
fwrite(res[rg_cons <= max(SYN$rang_fin)], "sorties/points_retenus.csv")
```

Lecture : `pct_region` a 80 % dit que la region entiere ressemble a la station. A 2 % avec
`n_cohorte` de 3, l'analogie porte sur un sous-secteur precis — resultat valable, mais a formuler
comme tel.

### G.5 Palette commune a toutes les sorties

Une seule palette, reutilisee sur la carte, l'ACP, les profils mensuels et le tableau.

```r
PAL <- c(ancrage = "#111111",
         setNames(hcl.colors(nrow(SYN), "Dark 3"), SYN$region))
SYN[, couleur := PAL[region]]
```

### G.6 Carte — points reels, aucun centroide

| couche | contenu | symbole |
|---|---|---|
| station | `LOC_LON` / `LOC_LAT` | croix noire |
| point apparie | `ANC$site` | etoile noire |
| points retenus | cohortes du top-10 | cercles pleins, couleur `PAL` |
| autres points des memes regions | contexte | cercles vides, meme couleur |
| fond | ensemble de `WLD` | points gris fins |

```r
PT_TOP <- merge(res[rg_cons <= max(SYN$rang_fin), .(site, region, rg_cons)],
                unique(WLD[, .(site = get(COL_SITE), lon = get(COL_LON),
                               lat  = get(COL_LAT))]), by = "site")
PT_TOP[, couleur := PAL[region]]
fwrite(PT_TOP, "sorties/points_carte.csv")
```

Zoom par region : emprise = enveloppe des points de la region, plus marge. On lit directement si les
points retenus forment un secteur coherent ou un semis disperse.

### G.7 Plan ACP — axes 1 et 2

Meme palette. Le point cible, les points retenus, et une ellipse de dispersion par region.

```r
plan_acp <- function(res, DT, zb, SYN, PAL) {
  ctr <- DT[, lapply(.SD, mean), by = .(site, region), .SDcols = pcols]

  plot(ctr$PC1, ctr$PC2, pch = 16, cex = .25, col = "grey85",
       xlab = paste0("PC1 (", acp$var_pct[1], " %)"),
       ylab = paste0("PC2 (", acp$var_pct[2], " %)"),
       main = "Plan factoriel — regions analogues")

  for (rg in SYN$region) {
    sub <- ctr[region == rg]
    points(sub$PC1, sub$PC2, pch = 16, cex = .8, col = PAL[rg])
    if (nrow(sub) >= 5) {
      ee <- ellipse::ellipse(cov(as.matrix(sub[, .(PC1, PC2)])),
                             centre = colMeans(as.matrix(sub[, .(PC1, PC2)])),
                             level = 0.80)
      lines(ee, col = PAL[rg], lwd = 1.5)
    }
  }
  # point cible
  points(zb[1], zb[2], pch = 4, cex = 2.2, lwd = 3, col = "black")
  legend("topright", legend = c("cible (station)", SYN$region),
         col = c("black", PAL[SYN$region]),
         pch = c(4, rep(16, nrow(SYN))), bty = "n", cex = .7)
}
```

> Les axes 1 et 2 ne portent qu'une partie de la variance (voir `acp$cum_pct[2]`). Deux points
> proches dans ce plan peuvent etre eloignes sur les axes 3 a 6. Le plan est **illustratif**, le
> classement vient des K dimensions.

### G.8 Distribution des rangs par region

Le graphique qui repond a « ce point est-il le coeur ou l'extremite de sa region ? ».

```r
graph_dispersion <- function(res, SYN, PAL) {
  dat <- res[region %in% SYN$region]
  dat[, region := factor(region, levels = rev(SYN$region))]

  boxplot(rg_cons ~ region, data = dat, horizontal = TRUE, log = "x",
          las = 1, outline = FALSE, border = PAL[levels(dat$region)],
          xlab = "rang consensus (echelle log)", ylab = "",
          main = "Distribution des rangs, tous points de la region")

  # point de tete de chaque region
  tete <- res[region %in% SYN$region, .(mn = min(rg_cons)), by = region]
  tete[, y := match(region, levels(dat$region))]
  points(tete$mn, tete$y, pch = 18, cex = 1.6, col = "black")
  abline(v = ceiling(0.01 * nrow(res)), lty = 2, col = "red")
  legend("bottomright", c("point de tete", "seuil 1 %"),
         pch = c(18, NA), lty = c(NA, 2), col = c("black", "red"),
         bty = "n", cex = .7)
}
```

Boite etroite et a gauche : region homogene et globalement proche. Boite large avec losange
tres a gauche : l'analogie tient a un secteur, la region dans son ensemble est plus eloignee.

### G.9 Profils mensuels du top-10

```r
graph_profils <- function(SYN, res, D, yb, PAL, mois = 4:8) {
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  for (v in VARS4) {
    cols <- paste0(v, "_", sprintf("%02d", mois))
    yl <- range(unlist(D[, ..cols]), yb[cols], na.rm = TRUE)
    plot(mois, yb[cols], type = "n", ylim = yl, xlab = "mois", ylab = v, main = v)
    for (rg in SYN$region) {
      st <- res[region == rg & rg_cons <= max(SYN$rang_fin), site]
      if (!length(st)) st <- res[region == rg][order(rg_cons)][1, site]
      M <- D[site %in% st, lapply(.SD, median), .SDcols = cols]
      lines(mois, unlist(M), col = PAL[rg], lwd = 1.8)
    }
    lines(mois, yb[cols], col = "black", lwd = 3, lty = 1)
    points(mois, yb[cols], pch = 4, cex = 1.3, lwd = 2)
  }
  par(mfrow = c(1, 1))
}
```

Les courbes regionales sont les **medianes des points retenus**, pas de la region entiere : elles
decrivent le secteur analogue, pas un climat moyen inexistant.

### G.10 Rapport PDF

```r
rapport_pdf <- function(chemin = "sorties/rapport_analogues.pdf") {
  dir.create(dirname(chemin), showWarnings = FALSE, recursive = TRUE)
  pdf(chemin, width = 11, height = 8.5)

  # p.1 en-tete
  plot.new()
  text(0, 1, "Analogues climatiques de millesime", adj = 0, cex = 1.8, font = 2)
  txt <- c(
    paste0("Station : ", LOC_LON, " / ", LOC_LAT, "  |  millesime ", TARGET_YEAR),
    paste0("Point WLD apparie : ", ANC$site, " (", round(ANC$d_km, 2), " km)"),
    paste0("Region d'ancrage : ", ANC$region, " (", ANC$country, ")"),
    paste0("K = ", K, " CP, variance cumulee ", acp$cum_pct[K], " %  |  lambda = ", lambda),
    paste0("n_proches de la region d'ancrage (leave-one-out) : ", n_proches_anc),
    paste0("Fenetre : mois ", paste(range(MOIS), collapse = "-"),
           " — analogues de SAISON VEGETATIVE uniquement")
  )
  text(0, seq(0.85, 0.55, length.out = length(txt)), txt, adj = 0, cex = .9)

  # p.2 tableau de synthese
  plot.new(); title("Synthese — top 10 des regions")
  gridExtra::grid.table(TAB_SYN[, .(rang, region, rang_deb, rang_fin,
                                    n_cohorte, n_region, pct_region, profil)],
                        rows = NULL)

  # p.3+ graphiques
  graph_dispersion(res, SYN, PAL)
  plan_acp(res, DT, zb, SYN, PAL)
  graph_profils(SYN, res, D, yb, PAL)

  # diagnostics
  plot(grid, rowSums(LL), type = "b", xlab = "lambda",
       ylab = "log-vraisemblance CV", main = "Choix de lambda")
  plot(acp$var_pct[1:15], type = "b", pch = 19, xlab = "composante",
       ylab = "variance (%)", main = "Eboulis ACP"); abline(v = K, col = "red", lty = 2)

  dev.off()
  message("rapport ecrit : ", chemin)
}

rapport_pdf()
```

Packages : `gridExtra` pour le tableau, `ellipse` pour les ellipses ACP. Alternative plus soignee :
sortie Quarto (`.qmd`) rendue en PDF, si vous voulez du texte redige autour des figures.
## H. Generalisation — n'importe quel millesime, n'importe quelle region

Le pipeline ne comporte rien de specifique a Bordeaux. Pour diagnostiquer un autre millesime :

```r
diagnostiquer <- function(country, region, an, WLD, LOC_m = NULL) {
  LOC_COUNTRY <<- country; LOC_REGION <<- region; TARGET_YEAR <<- an
  ANC <- resoudre_ancrage(WLD)

  if (is.null(LOC_m)) {
    # cible tiree directement de WLD : aucun delta a appliquer
    px <- WLD[get(COL_COUNTRY) == country & get(COL_REGION) == region & year == an]
    if (!nrow(px)) stop("Millesime ", an, " absent de WLD pour ", region)
    yb <- unlist(px[, lapply(.SD, mean), .SDcols = VARS])
  } else {
    QD <- qualifier_delta(calibrer_delta(LOC_m, WLD, ANC))
    yb <- appliquer_delta(LOC_m, QD, an = an)
  }
  list(ancrage = ANC, yb = yb[VARS])
}
```

### H.1 Point de vigilance — exclusion de la cible

Si le millesime cible est **dans** `WLD`, il participe a l'ACP, a la metrique poolee `W` et au nuage
de sa propre region. Sur 940 000 lignes, l'effet sur l'ACP et sur `W` est negligeable. Sur le nuage
de sa region (20 millesimes dont lui), il ne l'est pas : la region d'ancrage est mecaniquement
avantagee.

Pour un controle d'ancrage honnete, recalculer le score de la region d'ancrage **en retirant l'annee
cible** de son propre nuage :

```r
score_ancrage_loo <- function(DT, ANC, an, zb, W, lambda, K) {
  M <- as.matrix(DT[region == ANC$region & year != an, ..pcols])
  S <- lambda * W + (1 - lambda) * cov(M)
  R <- chol(S)
  x <- backsolve(R, zb - colMeans(M), transpose = TRUE)
  -sum(log(diag(R))) - 0.5 * sum(x^2) - 0.5 * K * log(2 * pi)
}
```

### H.2 Interpretation croisee

Une lecture en deux axes, a produire systematiquement :

| | typique chez lui | atypique chez lui |
|---|---|---|
| **analogues nombreux ailleurs** | millesime ordinaire, region banale | millesime extreme, mais existant ailleurs → *analogue de futur* |
| **peu d'analogues** | region climatiquement isolee | millesime sans precedent connu |

La case en haut a droite est la plus interessante pour un article : un millesime hors norme
localement, mais qui correspond au regime habituel d'une autre region. C'est le cas 2026 a Bordeaux
ressortant en Mediterranee chaude et seche.

---

## I. Checklist avant interpretation

- [ ] `LOC_LON` / `LOC_LAT` corrects, ordre et signe verifies (section B.1)
- [ ] Distance d'appariement `< D_MAX_KM`, et coherente avec le semis (~1 km attendu)
- [ ] Noms de colonnes reels verifies (`COL_SITE`, `COL_REGION`, `COL_LON`, `COL_LAT`)
- [ ] Region du point apparie conforme a l'attendu (controle de bon sens)
- [ ] Voisinage 10 km inspecte : le point apparie n'est pas un accident de relief (B.2)
- [ ] Mois incomplets de `LOC` ecartes (section C.3)
- [ ] Au moins 5 annees communes, 10 de preference (section C.4)
- [ ] `QD` inspecte : aucun couple variable x mois en `ECARTER` ou `DERIVE` non documente
- [ ] Delta additif sur `Tmin`/`Tmax`, multiplicatif sur `P`/`ET0`
- [ ] Delta calibre sur le **point apparie**, pas sur la moyenne regionale
- [ ] Correction appliquee **vers** l'echelle du point WLD, `WLD` intacte
- [ ] Pipeline lance sur `yb_brut` **et** `yb_corr`, ecart mesure
- [ ] Validation deux sources faite sur une annee **mediane**, pas sur le millesime extreme
- [ ] Score du point d'ancrage calcule en leave-one-out (section H.1)
- [ ] `n_proches` du point d'ancrage releve comme etalon (section F.3)
- [ ] Cohortes verifiees : `pct_region` et `profil` coherents (section G.2-G.3)
- [ ] Aucune region du top-10 resumee par un centroide sur la carte (section G.6)
- [ ] Recouvrement des 4 colonnes de rang releve (§5.5-5.6 du document principal)

---

## J. Rappel — ce que la fenetre avril-aout ne dit pas

Ces analogues sont des **analogues de saison vegetative**, pas des analogues climatiques complets.
Aucune contrainte hivernale n'entre dans le calcul : pas de besoins en froid, pas de risque de gel
hivernal, pas de dormance. Un site tropical d'altitude peut donc ressortir comme excellent analogue
tout en etant viticolement inviable.

A mentionner explicitement dans toute publication issue de ce pipeline.
