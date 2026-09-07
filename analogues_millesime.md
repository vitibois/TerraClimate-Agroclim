# Analogues climatiques de millésime — spécification méthodologique

**Phase 1 : niveau parcellaire (local).** La phase 2 régionale est décrite au §9 mais n'est pas
implémentée ici.

**Question de recherche.** Dans quels vignobles du monde trouve-t-on *fréquemment* des millésimes
dont le profil climatique de saison végétative est très proche de celui de Bordeaux 2025 ?

Autrement dit : de quel site Bordeaux 2025 serait-il une **année typique** ?

---

## 0. Contexte pour l'agent

Ce document est une spécification, pas un tutoriel. Il contient le cadre statistique, les décisions
à respecter, le code de référence et les diagnostics obligatoires. Le code est fonctionnel mais doit
être adapté aux noms de colonnes réels.

### 0.1 Cadrage statistique

L'estimand est la **densité de la distribution des millésimes du site, évaluée au point
Bordeaux 2025** :

$$\text{score}(s) = f_s(x_0)$$

Tout le reste est un choix d'estimateur. Trois estimateurs sont calculés en parallèle et leur
concordance constitue l'argument de robustesse.

### 0.2 Non négociable

1. On travaille sur les **millésimes individuels**, jamais sur les moyennes 20 ans par site.
2. **Retrait de tendance obligatoire** avant tout calcul (§2.1). Sans lui, la variabilité intra-site
   est gonflée par le réchauffement et « millésime typique » mélange deux époques.
3. La métrique commune est la **covariance intra-site poolée** `W`, estimée sur ~893 000 degrés de
   liberté. C'est la métrique intra-classe de l'analyse discriminante linéaire.
4. La troncature ACP est une **régularisation**, pas une commodité.
5. Le classement principal est la **log-densite predictive bayesienne (NIW / Student)**. Le comptage
   `n_proches` est une statistique de **communication**, pas de classement.

### 0.3 Ouvert, à tester et rapporter

`K` (nb de composantes), `M0` (poids d'a priori, à **déclarer** et non optimiser),
`h` (fenêtre KDE), `m` (rang kNN), `r` (seuil de
comptage), la période de référence, la fenêtre mensuelle.

### 0.4 Avertissement sur le niveau local

Avec **20 millésimes en dimension 6**, l'estimation de densité est en régime très défavorable.
Conséquences à assumer et à documenter :

- Le **top-1 parcellaire n'a aucune signification**. Ne jamais le présenter comme un résultat.
- Les estimateurs non paramétriques (KDE, kNN) sont ici des **contrôles**, pas des classements.
- Le predictif bayesien est le seul estimateur exploitant correctement 20 points en dimension 6,
  et il n'y parvient qu'avec un a priori de poids substantiel (`M0` de l'ordre de 15).
- Le diagnostic de **cohérence spatiale intra-région** (§8.5) est le juge de paix : si une parcelle
  bien classée n'a pas de voisines bien classées dans sa propre région, c'est du bruit.

---

## 1. Données

### 1.1 Structure attendue

| objet | dimensions | description |
|---|---|---|
| `D` | ~940 000 × 23 | une ligne = un couple (parcelle, millésime) |
| `yb` | 20 | vecteur nommé du profil Bordeaux 2025 |

Colonnes de `D` :

- `site` : identifiant de parcelle (~47 000 modalités)
- `region` : identifiant de région viticole (~700 modalités)
- `year` : millésime
- 20 colonnes climatiques : `{Tmin,Tmax,P,ET0}_{04..08}`

### 1.2 Alignement hémisphérique

Hémisphère nord : avril–août de l'année *N*.

Hémisphère sud : le millésime *N* couvre **octobre *N−1* à février *N***. Ces cinq mois doivent être
rangés **dans les mêmes colonnes** `_04` à `_08`, dans l'ordre chronologique de la saison végétative
(octobre → `_04`, novembre → `_05`, … février → `_08`).

> **Vérification critique.** Confirmer que l'étiquetage des millésimes austraux suit la convention de
> récolte. Une erreur ici décale tout l'hémisphère sud d'un an et invalide silencieusement la moitié
> du résultat.

### 1.3 Période de référence

Prendre la période de 20 ans **la plus récente disponible** (p. ex. 2005–2024). Le retrait de
tendance (§2.1) atténue le problème, mais ne remplace pas des données récentes.

### 1.4 Homogénéité

- Source climatique unique pour Bordeaux 2025 et pour les 47 000 parcelles.
- ET0 : formulation identique partout (Penman-Monteith FAO-56 de préférence). Un mélange
  Hargreaves / Penman suffit à dominer le classement.

---

## 2. Préalables

### 2.1 Retrait de tendance et recalage au niveau climatique 2025

Les 20 millésimes ne sont pas échangeables : ils contiennent un réchauffement. On réexprime chaque
millésime **au niveau climatique de 2025**, par site et par variable :

$$x^{adj}_{t} = x_t - \hat{b}\,(t - 2025)$$

Le nuage représente alors la distribution interannuelle actuelle, et « millésime typique » redevient
bien défini. Mahony et al. (2017) calculent de même leur écart-type de référence après retrait de la
tendance.

### 2.2 Scores normaux plutôt que `log1p`

Les cumuls de précipitations d'avril–août sont fortement asymétriques et bornés à zéro en climat
méditerranéen ou semi-aride. Plutôt que de choisir une transformation ad hoc, on applique une
**transformation en scores normaux** (copule gaussienne) : les marges deviennent normales par
construction et la structure de dépendance est conservée. Appliquée à **toutes** les variables pour
la cohérence, et à l'identique aux 940 000 lignes et au vecteur cible.

---

## 3. Code — préparation

```r
library(data.table)

# ---------------------------------------------------------------------------
# ENTRÉES
#   D  : data.table ~940 000 lignes — site, region, year + 20 colonnes climatiques
#   yb : vecteur nommé de 20 valeurs — profil Bordeaux 2025, mêmes noms que D
# ---------------------------------------------------------------------------

TARGET_YEAR <- 2025L
MOIS <- 4:8
VARS <- as.vector(outer(c("Tmin", "Tmax", "P", "ET0"),
                        sprintf("%02d", MOIS), paste, sep = "_"))

setDT(D)
stopifnot(all(VARS %in% names(D)), all(VARS %in% names(yb)))
yb <- yb[VARS]

# --- contrôles d'intégrité --------------------------------------------------
stopifnot(!anyNA(D[, ..VARS]), !anyNA(yb))
nn <- D[, .N, by = site]$N
if (length(unique(nn)) != 1L)
  warning("nombre de millésimes non constant : ", paste(range(nn), collapse = "-"))
message("sites : ", uniqueN(D$site), " | régions : ", uniqueN(D$region),
        " | lignes : ", format(nrow(D), big.mark = " "))

# --- 2.1 retrait de tendance, recalage à TARGET_YEAR ------------------------
# x_adj = x - b*(year - TARGET_YEAR) ; pas besoin de l'intercept
D[, (VARS) := {
    yy  <- year
    ww  <- yy - mean(yy)
    Sww <- sum(ww^2)
    lapply(.SD, function(x) x - (sum(ww * x) / Sww) * (yy - TARGET_YEAR))
  }, by = site, .SDcols = VARS]

# --- 2.2 scores normaux (copule gaussienne) ---------------------------------
# ORDRE IMPORTANT : transformer la cible AVANT d'écraser la colonne de D,
# sinon l'ECDF appliquée à yb n'est plus celle des valeurs d'origine.
EPS <- 1e-6
for (v in VARS) {
  Fv    <- ecdf(D[[v]])
  yb[v] <- qnorm(pmin(pmax(Fv(yb[v]), EPS), 1 - EPS))
  set(D, j = v, value = qnorm(pmin(pmax(Fv(D[[v]]), EPS), 1 - EPS)))
}
```

---

## 4. Espace ACP et métrique poolée

```r
# ---------------------------------------------------------------------------
#  METRIQUE DE FOND : distance de Mahalanobis sur covariance intra-site poolee
#    Mahalanobis (1936), Proc. Natl. Inst. Sci. India 2(1), 49-55
#                       reed. Sankhya A doi:10.1007/s13171-019-00164-5
#  APPLICATION AUX ANALOGUES CLIMATIQUES
#    Mahony et al. (2017)  doi:10.1111/gcb.13645  [sigma dissimilarity]
#    Grenier et al. (2013) doi:10.1175/JAMC-D-12-0170.1  [comparaison de
#                          six metriques de dissimilarite climatique]
# ---------------------------------------------------------------------------
X   <- as.matrix(D[, ..VARS])
pca <- prcomp(X, center = TRUE, scale. = TRUE)

K <- 6L
message("K = ", K, " CP | variance retenue = ",
        round(100 * sum(pca$sdev[1:K]^2) / sum(pca$sdev^2), 1), " %")

# --- SORTIE DIAGNOSTIQUE DE L'ACP -------------------------------------------
# Une seule ACP pour tout le jeu mondial : ces pourcentages sont GLOBAUX,
# ils ne decrivent aucun site en particulier (ni Bordeaux ni un autre).
vp  <- pca$sdev^2
acp <- data.table(CP       = paste0("PC", seq_along(vp)),
                  var_pct  = round(100 * vp / sum(vp), 2),
                  cum_pct  = round(100 * cumsum(vp) / sum(vp), 2))
print(head(acp, 12))
message("K = ", K, " -> variance cumulee retenue : ", acp$cum_pct[K], " %")

# eboulis + seuils de reference
plot(acp$var_pct[1:15], type = "b", pch = 19,
     xlab = "composante", ylab = "variance (%)", main = "Eboulis")
abline(v = K, col = "red", lty = 2)

# interpretation climatique des axes : correlations variables / composantes
charg <- pca$rotation[, 1:K] * rep(pca$sdev[1:K], each = nrow(pca$rotation))
charg <- round(charg, 2)
message("--- correlations variables d'origine / composantes ---")
print(charg)

# pour chaque axe, les 5 variables les plus structurantes
for (j in 1:K) {
  o <- order(abs(charg[, j]), decreasing = TRUE)[1:5]
  message("PC", j, " (", acp$var_pct[j], " %) : ",
          paste0(rownames(charg)[o], " ", charg[o, j], collapse = " | "))
}

Z  <- pca$x[, 1:K, drop = FALSE]
zb <- drop(predict(pca, matrix(yb, nrow = 1, dimnames = list(NULL, VARS))))[1:K]
rm(X); gc()

# --- covariance intra-site poolée -------------------------------------------
pcols <- paste0("PC", 1:K)
DT <- data.table(site = D$site, region = D$region, year = D$year)
DT[, (pcols) := as.data.table(Z)]

ANO <- copy(DT)
ANO[, (pcols) := lapply(.SD, function(x) x - mean(x)), by = site, .SDcols = pcols]

df_pool <- nrow(ANO) - uniqueN(ANO$site)
W <- crossprod(as.matrix(ANO[, ..pcols])) / df_pool
message("ddl poolés : ", format(df_pool, big.mark = " "),
        " | conditionnement de W : ", round(kappa(W), 2))
rm(ANO); gc()

Rw <- chol(W)                                    # W = Rw'Rw

# --- blanchiment : euclidien dans Zw == Mahalanobis dans Z ------------------
Zw  <- t(backsolve(Rw, t(Z),  transpose = TRUE))
zbw <- backsolve(Rw, zb, transpose = TRUE)

wcols <- paste0("W", 1:K)
DTw <- data.table(site = D$site, region = D$region, year = D$year)
DTw[, (wcols) := as.data.table(Zw)]
DTw[, d := sqrt(rowSums(sweep(Zw, 2, zbw, "-")^2))]
```

`d` est la distance de chaque millésime du monde à Bordeaux 2025, exprimée en **écarts de
millésime** : unité commune, comparable entre Bordeaux, Napa et Mendoza.

---

## 5. Les trois estimateurs

### 5.1 Estimateur principal — predictif bayesien (Normal-inverse-Wishart)

#### 5.1.0 Explication pour lecteur non bayesien

> **A lire avant le code.** Cette sous-section n'apporte aucune formule, elle explique l'idee. Elle
> est destinee a etre reutilisable telle quelle dans un article ou une presentation.
>
> **Le probleme.** On veut savoir quelle est la *forme* du nuage des millesimes d'un vignoble : est-il
> allonge dans le sens chaud-sec, resserre, incline ? Cette forme, c'est la matrice de covariance. Or
> on ne dispose que de 20 millesimes pour l'estimer en dimension 6, soit 21 parametres a partir de
> 20 observations. L'estimation classique est donc tres instable : changez un millesime, la forme
> bascule.
>
> **L'idee bayesienne.** Avant meme de regarder ce vignoble precis, on sait deja quelque chose. La
> variabilite interannuelle du climat a une allure generale, commune a l'ensemble des vignobles du
> monde : c'est exactement la covariance poolee `W`, estimee sur ~893 000 degres de liberte. C'est ce
> qu'on appelle **l'a priori** : une croyance formulee avant l'observation locale.
>
> On regarde ensuite les 20 millesimes du vignoble, et on **met a jour** cette croyance. Le resultat,
> appele **a posteriori**, est un compromis automatique :
>
> - si les 20 millesimes sont peu informatifs ou incoherents, le resultat reste proche de la forme
>   generale `W` ;
> - s'ils racontent clairement autre chose, la forme locale prend le dessus.
>
> **Le lien avec le lambda de la version precedente.** C'est le meme mecanisme de compromis, mais la
> nature du reglage change. Auparavant, `lambda` etait choisi par validation croisee : un chiffre qui
> optimise un critere, sans signification directe. Desormais, le reglage est un **poids d'a priori**
> note `nu_0`, qui se lit comme une phrase : « ma forme generale `W` pese autant que `nu_0`
> millesimes fictifs ». C'est une hypothese **declaree**, qu'un relecteur peut discuter, et non un
> parametre ajuste.
>
> **La consequence pratique.** Au lieu d'une densite gaussienne, on obtient une **loi de Student
> multivariee**. Concretement, ses queues sont plus epaisses : la Student accorde plus de
> vraisemblance aux valeurs eloignees du centre que la gaussienne. C'est une facon honnete de
> reconnaitre qu'avec 20 millesimes on ne connait pas la forme du nuage avec certitude. Un millesime
> extreme est donc moins brutalement rejete — ce qui compte precisement pour un millesime hors norme
> comme Bordeaux 2026.
>
> **Ce que cette approche ne fait pas.** Elle ne rend pas le resultat « plus vrai ». Elle rend
> l'hypothese de regularisation explicite et le comportement en petit echantillon plus prudent.
> Le classement obtenu est proche de celui de la gaussienne regularisee (recouvrement attendu
> ~95 %) ; le gain est conceptuel et de robustesse en queue, pas un bouleversement des resultats.

#### 5.1.0bis Statut epistemologique de la construction — a lire avant de rediger

> Cette sous-section repond a une objection previsible en relecture : *melanger un estimateur
> bayesien et deux estimateurs non bayesiens, est-ce coherent ?* La reponse est oui, mais a
> condition de nommer correctement ce qu'on fait. Trois points.

**1. Ce que nous faisons s'appelle un *empirical Bayes*, pas un bayesien complet.**

L'a priori `W` n'est pas une croyance externe : il est **estime sur les memes donnees** que celles
qu'on analyse (les ~893 000 degres de liberte de la covariance intra-site poolee). C'est la
definition de l'*empirical Bayes*, ou *Bayes empirique*, formalise comme maximum de vraisemblance de
type II dans un modele hierarchique.

- Efron B. (2010). *Large-Scale Inference: Empirical Bayes Methods for Estimation, Testing, and
  Prediction*. IMS Monographs 1, Cambridge University Press. doi:10.1017/CBO9780511761362
- Robbins H. (1956). An empirical Bayes approach to statistics. *Proc. Third Berkeley Symposium on
  Mathematical Statistics and Probability*, vol. 1, 157–163. — texte fondateur.

**Il faut ecrire « empirical Bayes » ou « predictif bayesien a a priori empirique », jamais
« analyse bayesienne » tout court.** C'est la difference entre une formulation defendable et une
formulation attaquable.

**2. La critique connue, a assumer explicitement.**

L'empirical Bayes « utilise les donnees deux fois » : une fois pour construire l'a priori, une fois
pour l'inference. Consequence : les incertitudes a posteriori sont **sous-estimees**, et les
garanties d'optimalite du cadre bayesien strict ne tiennent plus. C'est la critique classique de
Lindley, toujours d'actualite dans la litterature methodologique recente.

*Pourquoi c'est acceptable ici, et comment le dire.* Nous n'utilisons **aucune** quantite
d'incertitude a posteriori : ni intervalle de credibilite, ni probabilite, ni test. Le predictif ne
sert qu'a produire un **score de classement ordinal** entre sites. La sous-estimation de
l'incertitude n'affecte donc pas nos conclusions, qui reposent sur des rangs. En revanche, cela
**interdit** de presenter le score comme une probabilite, ou de dire « ce vignoble a X pour cent de
chances de ». Formulation a retenir dans l'article : *l'a priori empirique regularise l'estimation
de forme en petit echantillon ; le score en resultant est utilise comme statistique de classement,
non comme mesure probabiliste calibree.*

**3. Le melange de paradigmes est une pratique etablie, sous le nom de comparaison multi-estimateurs.**

Nos trois estimateurs ne sont pas combines en une seule valeur : ils sont **compares**, et leur
desaccord est le diagnostic. C'est exactement l'usage recommande en ecologie.

- Dormann C.F., Calabrese J.M., Guillera-Arroita G., Matechou E., Bahn V., Bartoń K., et al. (2018).
  Model averaging in ecology: a review of Bayesian, information-theoretic, and tactical approaches
  for predictive inference. *Ecological Monographs* 88(4), 485–504. doi:10.1002/ecm.1309
  — **la reference a citer sur ce point.** Elle montre que l'agregation multi-modeles est utile
  surtout quand la covariance entre modeles est faible, ce qui justifie precisement notre choix de
  garder le kNN (non parametrique, faiblement correle) plutot qu'un quatrieme estimateur gaussien
  redondant. Elle recommande aussi de rapporter la dispersion entre modeles plutot que la seule
  moyenne — c'est notre `rg_span` et notre drapeau `solidite` (§5.5).

> **Consequence pour §5.5.** Notre consensus par rang median n'est **pas** un *Bayesian model
> averaging* : il n'y a pas de poids de vraisemblance de modele. C'est une agregation ordinale
> robuste, a decrire comme telle. Ne pas employer le terme « model averaging » sans le qualifier.

**4. Absence de precedent direct — a assumer aussi.**

Recherche faite : **aucune publication n'applique un predictif bayesien conjugue NIW a
l'identification d'analogues climatiques.** Les analogues climatiques emploient des distances
(Mahalanobis, euclidienne, Zech-Aslan) ; le NIW conjugue est employe ailleurs, en genomique et en
statistique multivariee. La construction proposee ici est donc **une combinaison nouvelle de briques
chacune bien etablies**, et non une methode attestee telle quelle.

C'est defendable, et c'est meme un argument de nouveaute, a condition de le presenter honnetement :

| brique | statut dans la litterature | reference d'ancrage |
|---|---|---|
| distance de Mahalanobis sur variabilite interannuelle detrendee | etabli en climatologie des analogues | Mahony et al. 2017 |
| covariance poolee intra-groupe comme metrique commune | etabli (metrique intra-classe de l'ADL) | Mardia, Kent & Bibby 1979 |
| shrinkage de covariance en petit echantillon | etabli, tres utilise | Schäfer & Strimmer 2005 |
| a priori conjugue NIW, predictif Student | etabli en statistique bayesienne | Gelman et al. 2013, ch. 3 |
| a priori estime sur les donnees (empirical Bayes) | etabli, avec reserves connues | Efron 2010 ; Robbins 1956 |
| KDE en espace environnemental issu d'une ACP | etabli en biogeographie | Broennimann et al. 2012 |
| densite kNN sur series climatiques | etabli en hydrologie | Lall & Sharma 1996 |
| comparaison de plusieurs estimateurs plutot qu'un seul | recommande en ecologie | Dormann et al. 2018 |

**Phrase type pour la section methodes :** *« Nous combinons trois estimateurs de densite dont les
hypotheses de forme differentes permettent d'evaluer la robustesse du classement : un predictif
bayesien a a priori empirique conjugue (NIW), un estimateur a noyau, et une densite par plus proches
voisins. A notre connaissance, cette construction n'a pas ete appliquee a l'identification
d'analogues de millesime ; chacune de ses composantes est cependant etablie [references]. »*

#### 5.1.1 Formulation

Modele conjugue sur les millesimes du site $s$, dans l'espace ACP a $K$ dimensions :

$$x \mid \mu, \Sigma \sim \mathcal{N}_K(\mu, \Sigma), \qquad
(\mu, \Sigma) \sim \text{NIW}(\mu_0, \kappa_0, \nu_0, \Psi_0)$$

Hyperparametres, tous ancres sur des quantites deja calculees :

| hyperparametre | valeur | lecture |
|---|---|---|
| $\mu_0$ | $\bar{x}_s$ (moyenne du site) | a priori non informatif sur la position |
| $\kappa_0$ | 0,01 | position quasi libre : c'est la *forme* qu'on regularise, pas le centre |
| $\nu_0$ | $K + 2 + m_0$, avec $m_0 = 15$ | `W` pese comme 15 millesimes fictifs |
| $\Psi_0$ | $(\nu_0 - K - 1)\, W$ | cale l'esperance a priori de $\Sigma$ sur exactement $W$ |

Le calage de $\Psi_0$ garantit $\mathbb{E}[\Sigma] = W$ a priori : en l'absence de donnees locales,
on retombe exactement sur la metrique poolee.

La densite predictive a posteriori est une **Student multivariee** en forme fermee — aucune
simulation MCMC :

$$x_0 \mid \text{data} \sim t_{\nu_n - K + 1}\!\left(\mu_n,\ \frac{\kappa_n + 1}{\kappa_n(\nu_n - K + 1)}\Psi_n\right)$$

avec les mises a jour standard :

$$\kappa_n = \kappa_0 + n, \quad \nu_n = \nu_0 + n, \quad
\mu_n = \frac{\kappa_0\mu_0 + n\bar{x}}{\kappa_n}$$
$$\Psi_n = \Psi_0 + (n-1)\,\text{cov}(Z_s) + \frac{\kappa_0 n}{\kappa_n}(\bar{x}-\mu_0)(\bar{x}-\mu_0)'$$

> **Continuite avec la version precedente.** $\Psi_n / (\nu_n - K - 1)$ est une moyenne ponderee de
> `W` et de `cov(Z_s)` : c'est structurellement le meme shrinkage que
> $S_s = \lambda W + (1-\lambda)\text{cov}(Z_s)$, avec $\lambda$ implicite
> $\approx m_0 / (m_0 + n - 1) = 15/34 \approx 0{,}44$ sur la forme, mais assorti de queues Student.
> Avec $m_0 = 15$ et $n = 20$, a priori et donnees pesent d'un poids comparable.

#### 5.1.2 Code

```r
# ===========================================================================
# 5.1  ESTIMATEUR PRINCIPAL — predictif bayesien NIW / Student
#      Cout identique au gaussien : forme fermee, une decomposition de
#      Cholesky K x K par site. Aucun MCMC.
#
#  SOURCE DE LA METHODE
#    Raiffa & Schlaifer (1961)  doi:10.1002/9781118625125   [cadre conjugue]
#    Gelman et al. (2013), Bayesian Data Analysis 3e, ch. 3
#                               doi:10.1201/b16018          [predictif Student, MAJ NIW]
#    Schafer & Strimmer (2005)  doi:10.2202/1544-6115.1175  [shrinkage de covariance]
#  EXEMPLE D'APPLICATION EN ENVIRONNEMENT
#    Mahony et al. (2017), Global Change Biology
#                               doi:10.1111/gcb.13645       [Mahalanobis regularise,
#                                                            analogues climatiques]
# ===========================================================================

M0     <- 15        # poids d'a priori, en "millesimes fictifs"  <- A DECLARER
KAPPA0 <- 0.01      # a priori faible sur la position

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

sc_bay <- DT[, {
  M    <- as.matrix(.SD)
  n    <- nrow(M)
  xbar <- colMeans(M)

  kn <- KAPPA0 + n
  nn <- nu0    + n
  mn <- (KAPPA0 * xbar + n * xbar) / kn      # mu0 = xbar => mn = xbar
  d  <- xbar - xbar                          # nul par construction
  Pn <- Psi0 + (n - 1) * cov(M) + (KAPPA0 * n / kn) * tcrossprod(d)

  df    <- nn - K + 1
  Scale <- Pn * (kn + 1) / (kn * df)

  .(score_bay = ldstudent(zb, mn, Scale, df), df_post = df)
}, by = .(site, region), .SDcols = pcols]

message("degres de liberte a posteriori : ", sc_bay$df_post[1],
        "  (plus eleve = plus proche de la gaussienne)")
```

> **Note sur $\mu_0$.** Prendre $\mu_0 = \bar{x}_s$ rend le terme de position degenere (`d = 0`) :
> le predictif est alors centre sur la moyenne du site, comme dans la version gaussienne. C'est
> voulu — on ne regularise que la **forme**. Le code conserve le terme general pour permettre de
> tester un $\mu_0$ regional (centroide de la region) si vous voulez aussi regulariser la position.

#### 5.1.3 Sensibilite a rapporter

`M0` remplace `lambda` comme unique parametre libre de l'estimateur principal. Il doit etre
**declare, pas optimise**. Rapporter neanmoins la stabilite du classement :

```r
for (m in c(5, 15, 30, 60)) {
  # relancer sc_bay avec M0 = m, puis :
  # recouvrement du top-200 avec le run de reference (M0 = 15)
}
```

Attendu : recouvrement du top-200 > 0,85 entre `M0 = 5` et `M0 = 30`. Si le classement bascule,
c'est que les nuages locaux sont trop peu informatifs et que **seul** l'a priori parle — a signaler
explicitement.

**Comparaison de continuite obligatoire.** Faire tourner une fois l'ancienne gaussienne regularisee
(`lambda` par CV) et rapporter la correlation de Spearman avec `score_bay`. Attendu > 0,95. Un
ecart plus important signale une erreur d'implementation, pas une decouverte.

### 5.2 Contrôle 1 — noyau gaussien

Remplace le noyau boîte du comptage. Même logique, variance nettement plus faible, et l'arbitraire
du seuil disparaît au profit d'une fenêtre `h`.

```r
# ---------------------------------------------------------------------------
#  SOURCE DE LA METHODE  (noyau gaussien / estimation a noyau)
#    Rosenblatt (1956)  doi:10.1214/aoms/1177728190
#    Parzen (1962)      doi:10.1214/aoms/1177704472
#  EXEMPLES D'APPLICATION EN ENVIRONNEMENT
#    Broennimann et al. (2012), Global Ecol. Biogeogr.
#                       doi:10.1111/j.1466-8238.2011.00698.x
#                       [noyaux en espace climatique issu d'une ACP -- cas le
#                        plus proche du notre ; package ecospat]
#    Qiao et al. (2017), Global Ecol. Biogeogr.  doi:10.1111/geb.12492
#                       [GARDE-FOU : biais du KDE selon dimension et effectif]
# ---------------------------------------------------------------------------
h <- sqrt(qchisq(0.50, df = K)) / 2
sc_kde <- DTw[, .(score_kde = mean(exp(-0.5 * (d / h)^2))), by = .(site, region)]
```

### 5.3 Contrôle 2 — densité kNN adaptative

$$\hat f_s(x_0) \propto \frac{m}{n \cdot d_{(m)}^{K}}$$

Aucune hypothèse de forme, fenêtre adaptative. `m = 5` plutôt que 3 pour limiter la variance.

```r
# ---------------------------------------------------------------------------
#  SOURCE DE LA METHODE  (densite par k plus proches voisins)
#    Loftsgaarden & Quesenberry (1965)  doi:10.1214/aoms/1177700079
#  EXEMPLE D'APPLICATION EN ENVIRONNEMENT
#    Lall & Sharma (1996), Water Resources Research  doi:10.1029/95WR02966
#                       [kNN sur series climatiques multivariees]
# ---------------------------------------------------------------------------
m <- 5L
sc_knn <- DTw[, {
  dm <- sort(d)[m]
  .(score_knn = log(m) - log(.N) - K * log(dm), d_m = dm)
}, by = .(site, region)]
```

### 5.4 Statistique de communication

Imbattable pour dire « 7 millésimes sur 20 ». Descriptive, jamais utilisée pour classer.

```r
r <- sqrt(qchisq(0.50, df = K))          # écart médian entre deux millésimes
sc_cnt <- DTw[, .(n_proches = sum(d <= r), n_tot = .N,
                  d_min = min(d), an_proche = year[which.min(d)]),
              by = .(site, region)]
```

### 5.5 Synthèse — tableau de comparaison des trois méthodes

Le tableau de sortie principal a **quatre colonnes de rang** : une par méthode, puis le consensus.
Les trois premières servent de contrôle visuel de la quatrieme.

Consensus par **rang median** des trois methodes, plus robuste que la moyenne : une methode qui
decroche seule sur un site ne fait pas basculer la ligne. L'**etendue** des trois rangs (`rg_span`)
mesure la fragilite.

```r
res <- Reduce(function(a, b) merge(a, b, by = c("site", "region")),
              list(sc_bay, sc_kde, sc_knn, sc_cnt))

# --- un rang par methode (1 = meilleur) -------------------------------------
res[, `:=`(rg_bay = frank(-score_bay, ties.method = "min"),
           rg_kde = frank(-score_kde, ties.method = "min"),
           rg_knn = frank(-score_knn, ties.method = "min"))]

# --- consensus : rang median + fragilite ------------------------------------
res[, `:=`(rg_med  = apply(cbind(rg_bay, rg_kde, rg_knn), 1, median),
           rg_span = apply(cbind(rg_bay, rg_kde, rg_knn), 1,
                           function(x) diff(range(x))))]
res[, rg_cons := frank(rg_med, ties.method = "min")]

# --- drapeau de solidite ----------------------------------------------------
N1 <- ceiling(0.01 * nrow(res))
res[, accord := (rg_bay <= N1) + (rg_kde <= N1) + (rg_knn <= N1)]
res[, solidite := fifelse(accord == 3L, "3/3",
                  fifelse(accord == 2L, "2/3",
                  fifelse(accord == 1L, "1/3", "0/3")))]

# --- TABLEAU PRINCIPAL : 4 colonnes de rang + consensus ---------------------
setorder(res, rg_cons)
TAB <- res[, .(site, region,
               bayesien = rg_bay, kde = rg_kde, knn = rg_knn,
               consensus = rg_cons, span = rg_span, solidite,
               n_proches, n_tot)]
print(head(TAB, 30))
```

### 5.6 Les trois top-N cote a cote

```r
top_par_methode <- function(dat, N = 10L, unite = "region") {
  f <- function(col) {
    d <- dat[order(dat[[col]])]
    if (unite == "region") d <- d[!duplicated(region)]
    d[1:N, get(unite)]
  }
  data.table(rang      = 1:N,
             bayesien  = f("rg_bay"),
             kde       = f("rg_kde"),
             knn       = f("rg_knn"),
             consensus = f("rg_cons"))[]
}

cote_a_cote <- top_par_methode(res, N = 10L, unite = "region")
print(cote_a_cote)

freq <- sort(table(unlist(cote_a_cote[, -1])), decreasing = TRUE)
message("--- presence dans les top-10 (max 4) ---")
print(freq)
```

Lecture : une region presente dans les quatre colonnes est un analogue solide. Une region presente
dans une seule est un artefact d'estimateur, a ne pas mettre en avant.

## 6. Profil mensuel — le mois qui disqualifie

Une parcelle peut afficher 8 millésimes proches en agrégé tout en divergeant systématiquement en
avril, ce qui annule son intérêt sous l'angle du risque gel.

```r
profil_mensuel <- function(mo) {
  v   <- paste0(c("Tmin", "Tmax", "P", "ET0"), "_", sprintf("%02d", mo))
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

Dm <- vapply(MOIS, profil_mensuel, numeric(nrow(D)))
colnames(Dm) <- month.abb[MOIS]

prof <- data.table(site = D$site, region = D$region, Dm)
prof_site <- prof[, lapply(.SD, median), by = .(site, region),
                  .SDcols = month.abb[MOIS]]
prof_site[, mois_pire := month.abb[MOIS][apply(.SD, 1, which.max)],
          .SDcols = month.abb[MOIS]]

res <- merge(res, prof_site, by = c("site", "region"))
```

---

## 7. Sorties de la phase 1

1. `res` — tableau parcellaire complet : 3 rangs, `n_proches`, profil mensuel, `mois_pire`.
2. Distribution des rangs **par région** (boxplot), pas le top-1 parcellaire.
3. Carte des parcelles colorées par `rg_bay`, pour lire la cohérence spatiale à l'œil.
4. Profils mensuels des 10 meilleures régions au sens de la médiane des rangs, superposés à Bordeaux.
5. Tableau de recouvrement des top-N entre estimateurs et entre variantes de paramètres.

---

## 8. Diagnostics obligatoires

### 8.1 Contrôle interne — Bordeaux contre lui-même

Bordeaux figure parmi les 47 000 parcelles. Son `n_proches` dit combien de millésimes bordelais de
la période ressemblent à 2025. **Ce nombre calibre tous les autres.**

```r
res[grepl("bordeaux", region, ignore.case = TRUE)][order(rg_bay)][1:20]
res[grepl("bordeaux", region, ignore.case = TRUE), summary(n_proches)]
```

- Médiane ≈ 2 → 2025 est singulier chez lui ; une région étrangère à 8 est une trouvaille solide.
- Médiane ≈ 12 → 2025 n'a rien de singulier ; le cadrage de l'article doit changer.

### 8.2 Bootstrap sur millésimes — le diagnostic central au niveau local

Avec n = 20, c'est le test qui dit si le classement est reproductible ou s'il tient au bruit
d'échantillonnage.

```r
boot_top <- function(B = 30L, N = 200L) {
  ref <- res[order(rg_bay)][1:N, site]
  ov  <- numeric(B)
  for (b in seq_len(B)) {
    idx <- DT[, .I[sample(.N, .N, replace = TRUE)], by = site]$V1
    Db  <- DT[idx]
    sb  <- Db[, {
      M  <- as.matrix(.SD)
      n  <- nrow(M)
      kn <- KAPPA0 + n; nn <- nu0 + n
      Pn <- Psi0 + (n - 1) * cov(M)
      df <- nn - K + 1
      S  <- Pn * (kn + 1) / (kn * df)   # matrice d'echelle Student
      R  <- chol(S)
      x <- backsolve(R, zb - colMeans(M), transpose = TRUE)
      .(sc = -sum(log(diag(R))) - 0.5 * sum(x^2))
    }, by = site, .SDcols = pcols]
    ov[b] <- length(intersect(ref, sb[order(-sc)][1:N, site])) / N
  }
  ov
}
ov <- boot_top()
message("recouvrement bootstrap du top-200 : ",
        round(100 * mean(ov)), " % [", round(100 * min(ov)), "–",
        round(100 * max(ov)), "]")
```

Sous 50 % de recouvrement, ne publier aucun classement parcellaire : passer directement à la
phase 2 régionale.

### 8.3 Concordance entre estimateurs

```r
topN <- function(col, N = 200L) res[order(res[[col]])][1:N, site]
for (p in list(c("rg_bay","rg_knn"), c("rg_bay","rg_kde"), c("rg_kde","rg_knn")))
  message(p[1], " / ", p[2], " : ",
          length(intersect(topN(p[1]), topN(p[2]))), "/200")
```

Un écart fort entre le predictif bayesien et le kNN signale que l'hypothèse de forme mord quelque part : vérifier
l'asymétrie résiduelle et la multimodalité des nuages concernés.

### 8.4 Sensibilité aux paramètres

Rejouer avec `K ∈ {4, 6, 8}`, `M0 ∈ {5, 15, 30, 60}`, `m ∈ {3, 5, 8}`,
`r = sqrt(qchisq(q, K))` pour `q ∈ {0,35 ; 0,50 ; 0,65}`. Rapporter le recouvrement des top-200,
dans l'esprit du critère de Grenier et al. (2013), qui montrent que les métriques concordent à
grande échelle mais divergent sur l'identité des meilleurs analogues. Conclusion à publier : une
**région d'analogues**, pas une parcelle unique.

### 8.5 Cohérence spatiale intra-région — le juge de paix

Les ~67 parcelles d'une même région sont climatiquement voisines. Si une parcelle est bien classée
et ses consœurs ne le sont pas, c'est du bruit d'estimation, pas un signal.

```r
coh <- res[, .(n = .N,
               rg_med = median(rg_bay),
               rg_iqr = IQR(rg_bay),
               part_top5pct = mean(rg_bay <= 0.05 * nrow(res))),
           by = region][order(rg_med)]
head(coh, 25)

# une région crédible : rg_med bas ET rg_iqr modéré ET part_top5pct élevée
plot(coh$rg_med, coh$rg_iqr, log = "xy",
     xlab = "rang médian", ylab = "IQR des rangs")
```

### 8.6 Calibration gaussienne

```r
qqplot(qchisq(ppoints(1e5), df = K), sample(rowSums(Zw^2), 1e5),
       xlab = "chi2 théorique", ylab = "distances au carré observées")
abline(0, 1, col = "red")
```

Si la queue décroche, présenter les résultats en rangs et percentiles plutôt qu'en probabilités.

---

## 9. Phase 2 — passage au niveau régional

À faire **après** validation de la phase 1, et seulement si §8.2 ou §8.5 montrent que le niveau
parcellaire est trop bruité pour être publié tel quel.

Le changement est conceptuel autant que technique : on estime la densité sur la **région** comme
unité, soit ~67 parcelles × 20 millésimes ≈ 1 340 points. Les parcelles étant fortement corrélées,
l'échantillon effectif est bien inférieur à 1 340, mais clairement supérieur à 20. On passe d'un
régime désespéré à un régime praticable pour l'estimation non paramétrique en dimension 6.

L'objet estimé devient aussi plus juste : « l'enveloppe climatique de la région, variabilité spatiale
interne comprise ». C'est ce qu'on veut dire par « à Napa, on trouve souvent des millésimes comme
Bordeaux 2025 ».

Modifications à prévoir :

- `by = .(site, region)` → `by = region` dans les §5.1 à 5.4.
- `M0` peut descendre nettement (l'a priori pèse moins face à un `n` régional élevé) : `cov` régionale sur ~1 340 points est stable.
- KDE et kNN redeviennent des estimateurs de plein droit, plus seulement des contrôles.
- `W` reste la métrique intra-site poolée : ne pas la recalculer en intra-région, sous peine
  d'absorber la variabilité spatiale dans l'unité de mesure.
- Conserver le niveau parcellaire pour la carte de variabilité intra-région, qui reste un résultat
  en soi si les corrections topographiques sont fines.

---

## 10. Notes de performance

- 940 000 × 20 en double ≈ 150 Mo. `prcomp` passe sans difficulté. Libérer `X` et `ANO` après usage.
- Le retrait de tendance (§3) fait 47 000 groupes × 20 colonnes : compter quelques dizaines de
  secondes en `data.table`.
- Le score bayesien fait 47 000 `chol` sur des matrices 6 × 6 : de l'ordre de la minute.
- Le bootstrap §8.2 est le poste le plus coûteux : le limiter à B = 30 et l'exécuter une fois.
- Si contrainte mémoire forte, estimer l'ACP sur un échantillon de 200 000 lignes puis projeter le
  reste avec `predict`.

---

## 11. Écarté après examen

- **Prédiction conforme.** Donnerait une p-valeur de typicalité valide en échantillon fini sous
  échangeabilité, mais ne répond qu'à l'estimand « typicalité relative », avec une résolution de
  1/21. À citer, pas à substituer.
- **DTW sur le profil mensuel.** Le décalage temporel n'a pas de sens dans une saison à calendrier
  fixe ; l'alignement hémisphérique est déjà traité par permutation de colonnes.
- **Analyse de données fonctionnelles.** Sur 5 pas de temps, l'ACP fonctionnelle dégénère vers l'ACP
  multivariée.
- **Classification une-classe (SVM, forêt aléatoire).** 20 observations par classe, aucune
  interprétabilité agronomique, aucun avantage sur le predictif bayesien.

---

## 12. Références

### 12.1 Origine des trois estimateurs — sources primaires

**Toutes les references ci-dessous ont ete verifiees par recherche web : le DOI resout vers la
notice de l'article. Aucune n'est reconstituee de memoire.**

**Estimateur 1 — predictif bayesien conjugue (Normal-inverse-Wishart / Student)**

- *Origine du cadre conjugue.* Raiffa H., Schlaifer R. (1961). *Applied Statistical Decision
  Theory*. Harvard Business School. — reedition Wiley Classics, doi:10.1002/9781118625125
- *Formulation moderne du predictif Student et des mises a jour NIW.* Gelman A., Carlin J.B.,
  Stern H.S., Dunson D.B., Vehtari A., Rubin D.B. (2013). *Bayesian Data Analysis*, 3e ed.,
  chapitre 3. Chapman & Hall/CRC. doi:10.1201/b16018
- *A priori estime sur les donnees — empirical Bayes.* Robbins H. (1956). An empirical Bayes
  approach to statistics. *Proc. Third Berkeley Symposium on Mathematical Statistics and
  Probability* 1, 157–163. ; Efron B. (2010). *Large-Scale Inference*. IMS Monographs 1, Cambridge
  University Press. doi:10.1017/CBO9780511761362
  — **le nom exact de notre construction ; voir §5.1.0bis pour la critique associee.**
- *Regularisation de covariance en petit echantillon, base du choix de `W` comme a priori.*
  Schäfer J., Strimmer K. (2005). A shrinkage approach to large-scale covariance matrix estimation
  and implications for functional genomics. *Statistical Applications in Genetics and Molecular
  Biology* 4(1), art. 32. doi:10.2202/1544-6115.1175

**Estimateur 2 — noyau gaussien (KDE)**

- Rosenblatt M. (1956). Remarks on some nonparametric estimates of a density function.
  *The Annals of Mathematical Statistics* 27(3), 832–837. doi:10.1214/aoms/1177728190
- Parzen E. (1962). On estimation of a probability density function and mode.
  *The Annals of Mathematical Statistics* 33(3), 1065–1076. doi:10.1214/aoms/1177704472

**Estimateur 3 — densite kNN adaptative**

- Loftsgaarden D.O., Quesenberry C.P. (1965). A nonparametric estimate of a multivariate density
  function. *The Annals of Mathematical Statistics* 36(3), 1049–1051. doi:10.1214/aoms/1177700079

**Metrique de fond — distance de Mahalanobis**

- Mahalanobis P.C. (1936). On the generalised distance in statistics. *Proceedings of the National
  Institute of Sciences of India* 2(1), 49–55. — reedition annotee : *Sankhya A* 80 (2018),
  doi:10.1007/s13171-019-00164-5

### 12.2 Applications en sciences de l'environnement

**KDE en espace climatique / biogeographie**

- Broennimann O., Fitzpatrick M.C., Pearman P.B., Petitpierre B., Pellissier L., Yoccoz N.G.,
  Thuiller W., Fortin M.-J., Randin C., Zimmermann N.E., Graham C.H., Guisan A. (2012). Measuring
  ecological niche overlap from occurrence and spatial environmental data. *Global Ecology and
  Biogeography* 21(4), 481–497. doi:10.1111/j.1466-8238.2011.00698.x
  — **le precedent le plus proche de notre usage** : noyaux appliques a des densites en espace
  environnemental issu d'une ACP, exactement notre construction. Implemente dans le package
  `ecospat` (fonction `ecospat.grid.clim.dyn`).
- Qiao H., Escobar L.E., Saupe E.E., Ji L., Soberón J. (2017). A cautionary note on the use of
  hypervolume kernel density estimators in ecological niche modelling. *Global Ecology and
  Biogeography* 26(9), 1066–1070. doi:10.1111/geb.12492
  — **a citer comme garde-fou** : montre que le KDE sur- ou sous-estime les volumes selon la
  dimensionnalite et le nombre d'observations, et conclut qu'il n'est viable qu'avec de grands
  effectifs et peu de dimensions. Justifie directement notre §0.4 et le statut de *controle*
  donne au KDE.

**Comparaison de plusieurs estimateurs plutot qu'un seul**

- Dormann C.F., Calabrese J.M., Guillera-Arroita G., Matechou E., Bahn V., Bartoń K., et al. (2018).
  Model averaging in ecology: a review of Bayesian, information-theoretic, and tactical approaches
  for predictive inference. *Ecological Monographs* 88(4), 485–504. doi:10.1002/ecm.1309
  — justifie de conserver des estimateurs faiblement correles entre eux, et de rapporter leur
  dispersion plutot que la seule moyenne (notre `rg_span`, §5.5).

**Metrique intra-classe**

- Mardia K.V., Kent J.T., Bibby J.M. (1979). *Multivariate Analysis*. Academic Press. — covariance
  intra-groupe poolee comme metrique de reference (fondement de l'analyse discriminante lineaire).

**kNN en climatologie / hydrologie**

- Lall U., Sharma A. (1996). A nearest neighbor bootstrap for resampling hydrologic time series.
  *Water Resources Research* 32(3), 679–693. doi:10.1029/95WR02966
  — reference fondatrice de l'usage du kNN sur series climatiques multivariees ; l'estimation de
  densite kNN y sert de base au reechantillonnage, sur de gros corpus.

**Shrinkage / Mahalanobis regularise en climatologie des analogues**

- Mahony C.R., Cannon A.J., Wang T., Aitken S.N. (2017). A closer look at novel climates: new
  methods and insights at continental to landscape scales. *Global Change Biology* 23(9),
  3934–3955. doi:10.1111/gcb.13645
  — *sigma dissimilarity* : Mahalanobis sur variabilite interannuelle detrendee, puis conversion en
  quantiles par la loi du chi. Notre §2.1 et notre metrique poolee en derivent.

### 12.3 Analogues climatiques — contexte general

- Grenier P., Parent A.-C., Huard D., Anctil F., Chaumont D. (2013). An assessment of six
  dissimilarity metrics for climate analogs. *J. Appl. Meteorol. Climatol.* 52(4), 733–752.
  doi:10.1175/JAMC-D-12-0170.1
- Mahony C.R., Cannon A.J., Wang T., Aitken S.N. (2017). A closer look at novel climates.
  *Global Change Biology* 23. doi:10.1111/gcb.13645
- Williams J.W., Jackson S.T., Kutzbach J.E. (2007). Projected distributions of novel and
  disappearing climates by 2100 AD. *PNAS* 104, 5738–5742.
- Fitzpatrick M.C., Dunn R.R. (2019). Contemporary climatic analogs for 540 North American urban
  areas. *Nature Communications* 10. doi:10.1038/s41467-019-08540-3
- Hallegatte S., Hourcade J.-C., Ambrosi P. (2007). Using climate analogues for assessing climate
  change economic impacts in urban areas. *Climatic Change* 82, 47–60.
- Ramírez-Villegas J. et al. (2011). Climate analogues: finding tomorrow's agriculture today.
  CCAFS Working Paper 12.
- Allaman H., Goyette S., Dubuis P.-H., Kasparian J. (2025). Future viability of European vineyards
  using bioclimatic climate analogues. *Agric. For. Meteorol.* doi:10.1016/j.agrformet.2025.110978
- Crombie J., Brown L., Lizzio J., Hood G. (2008). *Climatch user manual*. Australian Bureau of
  Rural Sciences.
- Hubbard J.A.G., Drake D.A.R., Mandrak N.E. (2025). Euclimatch: an R package for climate matching
  with Euclidean distance metrics. *Ecography*. doi:10.1111/ecog.07614

### Packages R

`data.table` (obligatoire), `FNN` (si rejeu multi-cibles), `corpcor` (shrinkage alternatif vers une
cible diagonale), `Euclimatch` (CRAN, critère de la pire variable, en comparaison de robustesse),
`ClimaRep` (CRAN).
