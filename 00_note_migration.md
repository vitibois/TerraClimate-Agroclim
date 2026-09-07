# Note de migration — a lire en premier

**Destinataire : agent de codage (Claude Code) charge de mettre a jour un projet R existant.**
**Auteur du projet : Benjamin Bois — analogues climatiques de millesime.**
**Date : septembre 2026.**

---

## 0. Comment utiliser ce jeu de trois documents

Vous recevez trois fichiers markdown et un projet R existant.

| fichier | role | statut |
|---|---|---|
| `00_note_migration.md` | **ce document** — ce qui change par rapport au code existant | a lire en premier |
| `analogues_millesime.md` | specification de la methode, phase 1 niveau parcellaire | etat cible de reference |
| `note_station_delta_carte.md` | ancrage station, calibration, agregation regionale, cartographie, rapport | etat cible de reference |

**Ordre de lecture imperatif :** ce document, puis `analogues_millesime.md`, puis
`note_station_delta_carte.md`.

**Regle centrale.** Les deux documents de specification decrivent **l'etat cible**. Ils ne disent pas
ce qui change. Ce document-ci dit ce qui change. Deux erreurs a eviter absolument :

1. **Ne pas repartir de zero.** Le projet existant contient du code fonctionnel : chargement des
   donnees, alignement hemispherique, ACP, boucles de scores. Il doit etre **modifie**, pas reecrit.
2. **Ne pas laisser cohabiter l'ancienne et la nouvelle logique.** Plusieurs elements de l'ancienne
   version doivent etre **supprimes**, pas conserves en parallele. Ils sont listes en §2.

---

## 1. Historique du projet en une page

### La question de recherche

A partir d'une base mondiale (`WLD`) d'environ 47 000 points de vigne repartis dans environ
700 regions viticoles, avec des donnees climatiques mensuelles sur 20 ans (environ 940 000 lignes),
determiner **dans quels vignobles du monde on trouve frequemment des millesimes au profil climatique
tres proche d'un millesime cible**. Autrement dit : de quel vignoble ce millesime serait-il une annee
*ordinaire* ?

Variables : Tmin, Tmax, P, ET0, sur avril a aout, soit 20 colonnes.

Cible de travail : **Bordeaux 2026**, attendu comme un « ovni climatique » devant trouver ses
analogues en Mediterranee chaude et seche. Les premiers resultats ont confirme cette attente
qualitativement. La methode doit ensuite se generaliser a n'importe quel millesime de n'importe
quelle region.

### Ce que ce n'est pas

**Ce n'est pas une classification.** Pas de classification ascendante hierarchique, pas de
partition, pas de clustering. L'estimand est un **score de plausibilite site par site** : la densite
de la distribution des millesimes du site, evaluee au point cible.

### Les etapes deja validees en session

- Synthese de litterature (climatologie de l'adaptation, ecologie, biosecurite, agronomie).
- Architecture statistique arretee : trois estimateurs de densite concurrents, compares et non
  fusionnes.
- Retrait de tendance obligatoire, covariance intra-site poolee comme metrique commune, troncature
  ACP comme regularisation.
- Passage de l'estimateur gaussien regularise a un **predictif bayesien a a priori empirique**.
- Verification par recherche web de l'integralite des references bibliographiques.
- Clarification du **statut epistemologique** de la construction (empirical Bayes, et non bayesien
  strict) — point critique pour la publication.

---

## 2. Ce qui doit etre SUPPRIME du code existant

**Ces elements ne doivent plus figurer nulle part apres migration.** Les laisser en commentaire est
acceptable ; les laisser actifs est un bug.

### 2.1 L'estimateur gaussien regularise comme estimateur principal

- Supprimer la fonction de validation croisee `loglik_lovo()` et la recherche de `lambda` sur grille.
- Supprimer la variable `lambda` et tous ses usages en production.
- **Exception unique :** conserver la possibilite de faire tourner l'ancien estimateur **une seule
  fois**, dans un script de controle separe, pour le test de continuite (§3.1). Il ne fait plus
  partie du pipeline.

### 2.2 L'ancrage de la region cible par son nom

L'ancienne version identifiait la region d'ancrage par correspondance de chaine sur le nom de region.
**Supprime.** L'ancrage se fait desormais par coordonnees geographiques (§3.3).

### 2.3 Tout calcul de centroide regional

- Supprimer tout `mean(lon)`, `mean(lat)` par region.
- Supprimer toute moyenne de variables climatiques par region utilisee comme representant de la
  region.
- Supprimer tout affichage cartographique de centroides.

**Justification, posee par Benjamin :** la diversite climatique interne d'une region est une
propriete reelle, pas du bruit. Un analogue est **un lieu**, pas une moyenne. Un centroide regional
peut tomber dans une zone climatiquement inexistante.

### 2.4 Le classement presente comme une liste de points

L'ancienne sortie etait un classement de points individuels. **Remplace** par la logique de cohortes
regionales sequentielles (§3.4).

### 2.5 Les transformations `log1p` sur les precipitations

Remplacees par des scores normaux via copule gaussienne, appliques a **toutes** les variables.

---

## 3. Ce qui doit etre AJOUTE ou REMPLACE

### 3.1 Estimateur principal : predictif bayesien NIW / Student

**Reference : `analogues_millesime.md` §5.1.**

Remplace l'estimateur gaussien regularise. Cout de calcul identique : forme fermee, une
decomposition de Cholesky K par K par site, **aucun MCMC**.

Renommages a propager dans l'ensemble du projet :

| ancien nom | nouveau nom |
|---|---|
| `sc_gau` | `sc_bay` |
| `score_gau` | `score_bay` |
| `rg_gau` | `rg_bay` |
| colonne de sortie `gaussien` | `bayesien` |

Nouveau parametre : `M0 = 15`, poids d'a priori exprime en « millesimes fictifs ». Il est
**declare, pas optimise**. Ne pas ecrire de validation croisee pour le choisir.

**Test de continuite obligatoire.** Faire tourner l'ancien estimateur gaussien une fois et rapporter
la correlation de Spearman avec `score_bay`. **Attendu superieur a 0,95.** En dessous, c'est une
erreur d'implementation, pas une decouverte. Signaler et s'arreter.

**Test de sensibilite.** `M0` dans {5, 15, 30, 60}, recouvrement du top-200 attendu superieur a 0,85.

### 3.2 Section pedagogique et epistemologique

**Reference : `analogues_millesime.md` §5.1.0 et §5.1.0bis.**

Ces deux sous-sections sont du texte, pas du code. Elles ne demandent aucune implementation, mais
**elles ne doivent pas etre perdues** : elles sont destinees a etre reprises dans un article.

Point critique a respecter dans tout commentaire de code, tout nom de variable et toute sortie
textuelle : la methode s'appelle **empirical Bayes**, ou **predictif bayesien a a priori empirique**.
**Jamais « analyse bayesienne » tout court.** L'a priori etant estime sur les memes donnees, les
incertitudes a posteriori sont sous-estimees ; le score doit donc etre presente comme une
**statistique de classement ordinal**, jamais comme une probabilite.

### 3.3 Ancrage par coordonnees geographiques

**Reference : `note_station_delta_carte.md` §B.**

La station meteo est reperee par `LOC_LON` et `LOC_LAT`. On l'apparie au point de `WLD` le
**plus proche geographiquement**, par distance de haversine, avec un seuil `D_MAX_KM = 15`.

**Precision importante :** `WLD` contient des points avec leurs coordonnees propres, a une densite
d'environ un point par kilometre carre. Ce n'est **pas** un raster a extraire. C'est un rapprochement
point a point.

Prevoir un garde-fou contre l'inversion longitude / latitude.

### 3.4 Cohortes regionales sequentielles

**Reference : `note_station_delta_carte.md` §G.2.**

Logique posee par Benjamin. On parcourt les points par rang croissant. Chaque nouvelle region
rencontree prend le rang regional suivant, et sa cohorte est constituee de ses points contigus dans
le classement.

Exemple : rangs 1 a 50 tous dans la region X, donc X est la region numero 1 avec une cohorte de 50 ;
rangs 51 a 53 dans la region Y, donc Y est la region numero 2, rang de fin 53 ; et ainsi de suite.

Colonnes a produire : `rang_deb`, `rang_fin`, `n_cohorte`, `n_region`, `pct_region`.

### 3.5 Tableau de comparaison a quatre colonnes de rang

**Reference : `analogues_millesime.md` §5.5 et §5.6.**

Colonnes `bayesien`, `kde`, `knn`, `consensus`. Consensus par **rang median** des trois. Ajouter
`rg_span`, l'etendue des trois rangs, et un drapeau `solidite` en trois niveaux : 3 sur 3, 2 sur 3,
1 sur 3.

**Ne pas appeler cela du « model averaging ».** Il n'y a pas de poids de vraisemblance de modele.
C'est une agregation ordinale robuste.

### 3.6 Calibration station contre grille

**Reference : `note_station_delta_carte.md` §C et §D.**

Delta calcule sur les annees communes, contre le **point apparie** et jamais contre une moyenne
regionale. **Additif** pour Tmin et Tmax, **multiplicatif** pour P et ET0, avec un plancher
`P_MIN = 5 mm`. Drapeau `fiabilite` en quatre niveaux. Conserver a la fois `yb_brut` et `yb_corr`.
La correction va toujours **vers** l'echelle de `WLD`.

### 3.7 Cartographie, plan ACP, rapport PDF

**Reference : `note_station_delta_carte.md` §G.5 a §G.10.**

Points reels uniquement, aucun centroide. Palette unique `PAL` partagee par la carte, le plan ACP,
les profils mensuels et le tableau. Rapport PDF assemble dans
`sorties/rapport_analogues.pdf`.

### 3.8 Bloc diagnostique ACP

**Reference : `analogues_millesime.md` §4.**

Variance par composante et variance cumulee, eboulis avec seuil `K` materialise, correlations
variables / composantes, cinq variables les plus structurantes par axe. **Avertissement a
conserver :** ces pourcentages sont globaux, pas locaux.

### 3.9 Commentaires bibliographiques dans le code

**Format demande par Benjamin, a respecter :** source de la methode, exemple d'application, DOI,
journal. Sans detail superflu. Deja redige aux §4, §5.1, §5.2 et §5.3 de
`analogues_millesime.md` — reprendre tel quel.

---

## 4. Ce qui NE CHANGE PAS

A conserver en l'etat, sans y toucher :

- Le chargement des donnees et la structure `data.table`.
- **L'alignement hemispherique** : pour l'hemisphere sud, millesime N = octobre N-1 a fevrier N,
  range dans les colonnes `_04` a `_08` par ordre chronologique. Une erreur ici decale silencieusement
  la moitie du corpus.
- Le retrait de tendance : `x_adj = x - b * (year - TARGET_YEAR)`.
- Le calcul de la covariance intra-site poolee `W`.
- La troncature ACP a `K = 6`.
- Les deux estimateurs de controle : noyau gaussien et densite kNN adaptative avec `m = 5`.
- La statistique `n_proches` avec `r = sqrt(qchisq(0.5, K))`, reservee a la communication.
- Le blanchiment par Cholesky, avec `chol` puis `backsolve(..., transpose = TRUE)`.

---

## 5. Pieges techniques deja identifies

- **Ordre dans la boucle de scores normaux :** transformer `yb` **avant** d'ecraser `D[[v]]`.
- **Non-independance des points :** la taille d'echantillon effective est d'environ 700 regions, pas
  47 000 points. Ne jamais raisonner sur 47 000 comme un effectif independant.
- **Homogeneite de ET0 :** verifier une formulation identique partout, Penman-Monteith FAO-56.
- **Noms de colonnes de `WLD` a verifier** avant execution : `CNT`, `WINE_REGION`, `site`, `lon`,
  `lat`.

---

## 6. Avertissement de fond a ne pas diluer

**Reference : `analogues_millesime.md` §0.4.**

Vingt millesimes en dimension six est un regime statistique defavorable. **Le top-1 parcellaire n'a
aucune signification.** Les diagnostics de robustesse ne sont pas optionnels :

- Bootstrap sur millesimes (§8.2) : **sous 50 % de recouvrement du top-200, ne rien publier au niveau
  parcellaire.**
- Coherence spatiale intra-region (§8.5) : c'est le juge de paix.

Ces deux diagnostics conditionnent le passage en phase 2 regionale.

---

## 7. Ordre de travail suggere

1. Lire les trois documents.
2. Inventorier le code existant et localiser chaque element de la §2 a supprimer.
3. Implementer le predictif bayesien (§3.1) et faire tourner le test de continuite (§3.1) **avant
   toute autre modification**. Si Spearman est inferieur a 0,95, s'arreter et signaler.
4. Basculer l'ancrage sur les coordonnees (§3.3).
5. Implementer les cohortes regionales (§3.4) et le tableau a quatre colonnes (§3.5).
6. Ajouter les diagnostics ACP (§3.8) et les commentaires bibliographiques (§3.9).
7. Cartographie, plan ACP et rapport PDF (§3.7).
8. Faire tourner les diagnostics obligatoires du §8 de `analogues_millesime.md`.

**Signaler tout point ou la specification est ambigue plutot que de trancher seul.**
