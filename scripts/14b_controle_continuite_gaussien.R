# 2026-09-04 B. Bois (with Claude Code)
#
# CONTROLE DE CONTINUITE -- ancien estimateur gaussien regularise.
#
# 00_note_migration.md section 3.1 : l'estimateur gaussien regularise
# (shrinkage lambda choisi par validation croisee leave-one-vintage-out) ne
# fait PLUS PARTIE DU PIPELINE. Il est conserve ici, dans un script separe, a
# la seule fin du test de continuite avec le nouvel estimateur principal
# (predictif bayesien a a priori empirique, NIW / Student).
#
#   A EXECUTER UNE SEULE FOIS, apres la migration.
#   Ne jamais l'appeler depuis 14_climate_analogues_millesime.R.
#
# CRITERE : correlation de Spearman entre score_gau et score_bay > 0,95.
# En dessous, c'est une ERREUR D'IMPLEMENTATION, pas une decouverte :
# signaler et s'arreter, ne pas interpreter l'ecart.
#
# ENTREE : le .rds ecrit par la section 15 du script principal
#          (DT_pc, W, zb, K, pcols, M0, score_bay).

library(data.table)

source("R/config.R")

## ---- 0. entree -------------------------------------------------------------

site_label   <- "Pessac Haut Lafue Nord_20000101_20260901"   # cf. script principal
station_name <- sub("_[0-9]{8}(_[0-9]{8})?$", "", site_label)
dir_out      <- file.path(data_root, "Analogues_Millesime", station_name)
target_year <- 2026L
in_path <- file.path(dir_out, paste0("AnaloguesMillesime_", site_label, "_",
                                     target_year, "_continuite_inputs.rds"))

if (!file.exists(in_path))
  stop("entree absente : ", in_path,
       "\nLancer d'abord 14_climate_analogues_millesime.R avec save_continuite <- TRUE")

IN <- readRDS(in_path)
DT <- IN$DT_pc; W <- IN$W; zb <- IN$zb; K <- IN$K; pcols <- IN$pcols
message("entree lue : ", nrow(DT), " lignes | ", uniqueN(DT$site), " sites | K = ", K,
        " | M0 du run bayesien = ", IN$M0)

n_ech_cv <- 400L    # sites echantillonnes pour la CV de lambda
set.seed(1)

## ---- 1. ANCIEN estimateur -- gaussienne regularisee, lambda par CV --------
# Conserve ici a titre historique uniquement. S <- lambda*W + (1-lambda)*cov.

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

ech  <- sample(unique(DT$site), min(n_ech_cv, uniqueN(DT$site)))
grid <- seq(0.1, 1, by = 0.1)
LL <- vapply(ech, function(s) {
  M <- as.matrix(DT[site == s, ..pcols])
  vapply(grid, function(l) loglik_lovo(M, l), numeric(1))
}, numeric(length(grid)))
lambda <- grid[which.max(rowSums(LL))]
message("lambda choisi par validation croisee leave-one-vintage-out : ", lambda)

pdf_path <- file.path(dir_out, paste0("AnaloguesMillesime_", site_label, "_",
                                      target_year, "_continuite.pdf"))
pdf(pdf_path)

plot(grid, rowSums(LL), type = "b", pch = 19, xlab = "lambda",
     ylab = "log-vraisemblance CV", main = "Choix de lambda (ancien estimateur)")

sc_gau <- DT[, {
  M <- as.matrix(.SD)
  S <- lambda * W + (1 - lambda) * cov(M)
  R <- chol(S)
  x <- backsolve(R, zb - colMeans(M), transpose = TRUE)
  .(score_gau = -sum(log(diag(R))) - 0.5 * sum(x^2) - 0.5 * K * log(2 * pi))
}, by = site, .SDcols = pcols]
sc_gau[, rg_gau := frank(-score_gau, ties.method = "min")]

## ---- 2. test de continuite -------------------------------------------------

CMP <- merge(sc_gau, IN$score_bay, by = "site")
rho <- cor(CMP$score_gau, CMP$score_bay, method = "spearman")
ov200 <- length(intersect(CMP[order(rg_gau)][1:200, site],
                          CMP[order(rg_bay)][1:200, site])) / 200

message("\n===== TEST DE CONTINUITE =====")
message("Spearman(score_gau, score_bay) = ", round(rho, 4), "   [attendu > 0,95]")
message("recouvrement des top-200        = ", round(100 * ov200), " %")

plot(CMP$rg_gau, CMP$rg_bay, pch = 16, cex = .2, col = "#00000033", log = "xy",
     xlab = "rang -- gaussienne regularisee (ancien)",
     ylab = "rang -- predictif bayesien (nouveau)",
     main = paste0("Continuite des classements, Spearman = ", round(rho, 3)))
abline(0, 1, col = "red")

dev.off()
message("wrote ", pdf_path)

fwrite(CMP, file.path(dir_out, paste0("AnaloguesMillesime_", site_label, "_",
                                      target_year, "_continuite.csv")))

if (rho <= 0.95)
  stop("CONTINUITE NON VERIFIEE (Spearman = ", round(rho, 4), " <= 0,95). ",
       "C'est une erreur d'implementation du predictif bayesien, pas une decouverte. ",
       "Verifier : Psi0 = (nu0 - K - 1) * W, nu0 = K + 2 + M0, la matrice d'echelle ",
       "Pn * (kn + 1) / (kn * df), et le fait que mu0 = xbar annule le terme de position.")

message("continuite verifiee : le nouvel estimateur reproduit l'ancien classement. ",
        "L'estimateur gaussien regularise peut etre definitivement abandonne.")
