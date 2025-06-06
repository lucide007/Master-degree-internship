# Code Refactored
```r
# Load necessary packages
library(gridExtra)
library(grid)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(stringr)
library(dryR)
library(rain)
library(Metacycle)
library(UpSetR)
library(GGally)

# Load silico datasets
silico_data_1 <- read.csv("path/to/CIS_data_2025-04-25.csv")
silico_data_2 <- read.csv("path/to/CIS_data_2025-04-17.csv")

# Preprocess the data
silico_data <- bind_rows(silico_data_1, silico_data_2) %>%
  select(-c(Cosine.Amplitude, Cosine.Periods, Cosine.Lag.Factors, X)) %>%
  na.omit() %>%
  rename_with(~ sub("^X", "", .), everything())

# Create replicates
temps_seq <- seq(1.5, 24, by = 1.5)
rep1 <- silico_data[, 1:16] %>% rename_with(~ paste0("Rep1_T", temps_seq), everything())
rep2 <- silico_data[, 17:32] %>% rename_with(~ paste0("Rep2_T", temps_seq), everything())
silico_data <- cbind(rep1, rep2)

# DryR analysis
silico_round <- round(silico_data) + 7
drylist_silico <- dryseq_single_multiperiod(silico_round, rep("silico", 32), time_minute, period_range = c(630, 810), period_step = c(6))

# Rain analysis
rain_silico <- rain(t(silico_data), deltat = 1.5, period = 12, period.delta = 1.5, adjp.method = "BH", peak.border = c(0.3, 0.7))
rain_silico_results <- rownames(rain_silico)[rain_silico$pVal < 0.05]

# Metacycle analysis
time_minute <- rep(seq(1.5, 24, by = 1.5) * 60, each = 2)
write.csv(silico_data, file = "silico_data_for_metacycle.csv", row.names = FALSE)
meta2d_silico <- meta2d(
  infile = "silico_data_for_metacycle.csv", filestyle = "csv",
  outdir = "metacycle_liver_results", timepoints = time_minute,
  minper = 720, maxper = 720, combinePvalue = "bonferroni",
  cycMethod = c("JTK", "LS", "ARS"), ARSdefaultPer = 7# =============================================================================
# Script R pour l'analyse de rythmes circadiens avec données in silico
# Comparaison de méthodes de détection de rythmes : RAIN, JTK, ARSER, Lomb-Scargle
# =============================================================================

# =============================================================================
# 1. CHARGEMENT DES PACKAGES
# =============================================================================

library(gridExtra)    # Arrangement de graphiques multiples
library(grid)         # Fonctions graphiques de base
library(ggplot2)      # Visualisation de données
library(patchwork)    # Composition de graphiques
library(dplyr)        # Manipulation de données
library(tidyr)        # Restructuration de données
library(stringr)      # Manipulation de chaînes de caractères
library(dryR)         # Analyse de rythmes circadiens
library(rain)         # Algorithme RAIN pour rythmes
library(Metacycle)    # Meta-analyse de méthodes de détection de rythmes
library(UpSetR)       # Visualisation d'intersections d'ensembles
library(GGally)       # Extension ggplot2 pour matrices de graphiques
library(data.table)   # Pour fwrite()

# =============================================================================
# 2. CHARGEMENT ET PRÉPARATION DES DONNÉES IN SILICO
# =============================================================================

cat("Chargement des données in silico...\n")

# Chargement de deux jeux de données CIS (Circadian In Silico)
# 24h, 1 point / 1.5h, 2 replicat, 30% de features rythmiques avec periode 12h min amp 1 max amp 6 out amp 0
silico_data_1 <- read.csv("path/to/CIS_data_2025-04-25.csv")
# 48h, 1 point / 1.5h, 1 replicat, 30% features rythmiques avec periode 10.5-13.5h min amp 1 max amp 3 out amp 0
silico_data_2 <- read.csv("path/to/CIS_data_2025-04-17.csv")

silico_data2 <- subset(silico_data, select = -c(Cosine.Amplitude, Cosine.Periods, Cosine.Lag.Factors))
silico_data3 <- subset(silico_data2, select = -X)
silico_data3 <- na.omit(silico_data3)

# =============================================================================
# 3. CRÉATION DES RÉPLICATS BIOLOGIQUES
# =============================================================================
# Pour le jeu de données CIS_data_2025-04-17.csv
cat("Création des réplicats biologiques...\n")

# Définition de la séquence temporelle (1.5h à 24h par pas de 1.5h)
temps_seq <- seq(1.5, 24, by = 1.5)

# Division en deux réplicats
# Réplicat 1 : colonnes 1-16 (premier cycle de 24h)
rep1 <- silico_data3[, 1:16] %>% 
  rename_with(~ paste0("Rep1_T", temps_seq), everything())

# Réplicat 2 : colonnes 17-32 (second cycle de 24h)  
rep2 <- silico_data3[, 17:32] %>% 
  rename_with(~ paste0("Rep2_T", temps_seq), everything())

# Fusion des réplicats
silico_data3 <- cbind(rep1, rep2)

# =============================================================================
# 4. ANALYSE AVEC DryR (DÉTECTION DE RYTHMES)
# =============================================================================

cat("Analyse DryR pour détection de rythmes...\n")

# Définition des timepoints en minutes pour DryR
time_minute <- rep(seq(1.5, 24, by = 1.5) * 60, each = 2)

# Transformation des données : arrondi + offset pour éviter les zéros
silico_round <- round(silico_data3) + 7

# Analyse DryR avec période de 12h (630-810 minutes)
drylist_silico <- dryseq_single_multiperiod(
  data = silico_round,
  group = rep("silico", 32),  # Tous les échantillons dans le même groupe
  time = time_minute,
  period_range = c(630, 810), # 10.5h à 13.5h (centré sur 12h)
  period_step = c(6)          # Pas de 6 minutes
)


# =============================================================================
# 5. ANALYSE AVEC RAIN (RHYTHMICITY ANALYSIS)
# =============================================================================

cat("Analyse RAIN...\n")

# Analyse RAIN avec période de 12h
rain_silico <- rain(
  t(silico_data3),           # Transposition : gènes en lignes
  deltat = 1.5,             # Intervalle de temps (1.5h)
  period = 12,              # Période d'intérêt (12h)
  period.delta = 0,       # Intervalle de période recherchée
  adjp.method = "BH",       # Correction de Benjamini-Hochberg
  peak.border = c(0.3, 0.7) # Largeur du pic acceptable
)

# Analyse RAIN avec période de 12+/-1.5h
rain_silico <- rain(
  t(silico_data3),           # Transposition : gènes en lignes
  deltat = 1.5,             # Intervalle de temps (1.5h)
  period = 12,              # Période d'intérêt (12h)
  period.delta = 1.5,       # Intervalle de période recherchée
  adjp.method = "BH",       # Correction de Benjamini-Hochberg
  peak.border = c(0.3, 0.7) # Largeur du pic acceptable
)

# Analyse RAIN avec période de 12+/-1.5h et 2 replicas par ZT
rain_silico <- rain(
  t(silico_data3),           # Transposition : gènes en lignes
  deltat = 1.5,             # Intervalle de temps (1.5h)
  period = 12,              # Période d'intérêt (12h)
  period.delta = 1.5,       # Intervalle de période recherchée
  adjp.method = "BH",       # Correction de Benjamini-Hochberg
  peak.border = c(0.3, 0.7), # Largeur du pic acceptable
  nr.series = 2             # Nombre de replicas par ZT
)

# Extraction des gènes rythmiques significatifs (p < 0.05)
rain_silico_results <- rownames(rain_silico)[rain_silico$pVal < 0.05]

# Extraction des gènes rythmiques significatifs (p < 0.01)
rain_silico_results <- rownames(rain_silico)[rain_silico$pVal < 0.01]

cat("RAIN : ", length(rain_silico_results), " gènes rythmiques détectés\n")

# Extraction des gènes rythmiques dans les raw data
matrix_silico_rain <- silico_data3[rownames(silico_data3) %in% rain_silico_results,]


# =============================================================================
# 6. ANALYSE AVEC METACYCLE (META-ANALYSE DE MÉTHODES)
# =============================================================================

cat("Analyse Metacycle (JTK, Lomb-Scargle, ARSER)...\n")

# Création du vecteur temps
time <- c(1.5,3,4.5,6,7.5,9,10.5,12,13.5,15,16.5,18,19.5,21,22.5,24)
time <- rep(time, each = 2)
time_minute <- time * 60

# Sauvegarde des données pour Metacycle
write.csv(silico_data2, 
          file = "silico_data_for_metacycle.csv", 
          row.names = FALSE)

# Analyse meta2d combinant plusieurs méthodes pour une période de 12h
meta2d_silico <- meta2d(
  infile = "silico_data_for_metacycle.csv",
  filestyle = "csv",
  outdir = "metacycle_liver_results_fix",
  timepoints = time_minute,
  minper = 720,                          # Période minimale : 12h
  maxper = 720,                          # Période maximale : 12h
  combinePvalue = "bonferroni",          # Correction de Bonferroni
  cycMethod = c("JTK", "LS", "ARS"),     # Méthodes : JTK-Cycle, Lomb-Scargle, ARSER
  ARSdefaultPer = 720                    # Période par défaut pour ARSER
)

# Analyse meta2d combinant plusieurs méthodes pour une période de 12+/-1.5h
meta2d_silico <- meta2d(
  infile = "silico_data_for_metacycle.csv",
  filestyle = "csv",
  outdir = "metacycle_liver_results_intervalle",
  timepoints = time_minute,
  minper = 630,                          # Période minimale : 10.5h
  maxper = 810,                          # Période maximale : 13.5h
  combinePvalue = "bonferroni",          # Correction de Bonferroni
  cycMethod = c("JTK", "LS", "ARS"),     # Méthodes : JTK-Cycle, Lomb-Scargle, ARSER
  ARSdefaultPer = 720                    # Période par défaut pour ARSER
)
               

# =============================================================================
# 7. EXTRACTION ET COMPARAISON DES RÉSULTATS
# =============================================================================

cat("Extraction des résultats Metacycle...\n")

# Chargement des résultats de meta2d
meta2d_silico_results <- read.table("metacycle_liver_results_fix/meta2d_silico_data_for_metacycle.csv", header = T, sep = ",")
# Ou
meta2d_silico_results <- read.table("metacycle_liver_results_intervalle/meta2d_silico_data_for_metacycle.csv", header = T, sep = ",")

# Création des ensembles de gènes rythmiques par méthode
gene_silico_sets <- list(
  "Methode RAIN" = rain_silico_results,
  "Methode Lomb-Scargle" = which(meta2d_silico_results$LS_pvalue < 0.05),
  "Methode JTK" = which(meta2d_silico_results$JTK_pvalue < 0.05),
  "Methode ARSER" = which(meta2d_silico_results$ARS_pvalue < 0.05),
  "Methode Metacycle 2d" = which(meta2d_silico_results$meta2d_pvalue < 0.05)
)


# =============================================================================
# 8. CRÉATION DE LA MATRICE D'INTERSECTION
# =============================================================================

cat("\nCréation de la matrice d'intersection...\n")

# Création de la matrice binaire d'appartenance
matrix_12_silico <- unique(unlist(gene_silico_sets)) %>% as.data.frame()
colnames(matrix_12_silico) <- "gene"
matrix_12_silico <- matrix_12_silico %>%
  mutate(rain = as.integer(gene %in% gene_silico_sets$`Methode RAIN`),
         ls = as.integer(gene %in% gene_silico_sets$`Methode Lomb-Scargle`),
         jtk = as.integer(gene %in% gene_silico_sets$`Methode JTK`),
         ars = as.integer(gene %in% gene_silico_sets$`Methode ARSER`),
         meta2d = as.integer(gene %in% gene_silico_sets$`Methode Metacycle 2d`))

matrix_12.1 <- matrix_12_silico[rowSums(matrix_12_silico[, c("rain", "ls", "jtk", "ars", "meta2d")]) >= 2,]
matrix_12.2 <- matrix_12_silico[matrix_12_silico$rain == 1,]
matrix_12_silico_selected <- unique(rbind(matrix_12.1, matrix_12.2))

# Sauvegarde de la matrice de gènes rythmiques
fwrite(matrix_12_silico_selected, file="matrix_silico_12h_rhythmic_genes2.csv", sep=";", row.names = TRUE, col.names = TRUE)


# =============================================================================
# 9. VISUALISATION DES INTERSECTIONS (UPSET PLOT)
# =============================================================================

cat("Génération du graphique UpSet...\n")



# Création du graphique UpSet
upsetplot_silico <- upset(
  matrix_12_silico, 
  sets = c("rain", "ls", "jtk", "ars", "meta2d"),  
  order.by = "freq",  # Tri par fréquence
  mainbar.y.label = "12h rhythmic genes in methods intersections",
  sets.x.label = "Set sizes (in silico)",
  point.size = 5,  
  line.size = 1.2,
  matrix.dot.alpha = 0.8,
  text.scale = 2)
upsetplot_silico

# =============================================================================
# 10. ANALYSE DES CORRÉLATIONS DE PHASE
# =============================================================================

# Création des matrices de phase
matrix_phase_silico_jtk <- subset(matrix_silico_jtk, select = c(CycID, JTK_adjphase))
matrix_phase_silico_ars <- subset(matrix_silico_ars, select = c(CycID, ARS_adjphase))
matrix_phase_silico_ls <- subset(matrix_silico_ls, select = c(CycID, LS_adjphase))
matrix_phase_silico_meta2d <- subset(matrix_silico_meta2d, select = c(CycID, meta2d_phase))
matrix_phase_silico_rain <- subset(rain_silico, pVal < 0.01, select = phase)
matrix_phase_silico_rain$CycID <- rownames(matrix_phase_silico_rain)
colnames(matrix_phase_silico_rain)[colnames(matrix_phase_silico_rain) == "phase"] <- "rain_phase"

matrix_phase_silico <- full_join(matrix_phase_silico_rain[, c("CycID", "rain_phase")],
                                 matrix_phase_silico_jtk[, c("CycID", "JTK_adjphase")], by = "CycID") %>%
  full_join(matrix_phase_silico_ars[, c("CycID", "ARS_adjphase")], by = "CycID") %>%
  full_join(matrix_phase_silico_ls[, c("CycID", "LS_adjphase")], by = "CycID") %>%
  full_join(matrix_phase_silico_meta2d[, c("CycID", "meta2d_phase")], by = "CycID")

matrix_phase_silico_ars$CycID <- as.character(matrix_phase_silico_ars$CycID)
matrix_phase_silico_rain$CycID <- as.character(matrix_phase_silico_rain$CycID)
matrix_phase_silico_jtk$CycID <- as.character(matrix_phase_silico_jtk$CycID)
matrix_phase_silico_ls$CycID <- as.character(matrix_phase_silico_ls$CycID)
matrix_phase_silico_meta2d$CycID <- as.character(matrix_phase_silico_meta2d$CycID)

cat("Analyse des corrélations de phase...\n")

# =============================================================================
# 11. VISUALISATION DES CORRÉLATIONS DE PHASE
# =============================================================================

cat("Génération du graphique de corrélations de phase...\n")

# Création du graphique de corrélations par paires
phase_correlation_plot <- ggpairs(
  matrix_phase_silico[, -1],  # Exclusion de la colonne CycID
  title = "Phases correlations : silico, 12h rhythmic genes",
  upper = list(continuous = wrap("cor", size = 6)),              # Corrélations en haut
  lower = list(continuous = wrap("points", alpha = 0.8, size = 1.5, color = "royalblue")), # Nuages de points en bas
  diag = list(continuous = wrap("densityDiag", size = 0.8))      # Densités sur la diagonale
) +
  theme_bw() +
  theme(axis.text = element_text(size = 12))

phase_correlation_plot

# =============================================================================
# 12. NETTOYAGE ET SAUVEGARDE
# =============================================================================

cat("Nettoyage de la mémoire...\n")

# Suppression des objets temporaires pour libérer la mémoire
rm(list = ls(pattern = "silico"))

cat("Analyse terminée avec succès !\n")

# =============================================================================
# RÉSUMÉ DE L'ANALYSE
# =============================================================================
# 1. Chargement et fusion de données in silico CIS
# 2. Création de réplicats biologiques (2 cycles de 24h)
# 3. Détection de rythmes avec 5 méthodes :
#    - DryR : analyse multi-période
#    - RAIN : détection de rythmes robuste
#    - JTK-Cycle : test de Jonckheere-Terpstra
#    - Lomb-Scargle : analyse spectrale
#    - ARSER : détection de rythmes par régression
# 4. Comparaison des ensembles de gènes rythmiques
# 5. Visualisation des intersections (UpSet plot)
# 6. Analyse des corrélations de phase entre méthodes
# 
# Fichiers de sortie :
# - silico_data_for_metacycle.csv : données pour Metacycle
# - matrix_silico_12h_rhythmic_genes2.csv : matrice d'intersection
# - metacycle_liver_results/ : résultats détaillés Metacycle
# =============================================================================
