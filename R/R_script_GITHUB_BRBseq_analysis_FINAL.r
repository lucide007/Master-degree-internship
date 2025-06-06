# =============================================================================
# Script R pour l'analyse de rythmes ultradiens à partir de séries temporelles RNAseq
# =============================================================================

# for more accuracy, we use separate dataset for each tissue (liver, soleus, gastrocnemius), with 2 replicates / ZT (0 to 22.5h)
# we work on log2 normalized counts with DESeq2::vst()
# we always set the random seed for reproductibility (set.seed(123)) before run a random calcul
# functions written for this analysis are in R_script_GITHUB_BRBseq_analysis_functions
# for the heatmaps, we use a dataset with mean gene expression of replicates at each ZT

setwd("path/to/repertory")

# =============================================================================
# CHARGEMENT DES PACKAGES ####
# =============================================================================

library(data.table)
library(Matrix)
library(readxl)
library(dplyr)
library(tidyverse)
library(gridExtra)
library(ggplot2)
library(viridis)
library(VennDiagram)
library(uwot)
library(pheatmap)
library(DESeq2)
library(rain)
library(UpSetR)
library(GGally)
library(patchwork)
library(cowplot)
library(enrichR)
# + fonction plot_enrichment()
# + fonction calculer_moyennes_ZT()
# + fonction circular_phase12H_histogram()
# + fonction plot_enrichment()

# =============================================================================
# CHARGEMENT DES DATASET ET OBJETS ####
# =============================================================================

# Datasets

# Foie
matrix_log_liver_rep <- read.delim("path/to/matrix_log_liver_nodedup_with_replicates.txt")
rownames(matrix_log_liver_rep) <-matrix_log_liver_rep$X
matrix_log_liver_rep <- subset(matrix_log_liver_rep, select = -X)

# Soleus
matrix_log_soleus_rep <- read.delim("path/to/matrix_log_soleus_nodedup_with_replicates.txt")
rownames(matrix_log_soleus_rep) <-matrix_log_soleus_rep$X
matrix_log_soleus_rep <- subset(matrix_log_soleus_rep, select = -X)

# Gastrocnemius
matrix_log_gastrocnemius_rep <- read.delim("path/to/matrix_log_gastrocnemius_nodedup_with_replicates.txt")
rownames(matrix_log_gastrocnemius_rep) <-matrix_log_gastrocnemius_rep$X
matrix_log_gastrocnemius_rep <- subset(matrix_log_gastrocnemius_rep, select = -X)

# Ajouter des prefixes aux noms de colonnes
add_prefix <- function(data, prefix) {
  colnames(data) <- paste0(prefix, colnames(data))
  return(data)
}

liver <- add_prefix(matrix_log_liver_rep, "liver_")
soleus <- add_prefix(matrix_log_soleus_rep, "soleus_")
gastrocnemius <- add_prefix(matrix_log_gastrocnemius_rep, "gastrocnemius_")

# Fusionner des dataframes
merged_df <- cbind(liver, soleus, gastrocnemius)

# Créer la matrice globale avec replicas, normalisée
prepare_log_matrix <- function(matrix_log, prefix) {
  matrix_log <- add_prefix(matrix_log, prefix)
  return(matrix_log)
}

matrix_log_rep <- cbind(
  prepare_log_matrix(matrix_log_liver_rep, "liver_"),
  prepare_log_matrix(matrix_log_soleus_rep, "soleus_"),
  prepare_log_matrix(matrix_log_gastrocnemius_rep, "gastrocnemius_")
)

# Supprimer les variables temporaires
rm(liver, soleus, gastrocnemius)

# =============================================================================
# REDUCTION DE DIMENSIONS SUR LE DATAFRAME GLOBAL ####
# =============================================================================

# Supprimer les lignes avec des NAs pour la PCA
matrix_log_rep <- na.omit(matrix_log_rep)

# Faire la PCA sur le dataframe centré sur le z-score
set.seed(123)
acp_all <- prcomp(
    t(matrix_log_rep),
    center = TRUE,
    scale. = TRUE
    )

# Afficher le résumé de la PCA
summary(acp_all)
# Conserver le nombre de composantes pour avoir 80% de variance au moins (si <10, en garder 10)

# Extraire les résultats de la PCA et ajouter la colonne avec noms des echantillons
pc_data <- as.data.frame(acp_all$x)
pc_data$ZT <- colnames(matrix_log_rep)

# Visualisation
ggplot(pc_data, aes(x = PC1, y = PC2, color = ZT, label = ZT)) +
  geom_text(size = 5) +
  labs(title = "PCA on all samples, log2-normalized counts",
       x = "Principal component 1 (54.37%)",
       y = "Principal component 2 (29.08%)") +
  theme_minimal()+
  theme(legend.position = "none")

# =============================================================================
# EXTRACTION DES RESULTATS PUBLIES DE ZHU 2020 ET REANALYSE ####
# =============================================================================

# data from table S3, used for figure S3 in "Regulation of mammalian 12-h clock by XBP1s" Zhu et al., 2020
# quick resume : RNAseq time serie in liver, every 2h during 48h in constant darkness condition on xbp1 flox male mice
# the team used RAIN algorithm with FDR cutoff set to 0.05 (a correction added manually, not the one in RAIN package)

# Rythmes ultradiens
rain_zhu <- read_excel("plos_zhu_2020_rain_genes.xlsx", sheet = "Sheet 1")
list_zhu_12 <- readLines("list_zhu_rain_12h.txt")
rain_zhu_12 <- rain_zhu[rain_zhu$...1 %in% list_zhu_12,]
rain_zhu_12 <- rain_zhu_12 %>%
  group_by(`...1`) %>%                  # groupe par le nom du gène si doublons
  slice_min(order_by = adjP, n = 1, with_ties = FALSE) %>% # garde la ligne avec la plus petite adjpVal
  ungroup()
rain_zhu_12 <- as.data.frame(rain_zhu_12)
rownames(rain_zhu_12) <- rain_zhu_12$...1
rain_zhu_12 <- subset(rain_zhu_12, select = -...1)
rain_zhu_12_sig <- rain_zhu_12[!grepl("rik", rownames(rain_zhu_12), ignore.case = TRUE), ]

# Rythmes circadiens
rain_zhu_24 <- read_excel("plos_zhu_2020_rain_genes.xlsx", sheet = "Sheet 3")
list_zhu_24 <- readLines("list_zhu_rain_24h.txt")
rain_zhu_24 <- rain_zhu_24[rain_zhu_24$...1 %in% list_zhu_24,]
rain_zhu_24 <- rain_zhu_24 %>%
  group_by(`...1`) %>%                  # groupe par le nom du gène si doublons
  slice_min(order_by = pVal, n = 1, with_ties = FALSE) %>% # garde la ligne avec la plus petite adjpVal
  ungroup()
rain_zhu_24 <- as.data.frame(rain_zhu_24)
rownames(rain_zhu_24) <- rain_zhu_24$...1
rain_zhu_24 <- subset(rain_zhu_24, select = -...1)
rain_zhu_24_sig <- rain_zhu_24[!grepl("rik", rownames(rain_zhu_24), ignore.case = TRUE), ]

# Visualisation de l'intersection entre les dataframes
zhu_list <- list(
    "Zhu 12h" = rownames(rain_zhu_12_sig),
    "Zhu 24h" = rownames(rain_zhu_24_sig)
    )

venn_zhu <- venn.diagram(x = zhu_list,
  filename = NULL,  # Ne pas sauvegarder en fichier
  col = "transparent",
  fill = viridis::viridis(length(zhu_list), option = "C"),
  alpha = 0.5,
  label.col = "black",
  cex = 2,
  fontfamily = "sans",
  cat.cex = 2,
  cat.fontfamily = "sans",
  cat.pos = c(-20, 20),  # Positions personnalisées pour les noms des ensembles
  cat.dist = 0.05)  # Distance des noms par rapport aux cercles
grid::grid.newpage()
grid::grid.draw(venn_zhu)

# Rechercher un gène en particulier dans le dataframe
any(grepl("Spry4", rownames(rain_zhu_12_sig)))

# Extraire les gènes détectés comme uniquement ultradiens (exclusion des gènes à l'intersection)
list_zhu_12f <- rownames(rain_zhu_12_sig[!rownames(rain_zhu_12_sig) %in% rownames(rain_zhu_24_sig),])

# Analyse d'enrichissement avec enrichR
result_enrichr_zhu12 <- plot_enrichment(list_zhu_12, dbs, num_terms = 10)

result_enrichr_zhu12f <- plot_enrichment(list_zhu_12f, dbs, num_terms = 10)

# -> les données publiées correspondent à l'enrichissement sur les gènes avec intersection inclue

# =============================================================================
# ANALYSE DES SERIES TEMPORELLES AVEC RAIN ####
# =============================================================================

# Foie - RAIN ####

# Gènes ultradiens
rain_liver_rep_12 <- rain(
    t(matrix_log_liver_rep),
    deltat = 1.5,
    period = 12,
    period.delta = 1.5,
    adjp.method = "BH",
    peak.border = c(0.3,0.7),
    nr.series = 2
    )

# Extraire les gènes rythmiques
rain_liver_rep_sig_12 <- rain_liver_rep_12[rain_liver_rep_12$pVal < 0.05,]

# Gènes circadiens
rain_liver_rep_24 <- rain(
    t(matrix_log_liver_rep),
    deltat = 1.5,
    period = 24,
    period.delta = 2,
    adjp.method = "BH",
    peak.border = c(0.3,0.7),
    nr.series = 2
    )

# Extraire les gènes rythmiques
rain_liver_rep_sig_24 <- rain_liver_rep_24[rain_liver_rep_24$pVal < 0.05,]

# Soleus - RAIN ####

# Gènes ultradiens
rain_soleus_rep_12 <- rain(
    t(matrix_log_soleus_rep),
    deltat = 1.5,
    period = 12,
    period.delta = 1.5,
    adjp.method = "BH",
    peak.border = c(0.3,0.7),
    nr.series = 2
    )

# Extraire les gènes rythmiques
rain_soleus_rep_sig_12 <- rain_soleus_rep_12[rain_soleus_rep_12$pVal < 0.05,]

# Gènes circadiens
rain_soleus_rep_24 <- rain(
    t(matrix_log_soleus_rep),
    deltat = 1.5,
    period = 24,
    period.delta = 2,
    adjp.method = "BH",
    peak.border = c(0.3,0.7),
    nr.series = 2
    )

# Extraire les gènes rythmiques
rain_soleus_rep_sig_24 <- rain_soleus_rep_24[rain_soleus_rep_24$pVal < 0.05,]

# Gastrocnemius - RAIN ####

# Gènes ultradiens
rain_gastrocnemius_rep_12 <- rain(
    t(matrix_log_gastrocnemius_rep),
    deltat = 1.5,
    period = 12,
    period.delta = 1.5,
    adjp.method = "BH",
    peak.border = c(0.3,0.7),
    nr.series = 2
    )

# Extraire les gènes rythmiques
rain_gastrocnemius_rep_sig_12 <- rain_gastrocnemius_rep_12[rain_gastrocnemius_rep_12$pVal < 0.05,]

# Gènes circadiens
rain_gastrocnemius_rep_24 <- rain(
    t(matrix_log_gastrocnemius_rep),
    deltat = 1.5,
    period = 24,
    period.delta = 2,
    adjp.method = "BH",
    peak.border = c(0.3,0.7),
    nr.series = 2
    )

# Extraire les gènes rythmiques
rain_gastrocnemius_rep_sig_24 <- rain_gastrocnemius_rep_24[rain_gastrocnemius_rep_24$pVal < 0.05,]

# Sauvegarder les résultats
dir.create("RAIN_results")

write.table(rain_liver_rep_12, file = "RAIN_results/rain_results_liver_10.5-13.5h.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(rain_soleus_rep_12, file = "RAIN_results/rain_results_soleus_10.5-13.5h.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(rain_gastrocnemius_rep_12, file = "RAIN_results/rain_results_gastrocnemius_10.5-13.5h.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

write.table(rain_liver_rep_24, file = "RAIN_results/rain_results_liver_21-27h.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(rain_soleus_rep_24, file = "RAIN_results/rain_results_soleus_21-27h.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(rain_gastrocnemius_rep_24, file = "RAIN_results/rain_results_gastrocnemius_21-27h.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

write.table(rownames(rain_liver_rep_sig_12), "list_rain_liver_12h.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(rownames(rain_soleus_rep_sig_12), "list_rain_soleus_12h.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(rownames(rain_gastrocnemius_rep_sig_12), "list_rain_gastrocnemius_12h.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# =============================================================================
# ANALYSE DES RESULTATS DE RAIN ####
# =============================================================================

# Chargement des dataset
rain_liver_rep_12 <- read.delim("C:/Users/lucid/OneDrive/Bureau/M1 BBC/02 stage/RNAseq_datas/brb-seq_analysis/RAIN_results/rain_results_liver_10.5-13.5h.txt")
rain_soleus_rep_12 <- read.delim("C:/Users/lucid/OneDrive/Bureau/M1 BBC/02 stage/RNAseq_datas/brb-seq_analysis/RAIN_results/rain_results_soleus_10.5-13.5h.txt")
rain_gastrocnemius_rep_12 <- read.delim("C:/Users/lucid/OneDrive/Bureau/M1 BBC/02 stage/RNAseq_datas/brb-seq_analysis/RAIN_results/rain_results_gastrocnemius_10.5-13.5h.txt")

rain_liver_rep_24 <- read.delim("C:/Users/lucid/OneDrive/Bureau/M1 BBC/02 stage/RNAseq_datas/brb-seq_analysis/RAIN_results/rain_results_liver_21-27h.txt")
rain_soleus_rep_24 <- read.delim("C:/Users/lucid/OneDrive/Bureau/M1 BBC/02 stage/RNAseq_datas/brb-seq_analysis/RAIN_results/rain_results_soleus_21-27h.txt")
rain_gastrocnemius_rep_24 <- read.delim("C:/Users/lucid/OneDrive/Bureau/M1 BBC/02 stage/RNAseq_datas/brb-seq_analysis/RAIN_results/rain_results_gastrocnemius_21-27h.txt")

rain_liver_rep_sig_12 <- rain_liver_rep_12[rain_liver_rep_12$pVal < 0.05,]
rain_liver_rep_sig_24 <- rain_liver_rep_24[rain_liver_rep_24$pVal < 0.05,]

rain_soleus_rep_sig_12 <- rain_soleus_rep_12[rain_soleus_rep_12$pVal < 0.05,]
rain_soleus_rep_sig_24 <- rain_soleus_rep_24[rain_soleus_rep_24$pVal < 0.05,]

rain_gastrocnemius_rep_sig_12 <- rain_gastrocnemius_rep_12[rain_gastrocnemius_rep_12$pVal < 0.05,]
rain_gastrocnemius_rep_sig_24 <- rain_gastrocnemius_rep_24[rain_gastrocnemius_rep_24$pVal < 0.05,]

# Avec un cut-off de pVal plus strict
rain_liver_rep_sig_12 <- rain_liver_rep_12[rain_liver_rep_12$pVal < 0.01,]
rain_soleus_rep_sig_12 <- rain_soleus_rep_12[rain_soleus_rep_12$pVal < 0.01,]
rain_gastrocnemius_rep_sig_12 <- rain_gastrocnemius_rep_12[rain_gastrocnemius_rep_12$pVal < 0.01,]

rain_liver_rep_sig_24 <- rain_liver_rep_24[rain_liver_rep_24$pVal < 0.01,]
rain_soleus_rep_sig_24 <- rain_soleus_rep_24[rain_soleus_rep_24$pVal < 0.01,]
rain_gastrocnemius_rep_sig_24 <- rain_gastrocnemius_rep_24[rain_gastrocnemius_rep_24$pVal < 0.01,]


# =============================================================================
# DONNEES MALES VS DONNEES FEMELLES ####
# =============================================================================

# Les données "males" sont celles publiées dans le Zhu et al 2020
# Les données "femelles" sont celles obtenues par l'équipe

# Gènes circadiens
list_24h <- list(
    "Male 24h (Zhu data)" = rownames(rain_zhu_24_sig),
    "Female 24h (Team data)" = rownames(rain_liver_rep_sig_24)
    )
venn_24h_liver <- venn.diagram(
    x = list_24h,
    filename = NULL,  # Ne pas sauvegarder en fichier
    col = "transparent",
    fill = viridis(length(list_24h), option = "C"),
    alpha = 0.5,
    label.col = "black",
    cex = 2,
    fontfamily = "sans",
    cat.cex = 2,
    cat.fontfamily = "sans",
    cat.pos = c(-10, 15),  # Positions personnalisées pour les noms des ensembles
    cat.dist = c(0.03, 0.07) # Distance des noms par rapport aux cercles
    )  
grid::grid.newpage()
grid::grid.draw(venn_24h_liver)

# Gènes ultradiens
list_12h_liver <- list(
    "Male 12h (Zhu data)" = rownames(rain_zhu_12_sig),
    "Female 12h \n(Team data)" = rownames(rain_liver_rep_sig_12)
    )
venn_12h_liver <- venn.diagram(
    x = list_12h_liver,
    filename = NULL,  # Ne pas sauvegarder en fichier
    col = "transparent",
    fill = viridis(length(list_12h_liver), option = "C"),
    alpha = 0.5,
    label.col = "black",
    cex = 2,
    fontfamily = "sans",
    cat.cex = 2,
    cat.fontfamily = "sans",
    cat.pos = c(-10, 10),  # Positions personnalisées pour les noms des ensembles
    cat.dist = c(0.03, 0.07) # Distance des noms par rapport aux cercles
    )  
grid::grid.newpage()
grid::grid.draw(venn_12h_liver)

# =============================================================================
# GENES ULTRADIENS VS CIRCADIENS PAR TISSU ####
# =============================================================================

# Foie - Gènes ultradiens VS gènes circadiens
liver_list <- list(
    "12h liver" = rownames(rain_liver_rep_sig_12),
    "24h liver" = rownames(rain_liver_rep_sig_24)
    )
venn_liver <- venn.diagram(
    x = liver_list,
    filename = NULL,  # Ne pas sauvegarder en fichier
    col = "transparent",
    fill = viridis(length(liver_list), option = "C"),
    alpha = 0.5,
    label.col = "black",
    cex = 2,
    fontfamily = "sans",
    cat.cex = 2,
    cat.fontfamily = "sans",
    cat.pos = c(-5, 9),  # Positions personnalisées pour les noms des ensembles
    cat.dist = c(-0.25, -0.43) # Distance des noms par rapport aux cercles
    ) 
grid::grid.newpage()
grid::grid.draw(venn_liver)

# Soleus - Gènes ultradiens VS gènes circadiens
soleus_list <- list(
    "12h soleus" = rownames(rain_soleus_rep_sig_12),
    "24h soleus" = rownames(rain_soleus_rep_sig_24)
    )
venn_soleus <- venn.diagram(
    x = soleus_list,
    filename = NULL,  # Ne pas sauvegarder en fichier
    col = "transparent",
    fill = viridis(length(soleus_list), option = "C"),
    alpha = 0.5,
    label.col = "black",
    cex = 2,
    fontfamily = "sans",
    cat.cex = 2,
    cat.fontfamily = "sans",
    cat.pos = c(-5, 9),  # Positions personnalisées pour les noms des ensembles
    cat.dist = c(-0.35, -0.43)
    )
grid::grid.newpage()
grid::grid.draw(venn_soleus)

# Gastrocnemius - Gènes ultradiens VS gènes circadiens
gastrocnemius_list <- list(
    "12h gastrocnemius" = rownames(rain_gastrocnemius_rep_sig_12),
    "24h gastrocnemius" = rownames(rain_gastrocnemius_rep_sig_24)
    )
venn_gastrocnemius <- venn.diagram(
    x = gastrocnemius_list,
    filename = NULL,  # Ne pas sauvegarder en fichier
    col = "transparent",
    fill = viridis(length(gastrocnemius_list), option = "C"),
    alpha = 0.5,
    label.col = "black",
    cex = 2,
    fontfamily = "sans",
    cat.cex = 2,
    cat.fontfamily = "sans",
    cat.pos = c(-5, 9),  # Positions personnalisées pour les noms des ensembles
    cat.dist = c(-0.4, -0.44)
    )
grid::grid.newpage()
grid::grid.draw(venn_gastrocnemius)

# =============================================================================
# GENES ULTRADIENS PAR TISSU, GENES CIRCADIENS PAR TISSU ####
# =============================================================================

# Gènes ultradiens
list_tissues_12 <- list(
    "liver 12h" = rownames(rain_liver_rep_sig_12),
    "soleus 12h" = rownames(rain_soleus_rep_sig_12),
    "gastrocnemius 12h" = rownames(rain_gastrocnemius_rep_sig_12)
    )
venn_tissues_12 <- venn.diagram(
    x = list_tissues_12,
    filename = NULL,  # Ne pas sauvegarder en fichier
    col = "transparent",
    fill = viridis(length(list_tissues_12), option = "C"),
    alpha = 0.5,
    label.col = "black",
    cex = 2,
    fontfamily = "sans",
    cat.cex = 2,
    cat.fontfamily = "sans"
    )
grid::grid.newpage()
grid::grid.draw(venn_tissues_12)

# Gènes circadiens
list_tissues_24 <- list(
    "liver 24h" = rownames(rain_liver_rep_sig_24),
    "soleus 24h" = rownames(rain_soleus_rep_sig_24),
    "gastrocnemius 24h" = rownames(rain_gastrocnemius_rep_sig_24)
    )
venn_tissues_24 <- venn.diagram(
    x = list_tissues_24,
    filename = NULL,  # Ne pas sauvegarder en fichier
    col = "transparent",
    fill = viridis(length(list_tissues_24), option = "C"),
    alpha = 0.5,
    label.col = "black",
    cex = 2,
    fontfamily = "sans",
    cat.cex = 2,
    cat.fontfamily = "sans"
    )
grid::grid.newpage()
grid::grid.draw(venn_tissues_24)

# Plus rapide, plus simple : upsetplot !
gene_sets <- list(
  "12h liver" = rownames(rain_liver_rep_sig_12),
  "24h liver" = rownames(rain_liver_rep_sig_24),
  "12h soleus" = rownames(rain_soleus_rep_sig_12),
  "24h soleus" = rownames(rain_soleus_rep_sig_24),
  "12h gastrocnemius" = rownames(rain_gastrocnemius_rep_sig_12),
  "24h gastrocnemius" = rownames(rain_gastrocnemius_rep_sig_24)
)
rhythmic_genes <- unique(unlist(gene_sets))
matrix_rhythmic_genes <- data.frame(gene = rhythmic_genes)
for (set_name in names(gene_sets)) {
  matrix_rhythmic_genes[[set_name]] <- ifelse(matrix_rhythmic_genes$gene %in% gene_sets[[set_name]], 1, 0)
}
rownames(matrix_rhythmic_genes) <- matrix_rhythmic_genes$gene
matrix_rhythmic_genes$gene <- NULL

# Visualisation
upset(
    matrix_rhythmic_genes,
    sets = names(gene_sets),
    order.by = "freq",
    point.size = 5,  
    line.size = 1.2,
    matrix.dot.alpha = 0.8,
    text.scale = 2,
    nintersects = NA,
    set_size.show = TRUE
    )

# =============================================================================
# EXTRACTION DES LISTES DE GENES RYTHMIQUES ####
# =============================================================================

# Gènes circadiens
list_all_tissues_24 <- Reduce(
    intersect,
    list(rownames(matrix_log_gastrocnemius_rep), # n importe quelle matrice qui a tous les genes dedans
        rownames(rain_liver_rep_sig_24),
        rownames(rain_soleus_rep_sig_24),
        rownames(rain_gastrocnemius_rep_sig_24)
    )
)
print(list_all_tissues_24)

# Gènes ultradiens
# Avec intersection gènes circadiens inclue
inter_12h <- Reduce(intersect, list_tissues_12)

# Sans l'intersetion
genes_24 <- unique(c(
    rownames(rain_liver_rep_sig_24),
    rownames(rain_soleus_rep_sig_24),
    rownames(rain_gastrocnemius_rep_sig_24)
    )
)

genes_12 <- Reduce(
    intersect,
    list(rownames(matrix_log_gastrocnemius_rep), # n importe quelle matrice qui a tous les genes dedans
        rownames(rain_liver_rep_sig_12),
        rownames(rain_soleus_rep_sig_12),
        rownames(rain_gastrocnemius_rep_sig_12)
    )
)

list_all_tissues_12 <- setdiff(genes_12, genes_24)
print(list_all_tissues_12)

# Sauvegarder la liste de gènes rythmiques
write.table(list_all_tissues_12, file = "rhythmic_common_genes_12h.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# =============================================================================
# HEATMAPS ####
# =============================================================================

# Foie - heatmaps ####

# Création des dataframe avec les valeurs moyennes pour chaque ZT
matrix_log_liver_mean <- calculer_moyennes_ZT(matrix_log_liver_rep)
matrix_log_liver_rain_12 <- matrix_log_liver_mean[rownames(matrix_log_liver_mean) %in% rownames(rain_liver_rep_sig_12),]
matrix_log_liver_rain_24 <- matrix_log_liver_mean[rownames(matrix_log_liver_mean) %in% rownames(rain_liver_rep_sig_24),]

# Gènes ultradiens - classés par corrélation de Pearson
heatmap_liver_12 <- pheatmap(
  matrix_log_liver_rain_12,
  cluster_rows = TRUE, # Cluster the rows of the heatmap (genes in this case)
  cluster_cols = FALSE, # keep samples in order
  show_rownames = FALSE, # There are too many genes to clearly show the labels
  main = "Non-Annotated Heatmap : Liver, ~12h rhythmic genes",
  clustering_distance_rows = "correlation", # method for hclust
  color = viridis(9),
  scale = "row" # Scale values in the direction of genes (rows)
) 

# Gènes circadiens - classés par corrélation de Pearson
heatmap_liver_24 <- pheatmap(
  matrix_log_liver_rain_24,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  main = "Non-Annotated Heatmap : Liver, ~24h rhythmic genes",
  clustering_distance_rows = "correlation",
  color = viridis(9),
  scale = "row"
)

# Groupement des heatmaps
heatmaps_liver <- plot_grid(
    heatmap_liver_12$gtable, heatmap_liver_24$gtable, 
    ncol = 2, 
    labels = c("A", "B"), 
    align = "h"
)
heatmaps_liver

# Test de clusterisation forcée
liver_test <- cbind(matrix_log_liver_rain_12, rain_liver_rep_sig_12)
matrix_liver_mean_ord <- liver_test %>%
  arrange(phase, period, pVal)
matrix_liver_mean_ord <- subset(matrix_liver_mean_ord, select = -c(period, phase, pVal, peak.shape))

heatmap_liver_12 <- pheatmap(
  matrix_liver_mean_ord,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  main = "Non-Annotated Heatmap : liver, ~12h rhythmic genes",
  color = viridis(9),
  scale = "row"
)

# Soleus - heatmaps ####

# Création des dataframe avec les valeurs moyennes pour chaque ZT
matrix_log_soleus_mean <- calculer_moyennes_ZT(matrix_log_soleus_rep)
matrix_log_soleus_rain_12 <- matrix_log_soleus_mean[rownames(matrix_log_soleus_mean) %in% rownames(rain_soleus_rep_sig_12),]
matrix_log_soleus_rain_24 <- matrix_log_soleus_mean[rownames(matrix_log_soleus_mean) %in% rownames(rain_soleus_rep_sig_24),]

# Gènes ultradiens - classés par corrélation de Pearson
heatmap_soleus_12 <- pheatmap(
  matrix_log_soleus_rain_12,
  cluster_rows = TRUE, # Cluster the rows of the heatmap (genes in this case)
  cluster_cols = FALSE, # keep samples in order
  show_rownames = FALSE, # There are too many genes to clearly show the labels
  main = "Non-Annotated Heatmap : soleus, ~12h rhythmic genes",
  clustering_distance_rows = "correlation", # method for hclust
  color = viridis(9),
  scale = "row" # Scale values in the direction of genes (rows)
) 

# Gènes circadiens - classés par corrélation de Pearson
heatmap_soleus_24 <- pheatmap(
  matrix_log_soleus_rain_24,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  main = "Non-Annotated Heatmap : soleus, ~24h rhythmic genes",
  clustering_distance_rows = "correlation",
  color = viridis(9),
  scale = "row"
)

# Groupement des heatmaps
heatmaps_soleus <- plot_grid(
    heatmap_soleus_12$gtable, heatmap_soleus_24$gtable, 
    ncol = 2, 
    labels = c("A", "B"), 
    align = "h"
)
heatmaps_soleus

# Test de clusterisation forcée
soleus_test <- cbind(matrix_log_soleus_rain_12, rain_soleus_rep_sig_12)
matrix_soleus_mean_ord <- soleus_test %>%
  arrange(phase, period, pVal)
matrix_soleus_mean_ord <- subset(matrix_soleus_mean_ord, select = -c(period, phase, pVal, peak.shape))

heatmap_soleus_12 <- pheatmap(
  matrix_soleus_mean_ord,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  main = "Non-Annotated Heatmap : soleus, ~12h rhythmic genes",
  color = viridis(9),
  scale = "row"
)

# Gastrocnemius - heatmaps ####

# Création des dataframe avec les valeurs moyennes pour chaque ZT
matrix_log_gastrocnemius_mean <- calculer_moyennes_ZT(matrix_log_gastrocnemius_rep)
matrix_log_gastrocnemius_rain_12 <- matrix_log_gastrocnemius_mean[rownames(matrix_log_gastrocnemius_mean) %in% rownames(rain_gastrocnemius_rep_sig_12),]
matrix_log_gastrocnemius_rain_24 <- matrix_log_gastrocnemius_mean[rownames(matrix_log_gastrocnemius_mean) %in% rownames(rain_gastrocnemius_rep_sig_24),]

# Gènes ultradiens - classés par corrélation de Pearson
heatmap_gastrocnemius_12 <- pheatmap(
  matrix_log_gastrocnemius_rain_12,
  cluster_rows = TRUE, # Cluster the rows of the heatmap (genes in this case)
  cluster_cols = FALSE, # keep samples in order
  show_rownames = FALSE, # There are too many genes to clearly show the labels
  main = "Non-Annotated Heatmap : gastrocnemius, ~12h rhythmic genes",
  clustering_distance_rows = "correlation", # method for hclust
  color = viridis(9),
  scale = "row" # Scale values in the direction of genes (rows)
) 

# Gènes circadiens - classés par corrélation de Pearson
heatmap_gastrocnemius_24 <- pheatmap(
  matrix_log_gastrocnemius_rain_24,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  main = "Non-Annotated Heatmap : gastrocnemius, ~24h rhythmic genes",
  clustering_distance_rows = "correlation",
  color = viridis(9),
  scale = "row"
)

# Groupement des heatmaps
heatmaps_gastrocnemius <- plot_grid(
    heatmap_gastrocnemius_12$gtable, heatmap_gastrocnemius_24$gtable, 
    ncol = 2, 
    labels = c("A", "B"), 
    align = "h"
)
heatmaps_gastrocnemius

# Test de clusterisation forcée
gastrocnemius_test <- cbind(matrix_log_gastrocnemius_rain_12, rain_gastrocnemius_rep_sig_12)
matrix_gastrocnemius_mean_ord <- gastrocnemius_test %>%
  arrange(phase, period, pVal)
matrix_gastrocnemius_mean_ord <- subset(matrix_gastrocnemius_mean_ord, select = -c(period, phase, pVal, peak.shape))

heatmap_gastrocnemius_12 <- pheatmap(
  matrix_gastrocnemius_mean_ord,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  main = "Non-Annotated Heatmap : gastrocnemius, ~12h rhythmic genes",
  color = viridis(9),
  scale = "row"
)

# Visualisation des heatmaps des 3 tissus ensembles
resume_heatmaps <- plot_grid(heatmap_liver_12$gtable, heatmap_soleus_12$gtable, heatmap_gastrocnemius_12$gtable,
                             heatmap_liver_24$gtable, heatmap_soleus_24$gtable, heatmap_gastrocnemius_24$gtable,
                             ncol = 3, align = "h")
resume_heatmaps

# =============================================================================
# ANALYSE DES DISTRIBUTIONS DE PHASES ####
# =============================================================================

# Distribution de phase dans le foie ####

# Gènes ultradiens
matrix_phase_liver_rain <- subset(rain_liver_rep_12, pVal < 0.05, select = phase)
matrix_phase_liver_rain$gene <- rownames(matrix_phase_liver_rain)

hist_liver_12 <- circular_phase12H_histogram(
    matrix_phase_liver_rain,
    "phase",
    "gene",
    title_override = "~12h rhythmic genes",
    fill_color = "royalblue3",
    bg_color1 = "white",
    bg_color2 = "white"
)
hist_liver_12

# Gèens circadiens
matrix_phase_liver_rain_24 <- subset(rain_liver_rep_24, pVal < 0.05, select = phase)
matrix_phase_liver_rain_24$gene <- rownames(matrix_phase_liver_rain_24)

hist_liver_24 <- circular_phase24H_histogram(
    matrix_phase_liver_rain_24,
    "phase",
    "gene",
    title_override = "~24h rhythmic genes",
    fill_color = "royalblue3",
    bg_color1 = "lightgoldenrodyellow",
    bg_color2 = "lightsteelblue2"
)
hist_liver_24

# Distribution de phase dans le soleus ####

# Gènes ultradiens
matrix_phase_soleus_rain <- subset(rain_soleus_rep_12, pVal < 0.05, select = phase)
matrix_phase_soleus_rain$gene <- rownames(matrix_phase_soleus_rain)

hist_soleus_12 <- circular_phase12H_histogram(
    matrix_phase_soleus_rain,
    "phase",
    "gene",
    title_override = "~12h rhythmic genes",
    fill_color = "firebrick3",
    bg_color1 = "white",
    bg_color2 = "white"
)
hist_soleus_12

# Gèens circadiens
matrix_phase_soleus_rain_24 <- subset(rain_soleus_rep_24, pVal < 0.05, select = phase)
matrix_phase_soleus_rain_24$gene <- rownames(matrix_phase_soleus_rain_24)

hist_soleus_24 <- circular_phase24H_histogram(
    matrix_phase_soleus_rain_24,
    "phase",
    "gene",
    title_override = "~24h rhythmic genes",
    fill_color = "firebrick3",
    bg_color1 = "lightgoldenrodyellow",
    bg_color2 = "lightsteelblue2"
)
hist_soleus_24

# Distribution de phase dans le gastrocnemius ####

# Gènes ultradiens
matrix_phase_gastrocnemius_rain <- subset(rain_gastrocnemius_rep_12, pVal < 0.05, select = phase)
matrix_phase_gastrocnemius_rain$gene <- rownames(matrix_phase_gastrocnemius_rain)

hist_gastrocnemius_12 <- circular_phase12H_histogram(
    matrix_phase_gastrocnemius_rain,
    "phase",
    "gene",
    title_override = "~12h rhythmic genes",
    fill_color = "palevioletred3",
    bg_color1 = "white",
    bg_color2 = "white"
)
hist_gastrocnemius_12

# Gèens circadiens
matrix_phase_gastrocnemius_rain_24 <- subset(rain_gastrocnemius_rep_24, pVal < 0.05, select = phase)
matrix_phase_gastrocnemius_rain_24$gene <- rownames(matrix_phase_gastrocnemius_rain_24)

hist_gastrocnemius_24 <- circular_phase24H_histogram(
    matrix_phase_gastrocnemius_rain_24,
    "phase",
    "gene",
    title_override = "~24h rhythmic genes",
    fill_color = "palevioletred3",
    bg_color1 = "lightgoldenrodyellow",
    bg_color2 = "lightsteelblue2"
)
hist_gastrocnemius_24

# Groupement des histogrammes
hist_phases <- (hist_liver_12 + hist_soleus_12 + hist_gastrocnemius_12) / (hist_liver_24 + hist_soleus_24 + hist_gastrocnemius_24)+
  plot_annotation(title = "Phases distribution of rhythmic genes upon the 3 tissues")
hist_phases

# =============================================================================
# TOP12 DES GENES RYTHMIQUES PAR TISSU ####
# =============================================================================
# Extraction des 12 gènes avec la pVal la plus faible pour chaque tissu

# Top12 gènes ultradiens dans le foie ###
top_12h_liver <- rain_liver_rep_sig_12 %>%
  rownames_to_column("Gene") %>%
  slice_min(order_by = adjpval, n = 12, with_ties = FALSE) %>%
  pull(Gene)

# Transformation du dataframe en df long
gene_data_liver <- matrix_log_liver_rep[rownames(matrix_log_liver_rep), , drop = FALSE]
gene_long_liver <- gene_data_liver %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "ZT", values_to = "expression") %>%
  mutate(
    Zeitgeber = as.numeric(str_extract(ZT, "\\d+")),
    Replicate = str_extract(ZT, "[ab]"))

# Extraction du df long contenant uniquement les gènes du top12    
resume_donnees_liver <- gene_long_liver %>%
  filter(Gene %in% top_12h_liver) %>%
  group_by(Gene, Zeitgeber) %>%
  summarise(
    expression_moyenne = mean(expression, na.rm = TRUE),
    erreur_std = sd(expression, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

# Visualisation des expressions moyennes
plot_top12_liver <- ggplot(resume_donnees_liver, aes(x = Zeitgeber, y = expression_moyenne, group = Gene)) +
  geom_line(color="royalblue3", linewidth = 1.5) +
  geom_point(size = 3, color="royalblue3") +
  geom_errorbar(aes(ymin = expression_moyenne - erreur_std, 
                    ymax = expression_moyenne + erreur_std), 
                width = 0.4, color = "royalblue3") +
  labs(title = "Mean expression of top ~12h rhythmic genes in liver",
       x = "ZT (Zeitgeber Time)",
       y = "mean(log2 normalized counts)") +
  theme_light()+
  theme(strip.text = element_text(size = 12, color = "black"),
        strip.background = element_rect(fill = "gray85"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  facet_wrap(~ Gene, scales = "free_y")
plot_top12_liver

# Top12 gènes ultradiens dans le soleus ###
top_12h_soleus <- rain_soleus_rep_sig_12 %>%
  rownames_to_column("Gene") %>%
  slice_min(order_by = adjpval, n = 12, with_ties = FALSE) %>%
  pull(Gene)

# Transformation du dataframe en df long
gene_data_soleus <- matrix_log_soleus_rep[rownames(matrix_log_soleus_rep), , drop = FALSE]
gene_long_soleus <- gene_data_soleus %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "ZT", values_to = "expression") %>%
  mutate(
    Zeitgeber = as.numeric(str_extract(ZT, "\\d+")),
    Replicate = str_extract(ZT, "[ab]"))

# Extraction du df long contenant uniquement les gènes du top12    
resume_donnees_soleus <- gene_long_soleus %>%
  filter(Gene %in% top_12h_soleus) %>%
  group_by(Gene, Zeitgeber) %>%
  summarise(
    expression_moyenne = mean(expression, na.rm = TRUE),
    erreur_std = sd(expression, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

# Visualisation des expressions moyennes
plot_top12_soleus <- ggplot(resume_donnees_soleus, aes(x = Zeitgeber, y = expression_moyenne, group = Gene)) +
  geom_line(color="firebrick3", linewidth = 1.5) +
  geom_point(size = 3, color="firebrick3") +
  geom_errorbar(aes(ymin = expression_moyenne - erreur_std, 
                    ymax = expression_moyenne + erreur_std), 
                width = 0.4, color = "firebrick3") +
  labs(title = "Mean expression of top ~12h rhythmic genes in soleus",
       x = "ZT (Zeitgeber Time)",
       y = "mean(log2 normalized counts)") +
  theme_light()+
  theme(strip.text = element_text(size = 12, color = "black"),
        strip.background = element_rect(fill = "gray85"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  facet_wrap(~ Gene, scales = "free_y")
plot_top12_soleus

# Top12 des gènes ultradiens dans le gastrocnemius ####
top_12h_gastrocnemius <- rain_gastrocnemius_rep_sig_12 %>%
  rownames_to_column("Gene") %>%
  slice_min(order_by = adjpval, n = 12, with_ties = FALSE) %>%
  pull(Gene)

# Transformation du dataframe en df long
gene_data_gastrocnemius <- matrix_log_gastrocnemius_rep[rownames(matrix_log_gastrocnemius_rep), , drop = FALSE]
gene_long_gastrocnemius <- gene_data_gastrocnemius %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "ZT", values_to = "expression") %>%
  mutate(
    Zeitgeber = as.numeric(str_extract(ZT, "\\d+")),
    Replicate = str_extract(ZT, "[ab]"))

# Extraction du df long contenant uniquement les gènes du top12    
resume_donnees_gastrocnemius <- gene_long_gastrocnemius %>%
  filter(Gene %in% top_12h_gastrocnemius) %>%
  group_by(Gene, Zeitgeber) %>%
  summarise(
    expression_moyenne = mean(expression, na.rm = TRUE),
    erreur_std = sd(expression, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

# Visualisation des expressions moyennes
plot_top12_gastrocnemius <- ggplot(resume_donnees_gastrocnemius, aes(x = Zeitgeber, y = expression_moyenne, group = Gene)) +
  geom_line(color="palevioletred3", linewidth = 1.5) +
  geom_point(size = 3, color="palevioletred3") +
  geom_errorbar(aes(ymin = expression_moyenne - erreur_std, 
                    ymax = expression_moyenne + erreur_std), 
                width = 0.4, color = "palevioletred3") +
  labs(title = "Mean expression of top ~12h rhythmic genes in gastrocnemius",
       x = "ZT (Zeitgeber Time)",
       y = "mean(log2 normalized counts)") +
  theme_light()+
  theme(strip.text = element_text(size = 12, color = "black"),
        strip.background = element_rect(fill = "gray85"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  facet_wrap(~ Gene, scales = "free_y")
plot_top12_gastrocnemius


# =============================================================================
# GENES ULTRADIENS COMMUNS AUX 3 TISSUS ####
# =============================================================================

palette_tissue <- c("liver" = "royalblue3", "soleus" = "firebrick3", "gastrocnemius" = "palevioletred3")

# Extraction des expressions moyennes des genes communs dans chaque tissu
matrix_top_liver <- gene_long %>%
  filter(Gene %in% list_all_tissues_12) %>%
  mutate(tissue = "liver")
matrix_top_soleus <- gene_long_soleus %>%
  filter(Gene %in% list_all_tissues_12) %>%
  mutate(tissue = "soleus")
matrix_top_gastrocnemius <- gene_long_gastrocnemius %>%
  filter(Gene %in% list_all_tissues_12) %>%
  mutate(tissue = "gastrocnemius")

# Fusionner les dataframes
matrix_top <- bind_rows(matrix_top_liver, matrix_top_soleus, matrix_top_gastrocnemius)

# Retransformer en dataframe long
resume_donnees <- matrix_top %>%
  group_by(Gene, Zeitgeber, tissue) %>%
  summarise(
    expression_moyenne = mean(expression, na.rm = TRUE),
    erreur_std = sd(expression, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

# Visualisation
plot_top <- ggplot(resume_donnees, aes(x = Zeitgeber, y = expression_moyenne, group = Gene, color = tissue)) +
  geom_line(aes(group = tissue), linewidth = 1.5) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = expression_moyenne - erreur_std, 
                    ymax = expression_moyenne + erreur_std), 
                width = 0.4) +
  scale_color_manual(values = palette_tissue)+
  labs(title = "Mean expression of the 8 common 12h rhythmic genes",
       x = "ZT (Zeitgeber Time)",
       y = "mean(log2 normalized counts)") +
  theme_light(base_size = 16)+
  theme(strip.text = element_text(size = 16, color = "black"),
        strip.background = element_rect(fill = "gray85"),
        legend.position = "bottom",
        legend.text = element_text(size=16))+
  facet_wrap(~ Gene, scales = "free_y", nrow = 2)
plot_top

# =============================================================================
# ANALYSE D'ENRICHISSEMENT ####
# =============================================================================

# Objets pour enrichR
websiteLive <- getOption("enrichR.live")

if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human and mouse genes
}

dbs <- c(
    "GO_Biological_Process_2025",
    "GO_Cellular_Component_2025",
    "GO_Molecular_Function_2025",
    "KEGG_2019_Mouse",
    "Pfam_Domains_2023",
    "WikiPathways_2024_Mouse",
    "Reactome_Pathways_2024",
    "Mouse_Gene_Atlas"
    )

# Enrichissement d'annotations dans le foie ####
list_liver_rain_12 <- rownames(rain_liver_rep_sig_12)

# Analyse pour extraire les listes
enrichr_liver_rain_12 <- enrichr(list_liver_rain_12, dbs)
# Enregistrer les listes
printEnrich(enrichr_liver_rain_12, prefix = "liver_12_enrichment", outFile = "excel")
# Analyse pour extraire les plots
result_enrichr_liver <- plot_enrichment(list_liver_rain_12, dbs, num_terms = 10)


# Enrichissement d'annotations dans le soleus ####
list_soleus_rain_12 <- rownames(rain_soleus_rep_sig_12)

# Analyse pour extraire les listes
enrichr_soleus_rain_12 <- enrichr(list_soleus_rain_12, dbs)
# Enregistrer les listes
printEnrich(enrichr_soleus_rain_12, prefix = "soleus_12_enrichment", outFile = "excel")
# Analyse pour extraire les plots
result_enrichr_soleus <- plot_enrichment(list_soleus_rain_12, dbs, num_terms = 10)


# Enrichissement d'annotations dans le gastrocnemius ####
list_gastrocnemius_rain_12 <- rownames(rain_gastrocnemius_rep_sig_12)

# Analyse pour extraire les listes
enrichr_gastrocnemius_rain_12 <- enrichr(list_gastrocnemius_rain_12, dbs)
# Enregistrer les listes
printEnrich(enrichr_gastrocnemius_rain_12, prefix = "gastrocnemius_12_enrichment", outFile = "excel")
# Analyse pour extraire les plots
result_enrichr_gastrocnemius <- plot_enrichment(list_gastrocnemius_rain_12, dbs, num_terms = 10)

# Visualiser une liste
View(enrichr_gastrocnemius_rain_12[["GO_Biological_Process_2025"]])

# =============================================================================
# FOCUS SUR LES MUSCLES SQUELETTIQUES ####
# =============================================================================

# Refaire une liste des listes de gènes rythmiques
list_genes <- list(
  liver_rep_sig_12 = rownames(rain_liver_rep_sig_12),
  liver_rep_sig_24 = rownames(rain_liver_rep_sig_24),
  soleus_rep_sig_12 = rownames(rain_soleus_rep_sig_12),
  soleus_rep_sig_24 = rownames(rain_soleus_rep_sig_24),
  gastrocnemius_rep_sig_12 = rownames(rain_gastrocnemius_rep_sig_12),
  gastrocnemius_rep_sig_24 = rownames(rain_gastrocnemius_rep_sig_24)
)

# Extraction des gègnes communs uniquement aux muscles squelettiques
genes_skm_12 <- intersect(rownames(rain_soleus_rep_sig_12), rownames(rain_gastrocnemius_rep_sig_12))
other_dfs <- setdiff(names(list_genes), c("soleus_rep_sig_12", "gastrocnemius_rep_sig_12"))
genes_in_others <- unique(unlist(list_genes[other_dfs]))
list_12h_skm <- setdiff(genes_skm_12, genes_in_others)

# Analyse d'enrichissement
enrichr_skm_rain_12 <- enrichr(list_12h_skm, dbs)
result_enrichr_skm_12h_genes <- plot_enrichment(list_12h_skm, dbs, num_terms = 6)

# =============================================================================
# FOCUS SUR UN SEUL GENE ####
# =============================================================================
palette_tissue <- c("liver" = "royalblue3", "soleus" = "firebrick3", "gastrocnemius" = "palevioletred3")

# Choisir un gene d'interet
gene_of_interest <- c("Xbp1")

# Extraire les valeurs moyennes d'expression des dataframes longs pour ce gene
matrix_interest_liver <- gene_long_liver %>%
  filter(Gene %in% gene_of_interest) %>%
  mutate(tissue = "liver")
matrix_interest_soleus <- gene_long_soleus %>%
  filter(Gene %in% gene_of_interest) %>%
  mutate(tissue = "soleus")
matrix_interest_gastrocnemius <- gene_long_gastrocnemius %>%
  filter(Gene %in% gene_of_interest) %>%
  mutate(tissue = "gastrocnemius")

# Fusionner les dataframes longs
matrix_interest <- bind_rows(matrix_interest_liver, matrix_interest_soleus, matrix_interest_gastrocnemius)

# Transformer le dataframe
# on calcule l'expression moyenne et l'ecart type à chaque ZT
resume_interest <- matrix_interest %>%
  group_by(Gene, Zeitgeber, tissue) %>%
  summarise(
    expression_moyenne = mean(expression, na.rm = TRUE),
    erreur_std = sd(expression, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

# Visualisation
plot_one_gene <- ggplot(resume_interest, aes(x = Zeitgeber, y = expression_moyenne, group = Gene, color = tissue)) +
  geom_line(aes(group = tissue), linewidth = 1.5) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = expression_moyenne - erreur_std, 
                    ymax = expression_moyenne + erreur_std), 
                width = 0.4) +
  scale_color_manual(values = palette_tissue)+
  labs(title = paste0("Mean expression of ", gene_of_interest, " upon time"),
       x = "ZT (Zeitgeber Time)",
       y = "mean(log2 normalized counts)") +
  theme_light(base_size = 16)+
  theme(strip.text = element_text(size = 16, color = "black"),
        strip.background = element_rect(fill = "gray85"),,
        legend.position = "bottom",
        legend.text = element_text(size=16))+
  facet_wrap(~ Gene, scales = "free_y", nrow = 2)
plot_one_gene

