# =============================================================================
# Script R pour le préprocessing des données BRB (Bulk RNA-seq)
# Analyse de données transcriptomiques multi-tissus avec timepoints circadiens
# =============================================================================

# =============================================================================
# 1. CHARGEMENT DES PACKAGES
# =============================================================================

library(data.table)          # Manipulation efficace des données
library(Matrix)              # Gestion des matrices sparses
library(dplyr)               # Manipulation des données
library(tidyverse)           # Collection de packages pour l'analyse de données
library(DESeq2)              # Analyse différentielle d'expression
library(SummarizedExperiment) # Structure de données pour l'omique

# =============================================================================
# 2. CHARGEMENT DES DONNÉES
# =============================================================================

# Définition du répertoire contenant les fichiers Matrix Market
matrix_dir <- "path/to/repertory"

# Chargement des métadonnées des gènes (IDs Ensembl et symboles)
features <- fread(file.path(matrix_dir, "features.tsv"), 
                  header = FALSE, 
                  data.table = FALSE)

# Chargement des codes-barres des échantillons
barcodes <- fread(file.path(matrix_dir, "barcodes.tsv"), 
                  header = FALSE, 
                  data.table = FALSE)

# Chargement des annotations des codes-barres (correspondance barcode -> échantillon)
barcode_annot <- read.delim("path/to/barcodes_row.txt")

# =============================================================================
# 3. CRÉATION DES MATRICES D'EXPRESSION
# =============================================================================

# --- Matrice sans déduplication des UMI ---
cat("Chargement de la matrice sans déduplication...\n")
matrix1 <- readMM(file.path(matrix_dir, "umiDedup-NoDedup.mtx")) %>%
  as.matrix() %>%
  as.data.frame()

# Attribution des noms de colonnes (codes-barres) et de lignes (gènes)
colnames(matrix1) <- barcodes$V1
rownames(matrix1) <- features$V1

# Sauvegarde de la matrice brute
fwrite(matrix1, 
       file = "raw_matrix_umi_no_dedup.txt", 
       sep = "\t", 
       quote = FALSE, 
       row.names = TRUE)

# --- Matrice avec 1 mismatch autorisé ---
cat("Chargement de la matrice avec 1 mismatch...\n")
matrix2 <- readMM(file.path(matrix_dir, "umiDedup-1MM_All.mtx")) %>%
  as.matrix() %>%
  as.data.frame()

# Attribution des noms
colnames(matrix2) <- barcodes$V1
rownames(matrix2) <- features$V1

# Sauvegarde
fwrite(matrix2, 
       file = "raw_matrix_umi_1_MM.txt", 
       sep = "\t", 
       quote = FALSE, 
       row.names = TRUE)

# =============================================================================
# 4. CONVERSION DES CODES-BARRES EN IDENTIFIANTS D'ÉCHANTILLONS
# =============================================================================

# Remplacement des codes-barres par les identifiants d'échantillons
matrix1_named <- matrix1 %>%
  rename_with(~ barcode_annot$sample_cycle[match(., barcode_annot$barcode)])

# =============================================================================
# 5. CONVERSION DES IDs ENSEMBL EN SYMBOLES DE GÈNES
# =============================================================================

# Extraction des correspondances ID Ensembl <-> Symbole
genes <- features[, c(1, 2)]  # V1: Ensembl ID, V2: symbole de gène
colnames(genes) <- c("ensembl_id", "symbol")

# Ajout des IDs Ensembl comme colonne
matrix1$ensembl_id <- rownames(matrix1)

# Fusion avec les symboles de gènes
matrix1 <- merge(matrix1, genes, by = "ensembl_id", all.x = TRUE)

# Utilisation des symboles comme noms de lignes
rownames(matrix1) <- matrix1$symbol

# =============================================================================
# 6. FILTRAGE DES GÈNES DE FOND (BACKGROUND)
# =============================================================================

# Chargement de la liste des gènes de fond
bg_genes <- scan("path/to/background_genes.txt", what = "", quiet = TRUE)

# Extraction de la sous-matrice des gènes de fond
bg_matrix <- matrix_count[rownames(matrix_count) %in% bg_genes, ]

# Calcul du seuil : maximum d'expression parmi les gènes de fond
seuil_max_bg <- max(bg_matrix)

# Filtrage : conservation des gènes avec expression > seuil de fond
matrix_count <- matrix_count[apply(matrix_count, 1, max) > seuil_max_bg, ]

# Sauvegarde de la matrice filtrée
fwrite(matrix_count, 
       file = "matrix_count_nodedup.txt", 
       sep = "\t", 
       quote = FALSE)

# =============================================================================
# 7. PRÉPARATION DES DONNÉES PAR TISSU
# =============================================================================

# Rechargement de la matrice filtrée
matrix_count <- read.delim("path/to/matrix_count_nodedup.txt")
rownames(matrix_count) <- matrix_count$X
matrix_count <- subset(matrix_count, select = -X)

# Nettoyage des noms de colonnes (remplacement des virgules par des points)
colnames(matrix_count) <- gsub(",", ".", colnames(matrix_count))

# Définition de l'ordre des échantillons (3 tissus x 32 timepoints)
ordre_voulu <- c(
  # Foie (liver) - ZT 0 à 46.5 par pas de 1.5h
  paste0("liv_", c(0, 1.5, 3, 4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15, 16.5, 
                   18, 19.5, 21, 22.5, 24, 25.5, 27, 28.5, 30, 31.5, 33, 
                   34.5, 36, 37.5, 39, 40.5, 42, 43.5, 45, 46.5)),
  
  # Muscle soléaire (soleus) - même série temporelle
  paste0("sol_", c(0, 1.5, 3, 4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15, 16.5, 
                   18, 19.5, 21, 22.5, 24, 25.5, 27, 28.5, 30, 31.5, 33, 
                   34.5, 36, 37.5, 39, 40.5, 42, 43.5, 45, 46.5)),
  
  # Muscle gastrocnémien (gastrocnemius) - même série temporelle
  paste0("gast_", c(0, 1.5, 3, 4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15, 16.5, 
                    18, 19.5, 21, 22.5, 24, 25.5, 27, 28.5, 30, 31.5, 33, 
                    34.5, 36, 37.5, 39, 40.5, 42, 43.5, 45, 46.5))
)

# Réorganisation des colonnes dans l'ordre voulu
matrix_count <- matrix_count[, ordre_voulu]

# --- Création des sous-matrices par tissu ---
cat("Création des matrices par tissu...\n")

# Foie : sélection des colonnes commençant par "liv"
matrix_count_liver <- matrix_count[, grepl("^liv", colnames(matrix_count))]

# Muscle soléaire : sélection des colonnes commençant par "sol"
matrix_count_soleus <- matrix_count[, grepl("^sol", colnames(matrix_count))]

# Muscle gastrocnémien : sélection des colonnes commençant par "gast"
matrix_count_gastrocnemius <- matrix_count[, grepl("^gast", colnames(matrix_count))]

# Sauvegarde des matrices par tissu
fwrite(matrix_count_liver, 
       file = "matrix_count_nodedup_liver.txt", 
       sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

fwrite(matrix_count_soleus, 
       file = "matrix_count_nodedup_soleus.txt", 
       sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

fwrite(matrix_count_gastrocnemius, 
       file = "matrix_count_nodedup_gastrocnemius.txt", 
       sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# =============================================================================
# 8. NORMALISATION GLOBALE (TOUS TISSUS)
# =============================================================================

cat("Normalisation globale avec DESeq2...\n")

# Rechargement des données
matrix_count <- fread("matrix_count_nodedup.txt", data.table = FALSE)
rownames(matrix_count) <- matrix_count$V1
matrix_count <- matrix_count[, -1]
colnames(matrix_count) <- gsub(",", ".", colnames(matrix_count))

# Création des métadonnées (design minimal)
meta <- as.data.frame(t(matrix_count))

# Création de l'objet DESeq2
dds <- DESeqDataSetFromMatrix(countData = matrix_count, 
                              colData = meta, 
                              design = ~1)

# Estimation des paramètres de dispersion
dds <- DESeq(dds)

# Transformation stabilisatrice de variance (VST)
matrix_log <- vst(dds, nsub = 10600) %>% 
  assay() %>% 
  as.data.frame()

# Sauvegarde de la matrice normalisée
fwrite(matrix_log, 
       file = "matrix_log_nodedup.txt", 
       sep = "\t", 
       quote = FALSE)

# =============================================================================
# 9. NORMALISATION PAR TISSU
# =============================================================================

# --- FOIE ---
cat("Normalisation du foie...\n")
matrix_count_liver <- fread("matrix_count_nodedup_liver.txt", data.table = FALSE)
rownames(matrix_count_liver) <- matrix_count_liver$V1
matrix_count_liver <- matrix_count_liver[, -1]
colnames(matrix_count_liver) <- gsub(",", ".", colnames(matrix_count_liver))

# Création de l'objet SummarizedExperiment
sel <- SummarizedExperiment(matrix_count_liver)
metal <- as.data.frame(t(matrix_count_liver))

# Normalisation DESeq2
ddsl <- DESeqDataSetFromMatrix(countData = matrix_count_liver, 
                               colData = metal, 
                               design = ~1)
vsdl <- vst(ddsl, nsub = 10000)
matrix_log_liver <- as.data.frame(assay(vsdl))

# Suppression du préfixe "liv_" des noms de colonnes
colnames(matrix_log_liver) <- gsub("^liv_", "", colnames(matrix_log_liver))

# Sauvegarde
fwrite(matrix_log_liver, 
       file = "matrix_log_liver_nodedup.txt", 
       sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# --- SOLÉAIRE ---
cat("Normalisation du muscle soléaire...\n")
matrix_count_soleus <- fread("matrix_count_nodedup_soleus.txt", data.table = FALSE)
rownames(matrix_count_soleus) <- matrix_count_soleus$V1
matrix_count_soleus <- matrix_count_soleus[, -1]
colnames(matrix_count_soleus) <- gsub(",", ".", colnames(matrix_count_soleus))

sel <- SummarizedExperiment(matrix_count_soleus)
metal <- as.data.frame(t(matrix_count_soleus))
ddsl <- DESeqDataSetFromMatrix(countData = matrix_count_soleus, 
                               colData = metal, 
                               design = ~1)
vsdl <- vst(ddsl, nsub = 10000)
matrix_log_soleus <- as.data.frame(assay(vsdl))

# Note: Bug dans le code original - devrait être "^sol_" au lieu de "^liv_"
colnames(matrix_log_soleus) <- gsub("^sol_", "", colnames(matrix_log_soleus))

fwrite(matrix_log_soleus, 
       file = "matrix_log_soleus_nodedup.txt", 
       sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# --- GASTROCNÉMIEN ---
cat("Normalisation du muscle gastrocnémien...\n")
matrix_count_gastrocnemius <- fread("matrix_count_nodedup_gastrocnemius.txt", data.table = FALSE)
rownames(matrix_count_gastrocnemius) <- matrix_count_gastrocnemius$V1
matrix_count_gastrocnemius <- matrix_count_gastrocnemius[, -1]
colnames(matrix_count_gastrocnemius) <- gsub(",", ".", colnames(matrix_count_gastrocnemius))

sel <- SummarizedExperiment(matrix_count_gastrocnemius)
metal <- as.data.frame(t(matrix_count_gastrocnemius))
ddsl <- DESeqDataSetFromMatrix(countData = matrix_count_gastrocnemius, 
                               colData = metal, 
                               design = ~1)
vsdl <- vst(ddsl, nsub = 10000)
matrix_log_gastrocnemius <- as.data.frame(assay(vsdl))

# Note: Bug dans le code original - devrait être "^gast_" au lieu de "^liv_"
colnames(matrix_log_gastrocnemius) <- gsub("^gast_", "", colnames(matrix_log_gastrocnemius))

fwrite(matrix_log_gastrocnemius, 
       file = "matrix_log_gastrocnemius_nodedup.txt", 
       sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# =============================================================================
# 10. RÉORGANISATION EN RÉPLICATS BIOLOGIQUES
# =============================================================================

# --- FOIE ---
cat("Réorganisation des réplicats - Foie...\n")

# Renommage des colonnes pour identifier les réplicats
# ZT 0-22.5 = réplicat A, ZT 24-46.5 = réplicat B (cycle de 24h répété)
matrix_log_liver_rep <- matrix_log_liver %>%
  rename(
    ZT_0a = X0, ZT_1.5a = X1.5, ZT_3a = X3, ZT_4.5a = X4.5, ZT_6a = X6,   
    ZT_7.5a = X7.5, ZT_9a = X9, ZT_10.5a = X10.5, ZT_12a = X12, ZT_13.5a = X13.5,
    ZT_15a = X15, ZT_16.5a = X16.5, ZT_18a = X18, ZT_19.5a = X19.5, ZT_21a = X21,
    ZT_22.5a = X22.5, 
    ZT_0b = X24, ZT_1.5b = X25.5, ZT_3b = X27, ZT_4.5b = X28.5,
    ZT_6b = X30, ZT_7.5b = X31.5, ZT_9b = X33, ZT_10.5b = X34.5, ZT_12b = X36,
    ZT_13.5b = X37.5, ZT_15b = X39, ZT_16.5b = X40.5, ZT_18b = X42, ZT_19.5b = X43.5,
    ZT_21b = X45, ZT_22.5b = X46.5
  )

# Réorganisation par timepoint avec réplicats adjacents
matrix_log_liver_rep <- matrix_log_liver_rep %>%
  select(
    ZT_0a, ZT_0b,
    ZT_1.5a, ZT_1.5b,
    ZT_3a, ZT_3b,
    ZT_4.5a, ZT_4.5b,
    ZT_6a, ZT_6b,
    ZT_7.5a, ZT_7.5b,
    ZT_9a, ZT_9b,
    ZT_10.5a, ZT_10.5b,
    ZT_12a, ZT_12b,
    ZT_13.5a, ZT_13.5b,
    ZT_15a, ZT_15b,
    ZT_16.5a, ZT_16.5b,
    ZT_18a, ZT_18b,
    ZT_19.5a, ZT_19.5b,
    ZT_21a, ZT_21b,
    ZT_22.5a, ZT_22.5b
  )

fwrite(matrix_log_liver_rep, 
       file = "matrix_log_liver_nodedup_with_replicates.txt", 
       sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# --- SOLÉAIRE ---
cat("Réorganisation des réplicats - Soléaire...\n")

matrix_log_soleus_rep <- matrix_log_soleus %>%
  rename(
    ZT_0a = X0, ZT_1.5a = X1.5, ZT_3a = X3, ZT_4.5a = X4.5, ZT_6a = X6,   
    ZT_7.5a = X7.5, ZT_9a = X9, ZT_10.5a = X10.5, ZT_12a = X12, ZT_13.5a = X13.5,
    ZT_15a = X15, ZT_16.5a = X16.5, ZT_18a = X18, ZT_19.5a = X19.5, ZT_21a = X21,
    ZT_22.5a = X22.5, 
    ZT_0b = X24, ZT_1.5b = X25.5, ZT_3b = X27, ZT_4.5b = X28.5,
    ZT_6b = X30, ZT_7.5b = X31.5, ZT_9b = X33, ZT_10.5b = X34.5, ZT_12b = X36,
    ZT_13.5b = X37.5, ZT_15b = X39, ZT_16.5b = X40.5, ZT_18b = X42, ZT_19.5b = X43.5,
    ZT_21b = X45, ZT_22.5b = X46.5
  )

matrix_log_soleus_rep <- matrix_log_soleus_rep %>%
  select(
    ZT_0a, ZT_0b,
    ZT_1.5a, ZT_1.5b,
    ZT_3a, ZT_3b,
    ZT_4.5a, ZT_4.5b,
    ZT_6a, ZT_6b,
    ZT_7.5a, ZT_7.5b,
    ZT_9a, ZT_9b,
    ZT_10.5a, ZT_10.5b,
    ZT_12a, ZT_12b,
    ZT_13.5a, ZT_13.5b,
    ZT_15a, ZT_15b,
    ZT_16.5a, ZT_16.5b,
    ZT_18a, ZT_18b,
    ZT_19.5a, ZT_19.5b,
    ZT_21a, ZT_21b,
    ZT_22.5a, ZT_22.5b
  )

fwrite(matrix_log_soleus_rep, 
       file = "matrix_log_soleus_nodedup_with_replicates.txt", 
       sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# --- GASTROCNÉMIEN ---
cat("Réorganisation des réplicats - Gastrocnémien...\n")

matrix_log_gastrocnemius_rep <- matrix_log_gastrocnemius %>%
  rename(
    ZT_0a = X0, ZT_1.5a = X1.5, ZT_3a = X3, ZT_4.5a = X4.5, ZT_6a = X6,   
    ZT_7.5a = X7.5, ZT_9a = X9, ZT_10.5a = X10.5, ZT_12a = X12, ZT_13.5a = X13.5,
    ZT_15a = X15, ZT_16.5a = X16.5, ZT_18a = X18, ZT_19.5a = X19.5, ZT_21a = X21,
    ZT_22.5a = X22.5, 
    ZT_0b = X24, ZT_1.5b = X25.5, ZT_3b = X27, ZT_4.5b = X28.5,
    ZT_6b = X30, ZT_7.5b = X31.5, ZT_9b = X33, ZT_10.5b = X34.5, ZT_12b = X36,
    ZT_13.5b = X37.5, ZT_15b = X39, ZT_16.5b = X40.5, ZT_18b = X42, ZT_19.5b = X43.5,
    ZT_21b = X45, ZT_22.5b = X46.5
  )

matrix_log_gastrocnemius_rep <- matrix_log_gastrocnemius_rep %>%
  select(
    ZT_0a, ZT_0b,
    ZT_1.5a, ZT_1.5b,
    ZT_3a, ZT_3b,
    ZT_4.5a, ZT_4.5b,
    ZT_6a, ZT_6b,
    ZT_7.5a, ZT_7.5b,
    ZT_9a, ZT_9b,
    ZT_10.5a, ZT_10.5b,
    ZT_12a, ZT_12b,
    ZT_13.5a, ZT_13.5b,
    ZT_15a, ZT_15b,
    ZT_16.5a, ZT_16.5b,
    ZT_18a, ZT_18b,
    ZT_19.5a, ZT_19.5b,
    ZT_21a, ZT_21b,
    ZT_22.5a, ZT_22.5b
  )

fwrite(matrix_log_gastrocnemius_rep, 
       file = "matrix_log_gastrocnemius_nodedup_with_replicates.txt", 
       sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

cat("Preprocessing terminé avec succès !\n")

# =============================================================================
# RÉSUMÉ DU PIPELINE
# =============================================================================
# 1. Chargement des matrices d'expression (Matrix Market format)
# 2. Conversion codes-barres -> identifiants d'échantillons
# 3. Conversion IDs Ensembl -> symboles de gènes
# 4. Filtrage basé sur les gènes de fond
# 5. Séparation par tissu (foie, soléaire, gastrocnémien)
# 6. Normalisation VST avec DESeq2 (globale et par tissu)
# 7. Réorganisation en réplicats biologiques pour analyse circadienne
# 
# Fichiers de sortie principaux :
# - matrix_log_liver_nodedup_with_replicates.txt
# - matrix_log_soleus_nodedup_with_replicates.txt  
# - matrix_log_gastrocnemius_nodedup_with_replicates.txt
# =============================================================================
