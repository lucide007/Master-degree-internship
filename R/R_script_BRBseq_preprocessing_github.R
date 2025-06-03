# R script for BRB data preprocessing #

# packages ####

library(data.table)
library(Matrix)
library(dplyr)
library(tidyverse)
library(DESeq2)
library(SummarizedExperiment)

# files ####
# Répertoire contenant les fichiers Matrix Market
matrix_dir <- "path/to/repertory"

# Chargement des noms de gènes et des codes barres
features <- fread(file.path(matrix_dir, "features.tsv"), header = FALSE, data.table = FALSE)
barcodes <- fread(file.path(matrix_dir, "barcodes.tsv"), header = FALSE, data.table = FALSE)

# Chargement des annotations de barcodes (échantillons)
barcode_annot <- read.delim("path/to/barcodes_row.txt")

# matrix no deduplication
matrix1 <- readMM(file.path(matrix_dir, "umiDedup-NoDedup.mtx")) %>%
  as.matrix() %>%
  as.data.frame()

colnames(matrix1) <- barcodes$V1
rownames(matrix1) <- features$V1
fwrite(matrix1, file = "raw_matrix_umi_no_dedup.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# matrix 1 mismatch
matrix2 <- readMM(file.path(matrix_dir, "umiDedup-1MM_All.mtx")) %>%
  as.matrix() %>%
  as.data.frame()

colnames(matrix2) <- barcodes$V1
rownames(matrix2) <- features$V1
fwrite(matrix2, file = "raw_matrix_umi_1_MM.txt", sep = "\t", quote = FALSE, row.names = TRUE)


# conversion barcodes sequences to sample id ####
matrix1_named <- matrix1 %>%
  rename_with(~ barcode_annot$sample_cycle[match(., barcode_annot$barcode)])

# conversion ensembl id to gene symbol ####
genes <- features[, c(1, 2)] # V1: Ensembl ID, V2: symbole de gène
colnames(genes) <- c("ensembl_id", "symbol")

matrix1$ensembl_id <- rownames(matrix1)
matrix1 <- merge(matrix1, genes, by = "ensembl_id", all.x = TRUE)
rownames(matrix1) <- matrix1$symbol


# filtering background genes ####
bg_genes <- scan("path/to/background_genes.txt", what = "", quiet = TRUE)

# Extraction de la matrice de fond
bg_matrix <- matrix_count[rownames(matrix_count) %in% bg_genes,]
seuil_max_bg <- max(bg_matrix)

# Garder les gènes exprimés plus que le max des gènes de fond
matrix_count <- matrix_count[apply(matrix_count, 1, max) > seuil_max_bg,]
fwrite(matrix_count, file = "matrix_count_nodedup.txt", sep = "\t", quote = FALSE)

# Faire une matrice pour chaque tissu
# Charger fichier
matrix_count <- read.delim("path/to/matrix_count_nodedup.txt")
rownames(matrix_count) <- matrix_count$X
matrix_count <- subset(matrix_count, select = -X)

colnames(matrix_count) <- gsub(",", ".", colnames(matrix_count))

ordre_voulu <- c("liv_0", "liv_1.5", "liv_3", "liv_4.5", "liv_6", "liv_7.5", "liv_9", "liv_10.5", "liv_12", "liv_13.5", "liv_15", 'liv_16.5', "liv_18", "liv_19.5", "liv_21", "liv_22.5", "liv_24", "liv_25.5", "liv_27", "liv_28.5", "liv_30", "liv_31.5", "liv_33", "liv_34.5", "liv_36", "liv_37.5", "liv_39", "liv_40.5", "liv_42", "liv_43.5", "liv_45", "liv_46.5",
                 "sol_0", "sol_1.5", "sol_3", "sol_4.5", "sol_6", "sol_7.5", "sol_9", "sol_10.5", "sol_12", "sol_13.5", "sol_15", 'sol_16.5', "sol_18", "sol_19.5", "sol_21", "sol_22.5", "sol_24", "sol_25.5", "sol_27", "sol_28.5", "sol_30", "sol_31.5", "sol_33", "sol_34.5", "sol_36", "sol_37.5", "sol_39", "sol_40.5", "sol_42", "sol_43.5", "sol_45", "sol_46.5",
                 "gast_0", "gast_1.5", "gast_3", "gast_4.5", "gast_6", "gast_7.5", "gast_9", "gast_10.5", "gast_12", "gast_13.5", "gast_15", 'gast_16.5', "gast_18", "gast_19.5", "gast_21", "gast_22.5", "gast_24", "gast_25.5", "gast_27", "gast_28.5", "gast_30", "gast_31.5", "gast_33", "gast_34.5", "gast_36", "gast_37.5", "gast_39", "gast_40.5", "gast_42", "gast_43.5", "gast_45", "gast_46.5")

matrix_count <- matrix_count[, ordre_voulu]

# subsets
matrix_count_liver <- matrix_count[, grepl("^liv", colnames(matrix_count))]
matrix_count_soleus <- matrix_count[, grepl("^sol", colnames(matrix_count))]
matrix_count_gastrocnemius <- matrix_count[, grepl("^gast", colnames(matrix_count))]

fwrite(matrix_count_liver, file="matrix_count_nodedup_liver.txt", sep="\t", quote=FALSE, row.names = TRUE, col.names = TRUE)
fwrite(matrix_count_soleus, file="matrix_count_nodedup_soleus.txt", sep="\t", quote=FALSE, row.names = TRUE, col.names = TRUE)
fwrite(matrix_count_gastrocnemius, file="matrix_count_nodedup_gastrocnemius.txt", sep="\t", quote=FALSE, row.names = TRUE, col.names = TRUE)


# normalization log2 all tissues ####
# Chargement
matrix_count <- fread("matrix_count_nodedup.txt", data.table = FALSE)
rownames(matrix_count) <- matrix_count$V1
matrix_count <- matrix_count[, -1]
colnames(matrix_count) <- gsub(",", ".", colnames(matrix_count))

# Création de l'objet DESeq2
meta <- as.data.frame(t(matrix_count))
dds <- DESeqDataSetFromMatrix(countData = matrix_count, colData = meta, design = ~1)
dds <- DESeq(dds)

# Normalisation
matrix_log <- vst(dds, nsub = 10600) %>% assay() %>% as.data.frame()
fwrite(matrix_log, file = "matrix_log_nodedup.txt", sep = "\t", quote = FALSE)

# per tissu normalization - liver ####
matrix_count_liver <- fread("matrix_count_nodedup_liver.txt", data.table = FALSE)
rownames(matrix_count_liver) <- matrix_count_liver$V1
matrix_count_liver <- matrix_count_liver[, -1]
colnames(matrix_count_liver) <- gsub(",", ".", colnames(matrix_count_liver))

sel <- SummarizedExperiment(matrix_count_liver)
metal <- as.data.frame(t(matrix_count_liver))
ddsl <- DESeqDataSetFromMatrix(countData = matrix_count_liver, colData = metal, design = ~ 1)
vsdl <- vst(ddsl, nsub = 10000)
matrix_log_liver <- as.data.frame(assay(vsdl))
colnames(matrix_log_liver) <- gsub("^liv_", "", colnames(matrix_log_liver))
fwrite(matrix_log_liver, file="matrix_log_liver_nodedup.txt", sep="\t", quote=FALSE, row.names = TRUE, col.names = TRUE)


# per tissu normalization - soleus ####
matrix_count_soleus <- fread("matrix_count_nodedup_soleus.txt", data.table = FALSE)
rownames(matrix_count_soleus) <- matrix_count_soleus$V1
matrix_count_soleus <- matrix_count_soleus[, -1]
colnames(matrix_count_soleus) <- gsub(",", ".", colnames(matrix_count_soleus))

sel <- SummarizedExperiment(matrix_count_soleus)
metal <- as.data.frame(t(matrix_count_soleus))
ddsl <- DESeqDataSetFromMatrix(countData = matrix_count_soleus, colData = metal, design = ~ 1)
vsdl <- vst(ddsl, nsub = 10000)
matrix_log_soleus <- as.data.frame(assay(vsdl))
colnames(matrix_log_soleus) <- gsub("^liv_", "", colnames(matrix_log_soleus))
fwrite(matrix_log_soleus, file="matrix_log_soleus_nodedup.txt", sep="\t", quote=FALSE, row.names = TRUE, col.names = TRUE)

# per tissu normalization - gastrocnemius ####
matrix_count_gastrocnemius <- fread("matrix_count_nodedup_gastrocnemius.txt", data.table = FALSE)
rownames(matrix_count_gastrocnemius) <- matrix_count_gastrocnemius$V1
matrix_count_gastrocnemius <- matrix_count_gastrocnemius[, -1]
colnames(matrix_count_gastrocnemius) <- gsub(",", ".", colnames(matrix_count_gastrocnemius))

sel <- SummarizedExperiment(matrix_count_gastrocnemius)
metal <- as.data.frame(t(matrix_count_gastrocnemius))
ddsl <- DESeqDataSetFromMatrix(countData = matrix_count_gastrocnemius, colData = metal, design = ~ 1)
vsdl <- vst(ddsl, nsub = 10000)
matrix_log_gastrocnemius <- as.data.frame(assay(vsdl))
colnames(matrix_log_gastrocnemius) <- gsub("^liv_", "", colnames(matrix_log_gastrocnemius))
fwrite(matrix_log_gastrocnemius, file="matrix_log_gastrocnemius_nodedup.txt", sep="\t", quote=FALSE, row.names = TRUE, col.names = TRUE)

# mettre les échantillons dans l'ordre des replicas - liver ####
matrix_log_liver_rep <- matrix_log_liver %>%
  rename(
    ZT_0a = X0, ZT_1.5a = X1.5, ZT_3a = X3, ZT_4.5a = X4.5, ZT_6a = X6,   
    ZT_7.5a = X7.5, ZT_9a = X9, ZT_10.5a = X10.5, ZT_12a = X12, ZT_13.5a = X13.5,
    ZT_15a = X15, ZT_16.5a = X16.5, ZT_18a = X18, ZT_19.5a = X19.5, ZT_21a = X21,
    ZT_22.5a = X22.5, ZT_0b = X24, ZT_1.5b = X25.5, ZT_3b = X27, ZT_4.5b = X28.5,
    ZT_6b = X30, ZT_7.5b = X31.5, ZT_9b = X33, ZT_10.5b = X34.5, ZT_12b = X36,
    ZT_13.5b = X37.5, ZT_15b = X39, ZT_16.5b = X40.5, ZT_18b = X42, ZT_19.5b = X43.5,
    ZT_21b = X45, ZT_22.5b = X46.5
  )
    
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
fwrite(matrix_log_liver_rep, file="matrix_log_liver_nodedup_with_replicates.txt", sep="\t", quote=FALSE, row.names = TRUE, col.names = TRUE)

# mettre les échantillons dans l'ordre des replicas - soleus ####
matrix_log_soleus_rep <- matrix_log_soleus %>%
  rename(
    ZT_0a = X0, ZT_1.5a = X1.5, ZT_3a = X3, ZT_4.5a = X4.5, ZT_6a = X6,   
    ZT_7.5a = X7.5, ZT_9a = X9, ZT_10.5a = X10.5, ZT_12a = X12, ZT_13.5a = X13.5,
    ZT_15a = X15, ZT_16.5a = X16.5, ZT_18a = X18, ZT_19.5a = X19.5, ZT_21a = X21,
    ZT_22.5a = X22.5, ZT_0b = X24, ZT_1.5b = X25.5, ZT_3b = X27, ZT_4.5b = X28.5,
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
fwrite(matrix_log_soleus_rep, file="matrix_log_soleus_nodedup_with_replicates.txt", sep="\t", quote=FALSE, row.names = TRUE, col.names = TRUE)

# mettre les échantillons dans l'ordre des replicas - gastrocnemius ####
matrix_log_gastrocnemius_rep <- matrix_log_gastrocnemius %>%
  rename(
    ZT_0a = X0, ZT_1.5a = X1.5, ZT_3a = X3, ZT_4.5a = X4.5, ZT_6a = X6,   
    ZT_7.5a = X7.5, ZT_9a = X9, ZT_10.5a = X10.5, ZT_12a = X12, ZT_13.5a = X13.5,
    ZT_15a = X15, ZT_16.5a = X16.5, ZT_18a = X18, ZT_19.5a = X19.5, ZT_21a = X21,
    ZT_22.5a = X22.5, ZT_0b = X24, ZT_1.5b = X25.5, ZT_3b = X27, ZT_4.5b = X28.5,
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
fwrite(matrix_log_gastrocnemius_rep, file="matrix_log_gastrocnemius_nodedup_with_replicates.txt", sep="\t", quote=FALSE, row.names = TRUE, col.names = TRUE)

