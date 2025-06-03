# Benchmarking methods on silico data ####

# packages

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
# + load function circular_phase12H_histogram()

# circa in silico used to generate silico dataset

# 24h, 1 point / 1.5h, 2 replicat, 30% de features rythmiques avec periode 12h min amp 1 max amp 6 out amp 0
silico_data <- read.csv("path/to/CIS_data_2025-04-25.csv")

# 48h, 1 point / 1.5h, 1 replicat, 30% features rythmiques avec periode 10.5-13.5h min amp 1 max amp 3 out amp 0
silico_data <- read.csv("path/to/CIS_data_2025-04-17.csv")

silico_data2 <- subset(silico_data, select = -c(Cosine.Amplitude, Cosine.Periods, Cosine.Lag.Factors))
silico_data3 <- subset(silico_data2, select = -X)
silico_data3 <- na.omit(silico_data3)


# à partir du fichier CIS_data_2025-04-17.csv pour faire 2 replicates et période 10.5h-13.5h
colnames(silico_data3) <- sub("^X", "", colnames(silico_data3))
# Créer deux sous-ensembles pour les réplicats
rep1 <- silico_data3[, 1:16]
rep2 <- silico_data3[, 17:32]
# Renommer les colonnes de chaque sous-ensemble (1 à 24)
temps_seq <- seq(1.5, 24, by = 1.5)
colnames(silico_data3)[1:16] <- paste0("Rep1_T", temps_seq)
colnames(silico_data3)[17:32] <- paste0("Rep2_T", temps_seq)
rep1_cols <- grep("^Rep1", colnames(silico_data3), value = TRUE)
rep2_cols <- grep("^Rep2", colnames(silico_data3), value = TRUE)
# Créer l’ordre intercalé : Rep1_T1.5, Rep2_T1.5, Rep1_T3, Rep2_T3, ...
interleaved_cols <- as.vector(rbind(rep1_cols, rep2_cols))
# Réordonner le data frame
silico_data3 <- silico_data3[, interleaved_cols]
silico_data2 <- silico_data3
silico_data2$X <- rownames(silico_data2)
silico_data2 <- silico_data2[, c("X", setdiff(names(silico_data2), "X"))]


# dryR
silico_round <- round(silico_data3)
silico_round7 <- silico_round + 7
silico <- c(rep("silico", 32))
drylist_silico = dryseq_single_multiperiod(silico_round7, silico, time_minute, period_range = c(630, 810), period_step = c(6))

# rain
# data with 2 replicates
rain_silico <- rain(t(silico_data3), deltat = 1.5, period = 12, period.delta = 1.5, adjp.method = "BH", peak.border = c(0.3,0.7), nr.series = 2)
rain_silico <- rain(t(silico_data3), deltat = 1.5, period = 12, period.delta = 1.5, adjp.method = "BH", peak.border = c(0.3,0.7))
# data with 1 replicate
rain_silico <- rain(t(silico_data3), deltat = 1.5, period = 12, period.delta = 1.5, adjp.method = "BH", peak.border = c(0.3,0.7))

rain_silico_results <- rownames(rain_silico)[rain_silico$pVal < 0.05] # 396
rain_silico_results <- rownames(rain_silico)[rain_silico$pVal < 0.01] # 333 ! 
matrix_silico_rain <- silico_data3[rownames(silico_data3) %in% rain_silico_results,]

# metacycle
time <- c(1.5,3,4.5,6,7.5,9,10.5,12,13.5,15,16.5,18,19.5,21,22.5,24)
time <- rep(time, each = 2)
time_minute <- time * 60
write.csv(silico_data2, file="silico_data_for_metacycle.csv", row.names = FALSE)
meta2d_silico <- meta2d(
    infile = "silico_data_for_metacycle.csv", filestyle = "csv",
    outdir = "metacycle_liver_results", timepoints = time_minute,
    minper = 720, maxper = 720, combinePvalue = "bonferroni",                    
    cycMethod = c("JTK", "LS", "ARS"),                    
    ARSdefaultPer = 720
)                    

meta2d_silico_results <- read.table("metacycle_liver_results/meta2d_silico_data_for_metacycle.csv", header = T, sep = ",")

# comparing rhythmic sets ####
matrix_silico_meta2d <- subset(meta2d_silico_results, meta2d_silico_results$meta2d_pvalue < 0.05) # 295
matrix_silico_jtk <- subset(meta2d_silico_results, meta2d_silico_results$JTK_pvalue < 0.05) # 286
matrix_silico_ls <- subset(meta2d_silico_results, meta2d_silico_results$LS_pvalue < 0.05) # 250
matrix_silico_ars <- subset(meta2d_silico_results, meta2d_silico_results$ARS_pvalue < 0.05) # 304

#gene_silico_dryseq_single <- rownames(matrix_silico_dryr)
gene_silico_rain <- rain_silico_results
gene_silico_ls <- matrix_silico_ls$CycID
gene_silico_jtk <- matrix_silico_jtk$CycID
gene_silico_ars <- matrix_silico_ars$CycID
gene_silico_meta2d <- matrix_silico_meta2d$CycID

gene_silico_sets <- list(
    "Methode RAIN" = gene_silico_rain,
    "Methode Lomb-Scargle" = gene_silico_ls,
    "Methode JTK" = gene_silico_jtk,
    "Methode ARSER" = gene_silico_ars,
    "Methode Metacycle 2d" = gene_silico_meta2d)

# convert to a matrix
matrix_12_silico <- unique(unlist(gene_silico_sets)) %>% as.data.frame()
colnames(matrix_12_silico) <- "gene"
matrix_12_silico <- matrix_12_silico %>%
  mutate(rain = as.integer(gene %in% gene_silico_sets$`Methode RAIN`),
         ls = as.integer(gene %in% gene_silico_sets$`Methode Lomb-Scargle`),
         jtk = as.integer(gene %in% gene_silico_sets$`Methode JTK`),
         ars = as.integer(gene %in% gene_silico_sets$`Methode ARSER`),
         meta2d = as.integer(gene %in% gene_silico_sets$`Methode Metacycle 2d`))

# extract matrix of selected rhythmic genes
matrix_12.1 <- matrix_12_silico[rowSums(matrix_12_silico[, c("rain", "ls", "jtk", "ars", "meta2d")]) >= 2,]
matrix_12.2 <- matrix_12_silico[matrix_12_silico$rain == 1,]
matrix_12_silico_selected <- unique(rbind(matrix_12.1, matrix_12.2))
# save it
fwrite(matrix_12_silico_selected, file="matrix_silico_12h_rhythmic_genes2.csv", sep=";", row.names = TRUE, col.names = TRUE)

# plot the results
upsetplot_silico <- upset(
  matrix_12_silico, 
  sets = c("rain", "ls", "jtk", "ars", "meta2d"),  
  order.by = "freq",
  mainbar.y.label = "12h rhythmic genes in methods intersections",
  sets.x.label = "Set sizes (in silico)",
  point.size = 5,  
  line.size = 1.2,
  matrix.dot.alpha = 0.8,
  text.scale = 2)
upsetplot_silico




# phase correlation ####
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

# plot the results
ggpairs(matrix_phase_silico[, -1], 
        title = "Phases correlations : silico, 12h rhythmic genes",
        upper = list(continuous = wrap("cor", size = 6)),
        lower = list(continuous = wrap("points", alpha = 0.8, size = 1.5, color = "royalblue")),
        diag = list(continuous = wrap("densityDiag", size = 0.8)))+
  theme_bw()+
  theme(axis.text = element_text(size = 12))

hist_meta2d <- circular_phase12H_histogram(matrix_phase_silico_meta2d, "meta2d_phase", "CycID", title_override = "meta2D method", fill_color = "royalblue3")
hist_rain <- circular_phase12H_histogram(matrix_phase_silico_rain, "rain_phase", "CycID", title_override = "RAIN method", fill_color = "royalblue3")
hist_jtk <- circular_phase12H_histogram(matrix_phase_silico_jtk, "JTK_adjphase", "CycID", title_override = "JTK method", fill_color = "royalblue3")
hist_ars <- circular_phase12H_histogram(matrix_phase_silico_ars, "ARS_adjphase", "CycID", title_override = "ARS method", fill_color = "royalblue3")
hist_ls <- circular_phase12H_histogram(matrix_phase_silico_ls, "LS_adjphase", "CycID", title_override = "LS method", fill_color = "royalblue3")

hist_phases <- ( hist_meta2d + hist_rain ) / ( hist_jtk + hist_ls + hist_ars)+
  plot_annotation(title = "silico 12h rhythmic genes - Phase distribution by method",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 16)))
hist_phases

# vider la mémoire pour la suite !
rm(list = ls(pattern = "silico"))