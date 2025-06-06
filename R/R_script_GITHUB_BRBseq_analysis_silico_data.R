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
  cycMethod = c("JTK", "LS", "ARS"), ARSdefaultPer = 720
)

# Compare rhythmic sets
meta2d_silico_results <- read.table("metacycle_liver_results/meta2d_silico_data_for_metacycle.csv", header = TRUE, sep = ",")
gene_silico_sets <- list(
  "Methode RAIN" = rain_silico_results,
  "Methode Lomb-Scargle" = meta2d_silico_results$LS_pvalue < 0.05,
  "Methode JTK" = meta2d_silico_results$JTK_pvalue < 0.05,
  "Methode ARSER" = meta2d_silico_results$ARS_pvalue < 0.05,
  "Methode Metacycle 2d" = meta2d_silico_results$meta2d_pvalue < 0.05
)

# Create a matrix of selected rhythmic genes
matrix_12_silico <- unique(unlist(gene_silico_sets)) %>%
  as.data.frame() %>%
  rename(gene = value) %>%
  mutate(across(everything(), ~ as.integer(gene %in% .)))

# Save selected genes
fwrite(matrix_12_silico, file = "matrix_silico_12h_rhythmic_genes2.csv", sep = ";", row.names = TRUE, col.names = TRUE)

# Plot results
upsetplot_silico <- upset(matrix_12_silico, 
                          sets = names(gene_silico_sets), 
                          order.by = "freq",
                          mainbar.y.label = "12h rhythmic genes in methods intersections",
                          sets.x.label = "Set sizes (in silico)",
                          point.size = 5,  
                          line.size = 1.2,
                          matrix.dot.alpha = 0.8,
                          text.scale = 2)

# Phase correlation
matrix_phase_silico <- full_join(matrix_phase_silico_rain, matrix_phase_silico_jtk, by = "CycID") %>%
  full_join(matrix_phase_silico_ars, by = "CycID") %>%
  full_join(matrix_phase_silico_ls, by = "CycID") %>%
  full_join(matrix_phase_silico_meta2d, by = "CycID")

# Plot phase correlations
ggpairs(matrix_phase_silico[, -1], 
        title = "Phases correlations : silico, 12h rhythmic genes",
        upper = list(continuous = wrap("cor", size = 6)),
        lower = list(continuous = wrap("points", alpha = 0.8, size = 1.5, color = "royalblue")),
        diag = list(continuous = wrap("densityDiag", size = 0.8))) +
  theme_bw() +
  theme(axis.text = element_text(size = 12))

# Clear memory
rm(list = ls(pattern = "silico"))
