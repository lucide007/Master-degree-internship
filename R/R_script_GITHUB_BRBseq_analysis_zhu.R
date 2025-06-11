# R script re analysis Zhu data GSE130890 only wild type ####


setwd("path/to/zhu_data")

# Preprocessing ####

# Load raw data
gene_counts <- read.delim("gene_counts.tsv", comment.char="#")

matrix_zhu <- subset(gene_counts, select = -c(Chr, Start, End, Strand, Length))

# Convert Emsembl IDs to gene symbols
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
ensembl_ids <- matrix_zhu$Geneid

results <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)

matrix_zhu_annotated <- merge(matrix_zhu, results,
                              by.x = "Geneid", by.y = "ensembl_gene_id",
                              all.x = TRUE)
# Exemple : placer mgi_symbol juste après Geneid
matrix_zhu_annotated <- matrix_zhu_annotated[, c("Geneid", "mgi_symbol", setdiff(names(matrix_zhu_annotated), c("Geneid", "mgi_symbol")))]

matrix_zhu <- subset(matrix_zhu_annotated, select = -Geneid) # 43433

library(tidyverse)
matrix_zhu <- matrix_zhu %>%
  group_by(mgi_symbol) %>%
  sample_n(1) %>%
  ungroup() # 38974

# Apply same filter than our data
matrix_zhu <- matrix_zhu[apply(matrix_zhu, 1, max) > 117,] # 11629
matrix_zhu <- matrix_zhu[!grepl("ENSMUS|Rik|Gm", rownames(matrix_zhu)),] # 11629
matrix_zhu <- matrix_zhu[!is.na(matrix_zhu$mgi_symbol), ] 
matrix_zhu <- as.data.frame(matrix_zhu)

columns_zhu <- colnames(matrix_zhu)
columns_zhu

metadata <- read.csv("GSE130890_SRR_Acc_List.csv", stringsAsFactors = FALSE)

library(stringr)
srr_ids <- str_extract(columns_zhu, "SRR\\d+")

srr_gsm_map <- metadata[metadata$Run %in% srr_ids, c("Run", "SampleName")]

gsm_names <- srr_gsm_map$SampleName[match(srr_ids, srr_gsm_map$Run)]

new_colnames <- ifelse(is.na(gsm_names), columns_zhu, gsm_names)

sample_cols <- grepl("SRR\\d+", colnames(matrix_zhu))
colnames(matrix_zhu)[sample_cols] <- new_colnames[sample_cols]

# Save count matrix
fwrite(matrix_zhu, file="matrix_count_zhu.txt", sep="\t", quote=FALSE, row.names = TRUE, col.names = TRUE)

# Reload ####
matrix_zhu <- read.delim("matrix_count_zhu.txt")
rownames(matrix_zhu) <- matrix_zhu$X
matrix_zhu <- subset(matrix_zhu, select = -X)

# Generate Zt names
times <- seq(0, 46, by = 2)
labels <- as.vector(t(outer(times, c("a", "b"), function(z, r) paste0("ZT_", z, r))))

colnames(matrix_zhu) <- labels

# Normalize with vst ####
library(DESeq2)
se_all <- SummarizedExperiment(matrix_zhu)
meta <- as.data.frame(t(matrix_zhu))
dds_all <- DESeqDataSetFromMatrix(countData = matrix_zhu, colData = meta, design = ~ 1)
dds_all <- DESeq(dds_all)
vsd <- vst(dds_all, nsub = 11600)
matrix_log_zhu <- as.data.frame(assay(vsd))

# Save normalized matrix
fwrite(matrix_log_zhu, file="matrix_log_zhu.txt", sep="\t", quote=FALSE, row.names = TRUE, col.names = TRUE)

# Reload ####
matrix_log_zhu <- read.delim("matrix_log_zhu.txt")
rownames(matrix_log_zhu) <- matrix_log_zhu$X
matrix_log_zhu <- subset(matrix_log_zhu, select = -X)

# RAIN analysis ####
library(rain)

rain_zhu_12 <- rain(t(matrix_log_zhu), deltat = 2, period = 12, period.delta = 1.5, adjp.method = "BH", peak.border = c(0.3,0.7), nr.series = 2)
rain_zhu_sig_12 <- rain_zhu_12[(rain_zhu_12$pVal < 0.05 & rain_zhu_12$period %in% c(11, 12, 13)),] # 3447

rain_zhu_12["Emc9",]

rain_zhu_24 <- rain(t(matrix_log_zhu), deltat = 2, period = 24, period.delta = 2, adjp.method = "BH", peak.border = c(0.3,0.7), nr.series = 2)
rain_zhu_sig_24 <- rain_zhu_24[rain_zhu_24$pVal < 0.01,]

xbp1_zhu <- plot_gene_expression(gene_name = "Xbp1",
                                 expression_matrix = matrix_log_zhu,
                                 rain_df_list = list(rain_zhu_12,
                                                     rain_zhu_24),
                                 df_names = c("zhu 12h", "zhu 24h"),
                                 color = "navyblue")

list_zhu <- list(
  "zhu liver ~12h" = rownames(rain_zhu_sig_12),
  "zhu liver ~24h" = rownames(rain_zhu_sig_24)
)

library(VennDiagram)
venn_tissues_12 <- venn.diagram(x = list_zhu,
                                filename = NULL,  # Ne pas sauvegarder en fichier
                                col = "transparent",
                                fill = viridis::viridis(length(list_zhu), option = "C"),
                                alpha = 0.5,
                                label.col = "black",
                                cex = 2,
                                fontfamily = "sans",
                                cat.cex = 2,
                                cat.fontfamily = "sans",
                                cat.pos = c(-5, -6),  # Positions personnalisées pour les noms des ensembles
                                cat.dist = c(-0.3, -0.43))  # Distance des noms par rapport aux cercles

grid::grid.newpage()
grid::grid.draw(venn_tissues_12)

list_zhu_liver <- list(
  "zhu liver ~12h" = rownames(rain_zhu_sig_12),
  "zhu liver ~24h" = rownames(rain_zhu_sig_24),
  "team liver ~12h" = rownames(rain_liver_rep_sig_12),
  "team liver ~24h" = rownames(rain_liver_rep_sig_24)
)

rhythmic_genes <- unique(unlist(list_zhu_liver))
matrix_rhythmic_genes <- data.frame(gene = rhythmic_genes)

for (set_name in names(list_zhu_liver)) {
  matrix_rhythmic_genes[[set_name]] <- ifelse(matrix_rhythmic_genes$gene %in% list_zhu_liver[[set_name]], 1, 0)
}
rownames(matrix_rhythmic_genes) <- matrix_rhythmic_genes$gene
matrix_rhythmic_genes$gene <- NULL

upset(matrix_rhythmic_genes,
      sets = names(list_zhu_liver),
      order.by = "freq",
      point.size = 5,  
      line.size = 1.2,
      matrix.dot.alpha = 0.8,
      text.scale = 2,
      nintersects = NA,
      set_size.show = TRUE)

test <- list("zhu" = rownames(rain_zhu_sig_12),
             "liver pval < 0.1" = rownames(rain_liver_rep_triche_12))
venn_tissues_12 <- venn.diagram(x = test,
                                filename = NULL,  # Ne pas sauvegarder en fichier
                                col = "transparent",
                                fill = viridis::viridis(length(test), option = "C"),
                                alpha = 0.5,
                                label.col = "black",
                                cex = 2,
                                fontfamily = "sans",
                                cat.cex = 2,
                                cat.fontfamily = "sans",
                                cat.pos = c(-5, -6),  # Positions personnalisées pour les noms des ensembles
                                cat.dist = c(-0.3, -0.43))  # Distance des noms par rapport aux cercles

grid::grid.newpage()
grid::grid.draw(venn_tissues_12)

# Compare with published results ####
# with different peak_borders

# 12h genes
rain_zhu <- read_excel("path/to/plos_zhu_2020_rain_genes.xlsx", sheet = "Sheet 1")
list_zhu_12 <- readLines("path/to/list_zhu_rain_12h.txt")
rain_zhu_12 <- rain_zhu[rain_zhu$...1 %in% list_zhu_12,]
rain_zhu_12 <- rain_zhu_12 %>%
  group_by(`...1`) %>%                  # groupe par le nom du gène si doublons
  slice_min(order_by = adjP, n = 1, with_ties = FALSE) %>% # garde la ligne avec la plus petite adjpVal
  ungroup()
rain_zhu_12 <- as.data.frame(rain_zhu_12)
rownames(rain_zhu_12) <- rain_zhu_12$...1
rain_zhu_12 <- subset(rain_zhu_12, select = -...1)
rain_zhu_12_sig <- rain_zhu_12[!grepl("rik", rownames(rain_zhu_12), ignore.case = TRUE), ]

# 24h genes
rain_zhu_24 <- read_excel("plos_zhu_2020_rain_genes.xlsx", sheet = "Sheet 3")
list_zhu_24 <- readLines("path/to/list_zhu_rain_24h.txt")
rain_zhu_24 <- rain_zhu_24[rain_zhu_24$...1 %in% list_zhu_24,]
rain_zhu_24 <- rain_zhu_24 %>%
  group_by(`...1`) %>%                  # groupe par le nom du gène si doublons
  slice_min(order_by = pVal, n = 1, with_ties = FALSE) %>% # garde la ligne avec la plus petite adjpVal
  ungroup()
rain_zhu_24 <- as.data.frame(rain_zhu_24)
rownames(rain_zhu_24) <- rain_zhu_24$...1
rain_zhu_24 <- subset(rain_zhu_24, select = -...1)
rain_zhu_24_sig <- rain_zhu_24[!grepl("rik", rownames(rain_zhu_24), ignore.case = TRUE), ]

test2 <- list("zhu_my_RAIN" = rownames(rain_zhu_sig_12),
             "zhu_published" = list_zhu_12)
venn_tissues_12 <- venn.diagram(x = test2,
                                filename = NULL,  # Ne pas sauvegarder en fichier
                                col = "transparent",
                                fill = viridis::viridis(length(test2), option = "C"),
                                alpha = 0.5,
                                label.col = "black",
                                cex = 2,
                                fontfamily = "sans",
                                cat.cex = 2,
                                cat.fontfamily = "sans",
                                cat.pos = c(-5, -6),  # Positions personnalisées pour les noms des ensembles
                                cat.dist = c(-0.3, -0.43))  # Distance des noms par rapport aux cercles

grid::grid.newpage()
grid::grid.draw(venn_tissues_12)

rain_zhu_12_peakborder28 <- rain(t(matrix_log_zhu), deltat = 2, period = 12, period.delta = 1.5, adjp.method = "BH", peak.border = c(0.2,0.8), nr.series = 2)
rain_zhu_sig_12 <- rain_zhu_12_peakborder28[(rain_zhu_12_peakborder28$pVal < 0.05 & rain_zhu_12_peakborder28$period %in% c(11, 12, 13)),] # 3447
test2 <- list("zhu_my_RAIN" = rownames(rain_zhu_sig_12),
              "zhu_published" = list_zhu_12)
venn_tissues_12 <- venn.diagram(x = test2,
                                filename = NULL,  # Ne pas sauvegarder en fichier
                                col = "transparent",
                                fill = viridis::viridis(length(test2), option = "C"),
                                alpha = 0.5,
                                label.col = "black",
                                cex = 2,
                                fontfamily = "sans",
                                cat.cex = 2,
                                cat.fontfamily = "sans",
                                cat.pos = c(-5, -6),  # Positions personnalisées pour les noms des ensembles
                                cat.dist = c(-0.3, -0.43))  # Distance des noms par rapport aux cercles

grid::grid.newpage()
grid::grid.draw(venn_tissues_12)

rain_zhu_12_peakborder19 <- rain(t(matrix_log_zhu), deltat = 2, period = 12, period.delta = 1.5, adjp.method = "BH", peak.border = c(0.1,0.9), nr.series = 2)
rain_zhu_sig_12 <- rain_zhu_12_peakborder19[(rain_zhu_12_peakborder19$pVal < 0.05 & rain_zhu_12_peakborder19$period %in% c(11, 12, 13)),] # 3447
test2 <- list("zhu_my_RAIN" = rownames(rain_zhu_sig_12),
              "zhu_published" = list_zhu_12)
venn_tissues_12 <- venn.diagram(x = test2,
                                filename = NULL,  # Ne pas sauvegarder en fichier
                                col = "transparent",
                                fill = viridis::viridis(length(test2), option = "C"),
                                alpha = 0.5,
                                label.col = "black",
                                cex = 2,
                                fontfamily = "sans",
                                cat.cex = 2,
                                cat.fontfamily = "sans",
                                cat.pos = c(-5, -6),  # Positions personnalisées pour les noms des ensembles
                                cat.dist = c(-0.3, -0.43))  # Distance des noms par rapport aux cercles

grid::grid.newpage()
grid::grid.draw(venn_tissues_12)

setdiff(rownames(rain_zhu_12_peakborder28), rownames(rain_zhu_12_peakborder19))
# 0 diff donc peakborder 0.2-0.8 plafond

# Test on 24h ####
# Do we have the same results of we analyse the first 24h, the second 24h and the total 48h ?
matrix_log_zhu1 <- matrix_log_zhu[, 1:24]
matrix_log_zhu2 <- matrix_log_zhu[, 25:48]

rain_zhu_12.1 <- rain(t(matrix_log_zhu1), deltat = 2, period = 12, period.delta = 1.5, adjp.method = "BH", peak.border = c(0.3,0.7), nr.series = 2)
rain_zhu_sig_12.1 <- rain_zhu_12.1[(rain_zhu_12.1$pVal < 0.05 & rain_zhu_12.1$period %in% c(11, 12, 13)),] # 3447

rain_zhu_12.1["Noc3",]

rain_zhu_24.1 <- rain(t(matrix_log_zhu1), deltat = 2, period = 24, period.delta = 2, adjp.method = "BH", peak.border = c(0.3,0.7), nr.series = 2)
rain_zhu_sig_24.1 <- rain_zhu_24.1[rain_zhu_24.1$pVal < 0.05,]

rain_zhu_12.2 <- rain(t(matrix_log_zhu2), deltat = 2, period = 12, period.delta = 1.5, adjp.method = "BH", peak.border = c(0.3,0.7), nr.series = 2)
rain_zhu_sig_12.2 <- rain_zhu_12.2[(rain_zhu_12.2$pVal < 0.05 & rain_zhu_12.2$period %in% c(11, 12, 13)),] # 3447

rain_zhu_12.2["Ncoa3",]

rain_zhu_24.2 <- rain(t(matrix_log_zhu2), deltat = 2, period = 24, period.delta = 2, adjp.method = "BH", peak.border = c(0.3,0.7), nr.series = 2)
rain_zhu_sig_24.2 <- rain_zhu_24.2[rain_zhu_24.2$pVal < 0.05,]


test3 <- list("zhu 12h (48h)" = rownames(rain_zhu_sig_12),
              "zhu 12h.1 (24.1h)" = rownames(rain_zhu_sig_12.1),
              "zhu 12h.2 (24.2h)" = rownames(rain_zhu_sig_12.2))
venn_tissues_12 <- venn.diagram(x = test3,
                                filename = NULL,  # Ne pas sauvegarder en fichier
                                col = "transparent",
                                fill = viridis::viridis(length(test3), option = "C"),
                                alpha = 0.5,
                                label.col = "black",
                                cex = 2,
                                fontfamily = "sans",
                                cat.cex = 2,
                                cat.fontfamily = "sans")
                                # cat.pos = c(-5, -6),  # Positions personnalisées pour les noms des ensembles
                                # cat.dist = c(-0.3, -0.43))  # Distance des noms par rapport aux cercles

grid::grid.newpage()
grid::grid.draw(venn_tissues_12)


rain_zhu_real_pval <- rain_zhu[rain_zhu$pVal < 0.05, ] # 5317
test4 <- list("zhu RAIN pval" = rain_zhu_real_pval$...1,
              "my zhu" = rownames(rain_zhu_sig_12))
venn_tissues_12 <- venn.diagram(x = test4,
                                filename = NULL,  # Ne pas sauvegarder en fichier
                                col = "transparent",
                                fill = viridis::viridis(length(test4), option = "C"),
                                alpha = 0.5,
                                label.col = "black",
                                cex = 2,
                                fontfamily = "sans",
                                cat.cex = 2,
                                cat.fontfamily = "sans")
# cat.pos = c(-5, -6),  # Positions personnalisées pour les noms des ensembles
# cat.dist = c(-0.3, -0.43))  # Distance des noms par rapport aux cercles

grid::grid.newpage()
grid::grid.draw(venn_tissues_12)


test5 <- list("zhu 12h (48h)" = rownames(rain_zhu_sig_12),
              "zhu 12h.1 (24.1h)" = rownames(rain_zhu_sig_12.1),
              "zhu 12h.2 (24.2h)" = rownames(rain_zhu_sig_12.2),
              "team's 12h liver" = rownames(rain_liver_rep_sig_12))
venn_tissues_12 <- venn.diagram(x = test5,
                                filename = NULL,  # Ne pas sauvegarder en fichier
                                col = "transparent",
                                fill = viridis::viridis(length(test5), option = "C"),
                                alpha = 0.5,
                                label.col = "black",
                                cex = 2,
                                fontfamily = "sans",
                                cat.cex = 2,
                                cat.fontfamily = "sans")
# cat.pos = c(-5, -6),  # Positions personnalisées pour les noms des ensembles
# cat.dist = c(-0.3, -0.43))  # Distance des noms par rapport aux cercles

grid::grid.newpage()
grid::grid.draw(venn_tissues_12)

# Answer is : no

# liste de genes publiés ####

genes_publi_2017_2020 <- c("Sec23b", "Hspa5", "Creld2", "Gck", "Acly", "Cpt1a", "Elovl6", "Tbp", "Gtf2h3", "Gtf2b", "Eif3a", "Srp72", "Txndc5", "Alg12", "Cnot3")

# Afficher les lignes correspondantes
rain_zhu_sig_12.1[genes_publi_2017_2020, ] # 3 / 15
rain_zhu_sig_12.2[genes_publi_2017_2020, ] # 4 / 15
rain_zhu_sig_12[genes_publi_2017_2020, ] # 10 / 15
rain_liver_rep_sig_12[genes_publi_2017_2020,] # 2 / 15
rain_soleus_rep_sig_12[genes_publi_2017_2020,] # 1/15
rain_gastrocnemius_rep_sig_12[genes_publi_2017_2020, ] # 1/15
genes_publi_2017_2020

# avec metacycle sur 48h, 24h.1 et 24h.2 VS liver data ####
library(MetaCycle)

time48 <- seq(0, 46, by = 2)
time48 <- rep(time48, each = 2)
time24 <- seq(0, 22, by = 2)
time24 <- rep(time24, each = 2)

# 48h
df_zhu48 <- matrix_log_zhu
df_zhu48$gene <- rownames(df_zhu48)
df_zhu48 <- df_zhu48[, c("gene", setdiff(names(df_zhu48), "gene"))]
write.csv(df_zhu48, file="matrix_log_zhu48_wgenes.csv", row.names = FALSE)

meta2d_zhu48 <- meta2d(infile = "matrix_log_zhu48_wgenes.csv", filestyle = "csv",
                       outdir = "metacycle_zhu48", timepoints = time48,
                       minper = 10, maxper = 14, combinePvalue = "bonferroni",
                       cycMethod = c("JTK", "LS", "ARS"),
                       ARSdefaultPer = 12)

# 24h.1
df_zhu24.1 <- matrix_log_zhu1
df_zhu24.1$gene <- rownames(df_zhu24.1)
df_zhu24.1 <- df_zhu24.1[, c("gene", setdiff(names(df_zhu24.1), "gene"))]
write.csv(df_zhu24.1, file="matrix_log_zhu24.1_wgenes.csv", row.names = FALSE)
meta2d_zhu24.1 <- meta2d(infile = "matrix_log_zhu24.1_wgenes.csv", filestyle = "csv",
                       outdir = "metacycle_zhu24.1", timepoints = time24,
                       minper = 10, maxper = 14, combinePvalue = "bonferroni",
                       cycMethod = c("JTK", "LS", "ARS"),
                       ARSdefaultPer = 12)

# 24h.2
df_zhu24.2 <- matrix_log_zhu2
df_zhu24.2$gene <- rownames(df_zhu24.2)
df_zhu24.2 <- df_zhu24.2[, c("gene", setdiff(names(df_zhu24.2), "gene"))]
write.csv(df_zhu24.2, file="matrix_log_zhu24.2_wgenes.csv", row.names = FALSE)
meta2d_zhu24.2 <- meta2d(infile = "matrix_log_zhu24.2_wgenes.csv", filestyle = "csv",
                         outdir = "metacycle_zhu24.2", timepoints = time24,
                         minper = 10, maxper = 14, combinePvalue = "bonferroni",
                         cycMethod = c("JTK", "LS", "ARS"),
                         ARSdefaultPer = 12)



meta2d_zhu48_results <- read.table("metacycle_zhu48/meta2d_matrix_log_zhu48_wgenes.csv", header = T, sep = ",")
meta2d_zhu24.1_results <- read.table("metacycle_zhu24.1/meta2d_matrix_log_zhu24.1_wgenes.csv", header = T, sep = ",")
meta2d_zhu24.2_results <- read.table("metacycle_zhu24.2/meta2d_matrix_log_zhu24.2_wgenes.csv", header = T, sep = ",")

meta2d_liver_results <- read.table("path/to/metacycle_liver_results_rep/meta2d_matrix_log_liver_rep_wgenes.csv", header = T, sep = ",")

meta2d_zhu48_sig <- meta2d_zhu48_results[meta2d_zhu48_results$meta2d_pvalue < 0.05,]

list_meta2d <- list("meta2d zhu on 48h" = meta2d_zhu48_results$CycID[meta2d_zhu48_results$meta2d_pvalue < 0.05],
                    "meta2d zhu on 24h.1" = meta2d_zhu24.1_results$CycID[meta2d_zhu24.1_results$meta2d_pvalue < 0.05],
                    "meta2d zhu on 24h.2" = meta2d_zhu24.2_results$CycID[meta2d_zhu24.2_results$meta2d_pvalue < 0.05],
                    "meta2d on liver data" = meta2d_liver_results$CycID[meta2d_liver_results$meta2d_pvalue < 0.05])
venn_meta2d <- venn.diagram(x = list_meta2d,
                                filename = NULL,  # Ne pas sauvegarder en fichier
                                col = "transparent",
                                fill = viridis::viridis(length(list_meta2d), option = "C"),
                                alpha = 0.5,
                                label.col = "black",
                                cex = 2,
                                fontfamily = "sans",
                                cat.cex = 2,
                                cat.fontfamily = "sans")
# cat.pos = c(-5, -6),  # Positions personnalisées pour les noms des ensembles
# cat.dist = c(-0.3, -0.43))  # Distance des noms par rapport aux cercles

grid::grid.newpage()
grid::grid.draw(venn_meta2d)

# check genes significatifs dans 48h qui ne sont pas dans les 24h.1 et 24h.2

setdiff(rownames(rain_zhu_sig_12.1), rownames(rain_zhu_sig_12.2))

# Enrichment analysis ####
library(enrichR)
library(gridExtra)
# + function plot_enrichment()

dbs <- c("GO_Biological_Process_2025",
         "GO_Cellular_Component_2025",
         "GO_Molecular_Function_2025",
         "KEGG_2019_Mouse")

enrichr_zhu12 <- plot_enrichment(rownames(rain_zhu_sig_12), dbs, max_cols = 1)
enrichr_zhu12.1 <- plot_enrichment(rownames(rain_zhu_sig_12.1), dbs, max_cols = 1)
enrichr_zhu12.2 <- plot_enrichment(rownames(rain_zhu_sig_12.2), dbs, max_cols = 1)

enrichr_liver <- plot_enrichment(rownames(rain_liver_rep_sig_12), dbs, max_cols = 1)
















