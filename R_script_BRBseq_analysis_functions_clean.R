# ==========================================================================
# CIRCADIAN ANALYSIS FUNCTIONS
# ==========================================================================
# These functions facilitate the analysis of circadian gene expression data,
# including circular histograms, differential expression analysis,
# enrichment analysis, and gene expression plotting.
# Author: Original code adapted and cleaned
# ==========================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(DESeq2)
library(enrichR)
library(gridExtra)
library(grid)

# ==========================================================================
# CIRCULAR HISTOGRAM FUNCTIONS
# ==========================================================================

#' Create a circular histogram for phase data with 720-minute (12-hour) scale
#'
#' @param df_original Dataframe containing phase data
#' @param phase_col Column name containing phase values
#' @param gene_col Column name containing gene identifiers (optional)
#' @param ttot Total time period in minutes (default: 720)
#' @param title_override Optional title override
#' @param fill_color Color for histogram bars (default: "dodgerblue4")
#' @return A ggplot2 object with circular histogram
circular_phase720min_histogram <- function(df_original,
                                           phase_col,
                                           gene_col = NULL,
                                           ttot = 720, 
                                           title_override = NULL,
                                           fill_color = "dodgerblue4") {
  # Extract phase data
  x <- df_original[[phase_col]]
  
  # Normalize phase values to the interval [0, ttot]
  x <- x %% ttot
  
  # Create histogram
  br <- seq(0, ttot, length.out = 25)
  h <- hist(x, breaks = br, plot = FALSE)
  
  # Create dataframe for ggplot
  df <- data.frame(x = as.numeric(h$mids), y = h$counts)
  
  # Parameters for hour display
  hours_breaks <- seq(0, ttot, 60)
  hours_labels <- seq(0, 12, 1)
  
  # Set title
  title <- title_override
  
  # Subtitle with number of genes
  n_genes <- length(unique(df_original[[gene_col]]))
  subtitle <- paste(n_genes, "genes")
  
  # Create the plot
  p <- ggplot(df, aes(x=x, y=y)) + 
    geom_bar(stat='identity', width = ttot/24, fill = fill_color) +
    coord_polar(start = -0.261799/2, direction=1) + 
    scale_x_continuous(breaks = hours_breaks, labels = hours_labels, 
                       limits = c(0, ttot)) +
    ylab("") + 
    xlab("Hours") + 
    theme_bw() + 
    ggtitle(title, subtitle = subtitle) +
    theme(aspect.ratio = 1, 
          axis.text=element_text(size=10), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks = element_blank(), 
          plot.title = element_text(hjust = 0.5, size=14),
          plot.subtitle = element_text(hjust = 0.5, size=12))
  
  return(p)
}

#' Create a circular histogram for phase data with 12-hour scale and background colors
#'
#' @param df_original Dataframe containing phase data
#' @param phase_col Column name containing phase values
#' @param gene_col Column name containing gene identifiers (optional)
#' @param ttot Total time period in hours (default: 12)
#' @param title_override Optional title override
#' @param fill_color Color for histogram bars (default: "dodgerblue4")
#' @param bg_color1 Background color for first half (default: "lightyellow")
#' @param bg_color2 Background color for second half (default: "lightblue")
#' @return A ggplot2 object with circular histogram
circular_phase12H_histogram <- function(df_original,
                                        phase_col,
                                        gene_col = NULL,
                                        ttot = 12, 
                                        title_override = NULL,
                                        fill_color = "dodgerblue4",
                                        bg_color1 = "lightyellow",
                                        bg_color2 = "lightblue") {
  # Extract phase data
  x <- df_original[[phase_col]]
  
  # Normalize phase values to the interval [0, ttot]
  x <- x %% ttot
  
  # Create histogram
  br <- seq(0, ttot, length.out = 25)
  h <- hist(x, breaks = br, plot = FALSE)
  
  # Create dataframe for ggplot
  df <- data.frame(x = as.numeric(h$mids), y = h$counts)
  
  # Parameters for hour display
  hours_breaks <- seq(0, ttot, 1)
  hours_labels <- seq(0, 12, 1)
  
  # Set title
  title <- title_override
  
  # Subtitle with number of genes
  n_genes <- length(unique(df_original[[gene_col]]))
  subtitle <- paste(n_genes, "genes")
  
  # Create data for background rectangles
  bg_data1 <- data.frame(
    start = 0,
    end = 6,
    ymin = 0,
    ymax = max(df$y) * 1.1  # Slightly higher than the tallest bar
  )
  
  bg_data2 <- data.frame(
    start = 6,
    end = 12,
    ymin = 0,
    ymax = max(df$y) * 1.1
  )
  
  # Create the plot
  p <- ggplot() +
    # Add colored backgrounds
    geom_rect(data = bg_data1, aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax),
              fill = bg_color1, alpha = 0.4) +
    geom_rect(data = bg_data2, aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax),
              fill = bg_color2, alpha = 0.4) +
    # Add histogram bars
    geom_bar(data = df, aes(x = x, y = y), 
             stat = 'identity', width = ttot/24, fill = fill_color) +
    coord_polar(start = -0.261799/2, direction = 1) + 
    scale_x_continuous(breaks = hours_breaks, labels = hours_labels, 
                       limits = c(0, ttot)) +
    ylab("") + 
    xlab("Hours") + 
    theme_bw() + 
    ggtitle(title, subtitle = subtitle) +
    theme(aspect.ratio = 1, 
          axis.text = element_text(size = 10), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks = element_blank(), 
          plot.title = element_text(hjust = 0.5, size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 12))
  
  return(p)
}

#' Create a circular histogram for phase data with 24-hour scale and background colors
#'
#' @param df_original Dataframe containing phase data
#' @param phase_col Column name containing phase values
#' @param gene_col Column name containing gene identifiers (optional)
#' @param ttot Total time period in hours (default: 24)
#' @param title_override Optional title override
#' @param fill_color Color for histogram bars (default: "dodgerblue4")
#' @param bg_color1 Background color for first half (default: "lightyellow")
#' @param bg_color2 Background color for second half (default: "lightblue")
#' @return A ggplot2 object with circular histogram
circular_phase24H_histogram <- function(df_original,
                                        phase_col,
                                        gene_col = NULL,
                                        ttot = 24, 
                                        title_override = NULL,
                                        fill_color = "dodgerblue4",
                                        bg_color1 = "lightyellow",
                                        bg_color2 = "lightblue") {
  # Extract phase data
  x <- df_original[[phase_col]]
  
  # Normalize phase values to the interval [0, ttot]
  x <- x %% ttot
  
  # Create histogram
  br <- seq(0, ttot, length.out = 25)
  h <- hist(x, breaks = br, plot = FALSE)
  
  # Create dataframe for ggplot
  df <- data.frame(x = as.numeric(h$mids), y = h$counts)
  
  # Parameters for hour display
  hours_breaks <- seq(0, ttot, 1)
  hours_labels <- seq(0, 24, 1)
  
  # Set title
  title <- title_override
  
  # Subtitle with number of genes
  n_genes <- length(unique(df_original[[gene_col]]))
  subtitle <- paste(n_genes, "genes")
  
  # Create data for background rectangles
  bg_data1 <- data.frame(
    start = 0,
    end = 12,
    ymin = 0,
    ymax = max(df$y) * 1.1  # Slightly higher than the tallest bar
  )
  
  bg_data2 <- data.frame(
    start = 12,
    end = 24,
    ymin = 0,
    ymax = max(df$y) * 1.1
  )
  
  # Create the plot
  p <- ggplot() +
    # Add colored backgrounds
    geom_rect(data = bg_data1, aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax),
              fill = bg_color1, alpha = 0.5) +
    geom_rect(data = bg_data2, aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax),
              fill = bg_color2, alpha = 0.5) +
    # Add histogram bars
    geom_bar(data = df, aes(x = x, y = y), 
             stat = 'identity', width = ttot/24, fill = fill_color) +
    coord_polar(start = -0.261799/2, direction = 1) + 
    scale_x_continuous(breaks = hours_breaks, labels = hours_labels, 
                       limits = c(0, ttot)) +
    ylab("") + 
    xlab("Hours") + 
    theme_bw() + 
    ggtitle(title, subtitle = subtitle) +
    theme(aspect.ratio = 1, 
          axis.text = element_text(size = 10), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks = element_blank(), 
          plot.title = element_text(hjust = 0.5, size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 12))
  
  return(p)
}

# ==========================================================================
# DIFFERENTIAL EXPRESSION ANALYSIS FUNCTIONS
# ==========================================================================

#' Perform differential expression analysis with period detection using DESeq2
#'
#' This function tests multiple periods to identify rhythmic genes in a dataset
#' 
#' @param countData Count matrix with genes in rows and samples in columns
#' @param group Vector identifying the group/condition of each sample
#' @param time Vector of time points for each sample
#' @param single Group to analyze (default: first group in 'group')
#' @param sample_name Vector of sample names (default: column names of countData)
#' @param period_range Range of periods to test in hours (default: c(20, 28))
#' @param period_step Step size between periods to test (default: 0.5)
#' @return List containing results for each tested period and best periods for each gene
dryseq_single_multiperiod <- function(countData,
                                      group,
                                      time,
                                      single = group[1],
                                      sample_name = names(countData),
                                      period_range = c(20, 28),
                                      period_step = 0.5) {
  
  # Filter data for the selected condition
  sel = group %in% single
  if(!any(sel)){
    warning("Your single condition is not included in 'group'")
  }
  time = time[sel]
  group = group[sel]
  countData = countData[,sel]
  sample_name = sample_name[sel]
  
  # Remove genes with zero counts
  countData = countData[rowSums(countData)!=0,]
  
  # Create sequence of periods to test
  periods = seq(period_range[1], period_range[2], by=period_step)
  
  # Initialize results list to store outputs for each period
  all_results = list()
  best_results = NULL
  min_p_values = rep(1, nrow(countData))
  best_periods = rep(NA, nrow(countData))
  gene_names = rownames(countData)
  
  # Loop through each period
  for(period in periods) {
    message("Testing period: ", period, " minutes")
    
    # Calculate sine and cosine components for this period
    s1 <- sin(2*pi*time/period)
    c1 <- cos(2*pi*time/period)
    conds = cbind(s1, c1)
    colnames(conds) = c("s1", "c1")
    colData <- data.frame(row.names=colnames(countData), conds)
    
    # Fit rhythms using DESeq2
    dds = DESeq2::DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ s1 + c1)
    dds = DESeq(dds, test="LRT", reduced=~1)
    
    # Extract results and calculate amplitude and phase
    res = as.data.frame(cbind(results(dds), coefficients(dds)))
    res = res[,c('pvalue', 'padj', 'Intercept', 's1', 'c1')]
    phase = period/(2*pi)*atan2(res$s1, res$c1)
    phase = phase%%period
    amp = 2*sqrt(res$s1^2 + res$c1^2)
    res = data.frame(res, phase, amp, period=period)
    
    # Store normalized counts (only need to do this once)
    if(is.null(best_results)) {
      ncounts = counts(dds, normalized = TRUE)
    }
    
    # Store period-specific results
    all_results[[as.character(period)]] = res
    
    # Update best period for each gene if p-value is lower
    for(i in 1:nrow(res)) {
      if(!is.na(res$pvalue[i]) && res$pvalue[i] < min_p_values[i]) {
        min_p_values[i] = res$pvalue[i]
        best_periods[i] = period
      }
    }
  }
  
  # Create best results dataframe using the optimal period for each gene
  best_results = data.frame(
    row.names = gene_names,
    best_period = best_periods,
    min_pvalue = min_p_values
  )
  
  # Add detailed results for each gene's best period
  for(i in 1:nrow(best_results)) {
    gene = gene_names[i]
    best_period = best_results$best_period[i]
    if(!is.na(best_period)) {
      period_key = as.character(best_period)
      best_results$padj[i] = all_results[[period_key]]$padj[i]
      best_results$Intercept[i] = all_results[[period_key]]$Intercept[i]
      best_results$s1[i] = all_results[[period_key]]$s1[i]
      best_results$c1[i] = all_results[[period_key]]$c1[i]
      best_results$phase[i] = all_results[[period_key]]$phase[i]
      best_results$amp[i] = all_results[[period_key]]$amp[i]
    }
  }
  
  # Combine with normalized counts
  global_df = cbind(ncounts, best_results)
  
  # Create output list
  out = list()
  out[["time"]] = time
  out[["period_range"]] = period_range
  out[["period_step"]] = period_step
  out[["results"]] = global_df
  out[["single"]] = single
  out[["all_period_results"]] = all_results
  
  message("finished!")
  return(out)
}

# ==========================================================================
# ENRICHMENT ANALYSIS FUNCTIONS
# ==========================================================================

#' Plot enrichment analysis results for multiple databases
#'
#' @param gene_list List of genes to test for enrichment
#' @param dbs Vector of database names to test
#' @param num_terms Number of terms to display per database (default: 10)
#' @param num_char Maximum number of characters for term names (default: 30)
#' @param max_cols Maximum number of columns in the grid plot (default: 3)
#' @return List containing individual plots and grid plot
plot_enrichment <- function(gene_list, dbs, num_terms = 10, num_char = 30, max_cols = 3) {
  # Perform enrichment analysis
  enrichr_results <- enrichr(gene_list, dbs)
  plots_list <- list()
  
  # For each database in the list
  for (db in dbs) {
    # Check if this database is present in the results
    if (db %in% names(enrichr_results)) {
      # Check if results exist and are not empty
      if (!is.null(enrichr_results[[db]]) && nrow(enrichr_results[[db]]) > 0) {
        # Create plot and add to list
        p <- plotEnrich(enrichr_results[[db]], 
                        showTerms = num_terms,  
                        numChar = num_char,   # Limit number of characters for better readability
                        y = "Count", 
                        orderBy = "Adjusted.P.value",
                        title = db)
        
        plots_list[[db]] <- p
      }
    }
  }
  
  # Determine number of plots
  num_plots <- length(plots_list)
  
  if (num_plots > 0) {
    # Determine optimal layout for the grid
    # Calculate number of columns and rows
    n_cols <- min(max_cols, num_plots)  # Maximum max_cols columns
    n_rows <- ceiling(num_plots / n_cols)
    
    # Create grid with all plots
    grid_plot <- do.call(grid.arrange, c(plots_list, ncol = n_cols))
    
    # Display grid
    print(grid_plot)
    
    # Return list of plots and grid for later use
    return(list(plots = plots_list, grid = grid_plot))
  } else {
    message("No enrichment results found for the specified databases.")
    return(NULL)
  }
}


#' Calcule les moyennes des mesures d'expression génique pour chaque point temporel ZT
#'
#' Cette fonction prend une matrice de données d'expression génique avec des réplicats 
#' techniques (a et b) pour chaque point temporel ZT (Zeitgeber Time) et calcule 
#' la moyenne entre ces réplicats.
#'
#' @param data Un data frame contenant les données d'expression. Les colonnes doivent 
#'             être nommées selon le format "ZT_Xa" et "ZT_Xb", où X est le temps ZT.
#' @return Un data frame contenant les moyennes calculées pour chaque temps ZT.
#'
#' @details Les colonnes du data frame d'entrée doivent suivre la convention de 
#'          nommage "ZT_Xa" et "ZT_Xb" où X est le temps Zeitgeber et "a"/"b" 
#'          indiquent les réplicats techniques.
#'
calculer_moyennes_ZT <- function(data) {
  # Extraire les noms de colonnes
  colonnes <- colnames(data)
  
  # Définir un pattern regex pour capturer le temps ZT sans la lettre a/b
  # Format attendu: "ZT_X[a|b]" où X peut être un nombre entier ou décimal
  zeitgeber_pattern <- "ZT_(\\d+(\\.\\d+)?)[ab]"
  
  # Extraire toutes les correspondances au pattern
  zeitgeber_matches <- regmatches(colonnes, regexpr(zeitgeber_pattern, colonnes))
  
  # Supprimer les suffixes a/b et obtenir les temps ZT uniques
  zeitgeber_times <- unique(gsub("[ab]$", "", zeitgeber_matches))
  
  # Créer un nouveau tableau pour stocker les moyennes calculées
  data_moyennes <- data.frame(matrix(nrow = nrow(data), ncol = length(zeitgeber_times)))
  rownames(data_moyennes) <- rownames(data)  # Conserver les identifiants de gènes
  colnames(data_moyennes) <- zeitgeber_times  # Utiliser les temps ZT comme noms de colonnes
  
  # Calculer les moyennes pour chaque point temporel ZT
  for (zt in gsub("ZT_", "", zeitgeber_times)) {
    # Construire les noms de colonnes pour les réplicats a et b
    col_a <- paste0("ZT_", zt, "a")
    col_b <- paste0("ZT_", zt, "b")
    
    # Calculer la moyenne entre les réplicats (ignorer NA si présent)
    data_moyennes[, paste0("ZT_", zt)] <- rowMeans(data[, c(col_a, col_b)], na.rm = TRUE)
  }
  
  return(data_moyennes)
}


#' Visualise l'expression d'un gène sur un cycle circadien et affiche les résultats RAIN
#'
#' Cette fonction crée un graphique montrant l'expression d'un gène spécifique au fil du temps,
#' avec la possibilité d'inclure les résultats des analyses de rythmicité RAIN.
#'
#' @param gene_name Le nom du gène à visualiser (doit correspondre à un nom dans l'expression_matrix)
#' @param expression_matrix Une matrice d'expression avec les gènes en lignes et les temps ZT en colonnes
#' @param rain_df_list Une liste de data frames contenant les résultats d'analyses RAIN (optionnel)
#' @param df_names Noms des jeux de données dans rain_df_list (optionnel)
#' @param color Couleur à utiliser pour les lignes et points du graphique
#' @return Une liste contenant le graphique (plot) et les résultats RAIN (rain_results)
#'
#' @details Le graphique généré inclut:
#'   - Une ligne d'expression moyenne sur le cycle circadien
#'   - Des barres d'erreur représentant l'erreur standard
#'   - Un fond coloré pour différencier jour (jaune) et nuit (bleu)
#'   - Des annotations avec les résultats de l'analyse RAIN (si fournie)
#'
plot_gene_expression <- function(gene_name, expression_matrix, 
                                 rain_df_list = list(), 
                                 df_names = NULL,
                                 color = "black") {
  # Extraire la ligne correspondant au gène recherché
  gene_data <- expression_matrix[rownames(expression_matrix) == gene_name, , drop = FALSE]
  
  # Vérifier que le gène existe dans nos données
  if (nrow(gene_data) == 0) {
    stop(paste("Le gène", gene_name, "n'a pas été trouvé dans la matrice d'expression"))
  }
  
  # Transformer les données au format long pour faciliter la visualisation avec ggplot2
  gene_long <- gene_data %>%
    rownames_to_column("Gene") %>%  # Convertir l'identifiant du gène en colonne
    pivot_longer(cols = -Gene, names_to = "ZT", values_to = "expression") %>%  # Format long
    mutate(
      Zeitgeber = as.numeric(str_extract(ZT, "\\d+")),  # Extraire le temps ZT
      Replicate = str_extract(ZT, "[ab]")              # Extraire l'identifiant du réplicat
    )
  
  # Calculer les moyennes et erreurs standard pour chaque point temporel
  resume_donnees <- gene_long %>%
    group_by(Gene, Zeitgeber) %>%
    summarise(
      expression_moyenne = mean(expression, na.rm = TRUE),
      erreur_std = sd(expression, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"  # Éviter le message d'avertissement sur les groupements
    )
  
  # Créer les rectangles de fond pour différencier jour/nuit
  # Jaune clair pour le jour (ZT0-ZT12) et bleu clair pour la nuit (ZT12-ZT24)
  background <- data.frame(
    xmin = c(0, 12),
    xmax = c(12, 22.5),  # Limité à 22.5 au lieu de 24 pour accommoder les données disponibles
    ymin = -Inf,
    ymax = Inf,
    fill = c("lightyellow", "lightblue")
  )
  
  # Préparer un data frame pour stocker les résultats des analyses RAIN
  rain_results <- data.frame(Dataset = character(),
                             pValue = numeric(),
                             Period = numeric(),
                             Significant = logical(),
                             stringsAsFactors = FALSE)
  
  # Variable pour stocker le texte des annotations
  annotation_text <- ""
  
  # Traiter les résultats RAIN si fournis
  if (length(rain_df_list) > 0) {
    # Si les noms des datasets ne sont pas fournis, créer des noms par défaut
    if (is.null(df_names)) {
      df_names <- paste0("Dataset ", 1:length(rain_df_list))
    }
    
    # Parcourir chaque jeu de données avec analyses RAIN
    for (i in 1:length(rain_df_list)) {
      df <- rain_df_list[[i]]
      df_name <- df_names[i]
      
      # Vérifier si le gène est présent dans ce jeu de données
      if (gene_name %in% rownames(df)) {
        # Extraire p-value et période
        p_val <- df[gene_name, "pVal"]
        period <- df[gene_name, "period"]
        significant <- p_val < 0.05  # Considérer significatif si p < 0.05
        
        # Ajouter au data frame de résultats
        rain_results <- rbind(rain_results, 
                              data.frame(Dataset = df_name,
                                         pValue = p_val,
                                         Period = period,
                                         Significant = significant))
        
        # Ajouter au texte d'annotation
        annotation_text <- paste0(annotation_text, 
                                  df_name, ": padj = ", sprintf("%.4f", p_val), 
                                  ", Période = ", sprintf("%.1f", period), "\n")
      } else {
        # Si le gène n'est pas trouvé dans ce jeu de données
        rain_results <- rbind(rain_results, 
                              data.frame(Dataset = df_name,
                                         pValue = NA,
                                         Period = NA,
                                         Significant = FALSE))
        
        annotation_text <- paste0(annotation_text, df_name, ": non trouvé\n")
      }
    }
    
    # Supprimer le dernier retour à la ligne
    annotation_text <- substr(annotation_text, 1, nchar(annotation_text) - 1)
  }
  
  # Essayer de déterminer le nom du tissu à partir du nom de la matrice d'expression
  tissue_name <- tryCatch({
    sub("matrix_log_(.+)_rep", "\\1", deparse(substitute(expression_matrix)))
  }, error = function(e) {
    "tissu"  # Nom par défaut si l'extraction échoue
  })
  
  # Créer le graphique
  plot <- ggplot(resume_donnees, aes(x = Zeitgeber, y = expression_moyenne, group = Gene)) +
    # Ajouter les rectangles de fond jour/nuit
    geom_rect(data = background, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
              inherit.aes = FALSE, alpha = 0.3) +
    # Ajouter la ligne d'expression
    geom_line(color = color, linewidth = 1.5) +
    # Ajouter les points pour chaque mesure
    geom_point(size = 3, color = color) +
    # Ajouter les barres d'erreur (erreur standard)
    geom_errorbar(aes(ymin = expression_moyenne - erreur_std, 
                      ymax = expression_moyenne + erreur_std), 
                  width = 0.4, color = color) +
    # Ajouter les labels et le titre
    labs(x = "ZT (Zeitgeber Time)",
         y = "mean(log2 normalized counts)",
         title = paste0(gene_name)) +
    # Utiliser un thème épuré
    theme_classic(base_size = 14) +
    # Utiliser les couleurs définies dans le data frame background
    scale_fill_identity()
  
  # Ajouter l'annotation avec les résultats RAIN si disponibles
  if (nchar(annotation_text) > 0) {
    # Calculer la plage des valeurs d'expression pour positionner l'annotation
    y_range <- diff(range(resume_donnees$expression_moyenne, na.rm = TRUE))
    y_pos <- min(resume_donnees$expression_moyenne, na.rm = TRUE) + y_range * 0.05
    
    # Ajouter le texte d'annotation
    plot <- plot +
      annotate("text", x = 30, 
               y = max(resume_donnees$expression_moyenne) + 1, 
               label = annotation_text, 
               hjust = 1, vjust = 1, 
               size = 4,
               color = "black")
  }
  
  # Afficher le graphique
  print(plot)
  
  # Afficher également les résultats dans la console
  if (nrow(rain_results) > 0) {
    cat("\nRésultats de l'analyse rythmique pour", gene_name, ":\n")
    cat("--------------------------------------------------\n")
    
    for (i in 1:nrow(rain_results)) {
      if (!is.na(rain_results$pValue[i])) {
        cat(sprintf("%s: p-value = %.4f, période = %.1f\n", 
                    rain_results$Dataset[i], 
                    rain_results$pValue[i], 
                    rain_results$Period[i]))
      } else {
        cat(sprintf("%s: %s non trouvé\n", rain_results$Dataset[i], gene_name))
      }
    }
  }
  
  # Retourner une liste avec le graphique et les résultats
  return(list(plot = plot, rain_results = rain_results))
}

