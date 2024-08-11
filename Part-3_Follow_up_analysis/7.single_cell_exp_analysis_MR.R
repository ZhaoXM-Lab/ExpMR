library(data.table)
library(stringi)
library(dplyr)
library(ggpubr)
library(export)
library(colorspace)
library(stringr)
library(corrplot)
library(ggplot2)
library(reshape2)
library(pheatmap)

# Load data
load_data <- function() {
  IDP_single_cell <- fread('IDP_per_cell_type_exp_end.csv', data.table = FALSE)
  disorder_single_cell <- fread('Dis_per_cell_type_exp_end.csv', data.table = FALSE)
  disorder_single_cell <- disorder_single_cell[which(disorder_single_cell$cell_type != 'Endothelial.cells' & disorder_single_cell$cell_type != 'Pericytes'),]
  disorder_single_cell$period_group <- disorder_single_cell$period
  disorder_single_cell[which(disorder_single_cell$period <= 7), 6] <- 1
  disorder_single_cell[which(disorder_single_cell$period > 7), 6] <- 2
  
  cell_type_anno <- data.frame(
    cell_name = c('Microglia', 'Astrocytes', 'Oligodendrocytes', 'OPCs...COPs', 'Excitatory.neurons', 'Inhibitory.neurons'),
    cell_class = c('Glia', 'Glia', 'Glia', 'Glia', 'Neuro', 'Neuro')
  )
  rownames(cell_type_anno) <- cell_type_anno$cell_name
  
  color_cell_type <- c('Astrocytes', 'Excitatory.neurons', 'Inhibitory.neurons', 'Microglia', 'Oligodendrocytes', 'OPCs...COPs')
  color_code <- c('#F0027F', '#BEAED4', '#FDC086', '#FFFF99', '#386CB0', '#4DC97F')
  cell_color_code <- data.frame(cell_type = color_cell_type, color_code = color_code)
  rownames(cell_color_code) <- cell_color_code$cell_type
  
  gwas_meta <- fread('/gwas_meta.csv', data.table = FALSE)
  image_pheno <- fread('./image_pheno_name.txt', sep = '\t', fill = TRUE, data.table = FALSE)
  image_pheno <- image_pheno[order(image_pheno$group),]
  rownames(image_pheno) <- image_pheno$x
  
  disorder_single_cell$pheno <- factor(disorder_single_cell$pheno, levels = gwas_meta$TraitAbbr)
  disorder_single_cell$period_group <- as.character(disorder_single_cell$period_group)
  disorder_single_cell$cell_class <- cell_type_anno[disorder_single_cell$cell_type, 2]
  rownames(gwas_meta) <- gwas_meta$TraitAbbr
  disorder_single_cell$pheno_class <- gwas_meta[disorder_single_cell$pheno, 1]
  
  list(
    IDP_single_cell = IDP_single_cell,
    disorder_single_cell = disorder_single_cell,
    cell_type_anno = cell_type_anno,
    cell_color_code = cell_color_code,
    gwas_meta = gwas_meta,
    image_pheno = image_pheno
  )
}

# Plot boxplot for cell classes
plot_cell_classes <- function(disorder_single_cell) {
  glia_cell_d <- disorder_single_cell[disorder_single_cell$cell_class == "Glia",]
  neuro_cell_d <- disorder_single_cell[disorder_single_cell$cell_class == "Neuro",]
  
  # Significance test
  wilcox.test(glia_cell_d$exp, neuro_cell_d$exp)
  
  # Boxplot
  ggboxplot(disorder_single_cell, x = 'pheno', y = 'exp', fill = 'cell_class', sort.val = "median", sort.by = "desc", x.text.angle = 30, outlier.shape = NA) +
    stat_compare_means(aes(group = cell_class), label = "p.signif") +
    scale_fill_manual(values = c('#ACA2C7', '#EFB000'))
}

# Plot boxplot for cell types
plot_cell_types <- function(disorder_single_cell, cell_color_code) {
  cell_order <- c('Astrocytes', 'Oligodendrocytes', 'Microglia', 'OPCs...COPs', 'Inhibitory.neurons', 'Excitatory.neurons')
  disorder_single_cell$cell_type <- factor(disorder_single_cell$cell_type, levels = cell_order)
  cell_order_color <- cell_color_code[cell_order, 2]
  
  ggboxplot(disorder_single_cell, x = 'cell_type', y = 'exp', sort.by = "desc", fill = 'cell_type', ylab = 'Expression', outlier.shape = NA, x.text.angle = 30) +
    scale_fill_manual(values = cell_order_color)
}

# Subclass analysis for different phenotype classes
subclass_analysis <- function(disorder_single_cell, gwas_meta) {
  psy_dis <- disorder_single_cell[which(disorder_single_cell$pheno_class == "Psychiatric disorder"),]
  neuro_dis <- disorder_single_cell[which(disorder_single_cell$pheno_class == "Neurological disorder"),]
  behaviour_dis <- disorder_single_cell[which(disorder_single_cell$pheno_class == "Behavioral-cognitive phenotype"),]
  
  p1 <- ggboxplot(psy_dis, x = 'pheno', y = 'exp', fill = 'period_group', xlab = 'Phenotypes', ylab = 'Expression', facet.by = 'cell_class') +
    stat_compare_means(aes(group = period_group), label = "p.signif") +
    scale_fill_manual(values = c('#D9D9D9', '#A6A6A6'))
  
  p2 <- ggboxplot(neuro_dis, x = 'pheno', y = 'exp', fill = 'period_group', xlab = 'Phenotypes', ylab = 'Expression', facet.by = 'cell_class', outlier.shape = NA) +
    stat_compare_means(aes(group = period_group), label = "p.signif") +
    scale_fill_manual(values = c('#D9D9D9', '#A6A6A6'))
  
  p3 <- ggboxplot(behaviour_dis, x = 'pheno', y = 'exp', fill = 'period_group', xlab = 'Phenotypes', ylab = 'Expression', facet.by = 'cell_class', outlier.shape = NA) +
    stat_compare_means(aes(group = period_group), label = "p.signif") +
    scale_fill_manual(values = c('#D9D9D9', '#A6A6A6'))
  
  ggarrange(p1, p2, p3, ncol = 1, nrow = 3)
}

# Plot trajectories for different cell types
plot_trajectories <- function(df_exp) {
  traject <- data.frame()
  for (i in unique(df_exp$pheno)) {
    temo_df <- df_exp[which(df_exp$pheno == i),]
    exp_perd <- tapply(temo_df$exp, temo_df$period, median)
    tem_re <- data.frame(pheno = rep(i, length(exp_perd)), period_name = as.numeric(names(exp_perd)), exp_d = exp_perd)
    traject <- rbind(traject, tem_re)
  }
  traject$period_name <- paste('P', traject$period_name, sep = "")
  traject$period_name <- factor(traject$period_name, levels = unique(traject$period_name))
  
  ggplot(data = traject, aes(x = period_name, y = exp_d, group = pheno, color = factor(pheno), fill = factor(pheno))) +
    geom_smooth(method = 'loess') +
    xlab('Periods') + ylab("Normalized expression") + theme_bw() +
    theme(panel.grid = element_blank())
}

# Plot trajectories for each cell type
plot_each_expression <- function(data_exp, cell_order) {
  plots <- lapply(cell_order, function(cell_type) {
    plot_trajectories(data_exp[which(data_exp$cell_type == cell_type),])
  })
  
  p <- ggarrange(plotlist = plots, ncol = 2, nrow = 3, labels = c("A", "B", "C", "D", "E", "F"))
  print(p)
  graph2ppt(file = 'disorder_prepare_new.pptx', height = 15, width = 8, append = TRUE)
}

# Preprocess data for correlation analysis
del_exp_data <- function(exp_df, all_cell_types) {
  save_df <- data.frame()
  for (i in unique(exp_df$pheno)) {
    rep_e <- data.frame(names = all_cell_types, exp = rep(0, length(all_cell_types)), pheno = rep(i, length(all_cell_types)))
    rownames(rep_e) <- rep_e$names
    tem_exp_df <- exp_df[which(exp_df$pheno == i),]
    emp_ex <- tapply(tem_exp_df$exp, tem_exp_df$cell_type, mean)
    rep_e[names(emp_ex), "exp"] <- emp_ex
    save_df <- rbind(save_df, rep_e)
  }
  return(save_df)
}

# Calculate correlation matrix between brain disorders and brain volume/white matter
calculate_correlation_matrix <- function(disorder_cell_exp, feature_all, all_cell_types) {
  disorder_ph_u <- unique(disorder_cell_exp$pheno)
  brain_ph_u <- unique(feature_all$pheno)
  
  cor_matrix <- matrix(0, length(disorder_ph_u), length(brain_ph_u))
  cor_p_matrix <- matrix(0, length(disorder_ph_u), length(brain_ph_u))
  
  for (i in seq_along(disorder_ph_u)) {
    temp_dis <- disorder_cell_exp[which(disorder_cell_exp$pheno == disorder_ph_u[i]),]
    temp_dis[is.na(temp_dis)] <- 0
    
    for (j in seq_along(brain_ph_u)) {
      temp_fea <- feature_all[which(feature_all$pheno == brain_ph_u[j]),]
      
      if (nrow(temp_fea[which(temp_fea$exp == 0),]) <= 6 & nrow(temp_dis[which(temp_dis$exp == 0),]) <= 6) {
        cor_matrix[i, j] <- cor(temp_dis$exp, temp_fea$exp, method = 'spearman')
        P_tem <- cor.test(temp_dis$exp, temp_fea$exp, method = 'spearman', exact = FALSE, monte.carlo = TRUE)
        cor_p_matrix[i, j] <- P_tem$p.value
      } else {
        cor_p_matrix[i, j] <- 1
        cor_matrix[i, j] <- 0
      }
    }
  }
  
  rownames(cor_matrix) <- disorder_ph_u
  colnames(cor_matrix) <- brain_ph_u
  rownames(cor_p_matrix) <- disorder_ph_u
  colnames(cor_p_matrix) <- brain_ph_u
  
  list(cor_matrix = cor_matrix, cor_p_matrix = cor_p_matrix)
}

# Adjust p-values
adjust_p_values <- function(cor_p) {
  cor_p_vector <- as.vector(cor_p)
  corrected_p_values <- p.adjust(cor_p_vector, method = "BH")
  corrected_p_values_matrix <- matrix(corrected_p_values, nrow = nrow(cor_p))
  rownames(corrected_p_values_matrix) <- rownames(cor_p)
  colnames(corrected_p_values_matrix) <- colnames(cor_p)
  return(corrected_p_values_matrix)
}

# Plot adjusted correlation matrix
plot_adjusted_correlation <- function(cor_matrix, cor_p_matrix, annotation_colors, bh_adj = TRUE) {
  if (bh_adj) {
    cor_p_matrix <- adjust_p_values(cor_p_matrix)
  }
  
  rows_to_keep_cor <- apply(cor_matrix, 1, function(row) any(row != 0))
  cols_to_keep_cor <- apply(cor_matrix, 2, function(col) any(col != 0))
  rows_to_keep_cor <- names(rows_to_keep_cor[which(rows_to_keep_cor == TRUE)])
  cols_to_keep_cor <- names(cols_to_keep_cor[which(cols_to_keep_cor == TRUE)])
  cor_matrix <- cor_matrix[rows_to_keep_cor, cols_to_keep_cor]
  cor_p_matrix <- cor_p_matrix[rows_to_keep_cor, cols_to_keep_cor]
  
  col_va <- colorRampPalette(c("#2B73B3", "white", "#B2182B"))(10)
  corrplot(corr = as.matrix(cor_matrix), p.mat = as.matrix(cor_p_matrix), outline = TRUE, method = "color", insig = "label_sig",
           sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.2, cl.cex = 0.8, tl.cex = 1, col = col_va, tl.col = "black")
}

# Main function to run the entire analysis
run_analysis <- function() {
  data <- load_data()
  
  plot_cell_classes(data$disorder_single_cell)
  plot_cell_types(data$disorder_single_cell, data$cell_color_code)
  
  subclass_analysis(data$disorder_single_cell, data$gwas_meta)
  
  cell_order <- c('Astrocytes', 'Oligodendrocytes', 'Microglia', 'OPCs...COPs', 'Inhibitory.neurons', 'Excitatory.neurons')
  plot_each_expression(data$disorder_single_cell, cell_order)
  
  # Calculate and plot correlation matrices
  all_cell_types <- unique(data$disorder_single_cell$cell_type)
  disorder_cell_exp <- del_exp_data(data$disorder_single_cell, all_cell_types)
  feature_all <- rbind(del_exp_data(data$IDP_single_cell, all_cell_types), del_exp_data(data$IDP_single_cell, all_cell_types))
  
  cor_results <- calculate_correlation_matrix(disorder_cell_exp, feature_all, all_cell_types)
  plot_adjusted_correlation(cor_results$cor_matrix, cor_results$cor_p_matrix, annotation_colors = list(
    Annotation = c("Neurological disorder" = "#E44D00", "Psychiatric disorder" = "#FEB500", "Behavioral-cognitive phenotype" = "#C58CA9")
  ))
}

# Run the analysis
run_analysis()
