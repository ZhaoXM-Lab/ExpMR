rm(list=ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)
library(stats)

f1 = fread('results/signif_genes.CellToIDP.fdr.txt', data.table = F)
f2 = fread('results/signif_genes.CellToDisorder.fdr.txt', data.table = F)

IDP_all = unique(f1$pheno)
disorder_all = unique(f2$pheno)

IDP = IDP_all[1]; disorder = disorder_all[2]
result0 = data.frame(stringsAsFactors = F)
for (IDP in IDP_all){
  for (disorder in disorder_all){
    sub_f1 = subset(f1, f1$pheno == IDP)
    sub_f2 = subset(f2, f2$pheno == disorder)
    
    gene_list1 = paste(sub_f1$tissue, sub_f1$gene_name, sep='-')
    gene_list2 = paste(sub_f2$tissue, sub_f2$gene_name, sep='-')
    inter_gene = intersect(gene_list1, gene_list2)
    if (length(inter_gene) == 0) p_hyper = 1
    
    if (length(inter_gene) > 0){
      q = length(inter_gene)
      m = length(gene_list1)
      n = 20356 * 8 - m
      k = length(gene_list2)
      p_hyper = phyper(q-1, m, n, k, lower.tail=F)
    }
    
    result0.tmp = data.frame(IDP = IDP, Disorder = disorder, Causal_genes1 = paste(gene_list1,collapse=','),
                             Causal_genes2 = paste(gene_list2, collapse=','),
                             Overlap_genes = paste(inter_gene,collapse=','), p_hyper = p_hyper, stringsAsFactors = F)
    result0 = rbind(result0, result0.tmp)
    cat(IDP, 'and', disorder, 'is OK.', '\n') 
  }
}
result0$fdr = p.adjust(result0$p_hyper, method = 'fdr', n = nrow(result0))
result = subset(result0, result0$p_hyper < 1)


##heatmap
trait_set1 = unique(result$IDP); trait_set2 = unique(result$Disorder)
matrix1 = matrix(1, length(trait_set1), length(trait_set2))
rownames(matrix1) = trait_set1  #IDP
colnames(matrix1) = trait_set2  #Brain disorder
for (i in 1:nrow(result)){
  trait1 = result[i, 'IDP']
  trait2 = result[i, 'Disorder']
  matrix1[trait1, trait2] = result[i, 'p_hyper']
}
matrix1 = -log10(matrix1)

matrix2 = matrix('', length(trait_set1), length(trait_set2))
rownames(matrix2) = trait_set1
colnames(matrix2) = trait_set2
for (i in 1:nrow(result)){
  trait1 = result[i, 'IDP']
  trait2 = result[i, 'Disorder']
  if (result[i,'fdr'] < 0.05) matrix2[trait1, trait2] = '*'
  if (result[i,'fdr'] < 0.005) matrix2[trait1, trait2] = '**'
}


trait_mapping = read.csv('gwas_meta_all.csv', header = T)
trait_mapping = trait_mapping[, c('TraitCategory1','TraitCategory2','TraitAbbr')]
annotation_col = data.frame(trait_mapping$TraitCategory2, trait_mapping$TraitCategory1)
rownames(annotation_col) = trait_mapping$TraitAbbr
names(annotation_col) = c('TraitCategory2', 'TraitCategory1')


ann_colors = list(TraitCategory1 = c(`Behavioral-cognitive phenotype` = "#E7298A", 
                                     `Psychiatric disorder` = "#66A61E",  `Neurological disorder` = "#D95F02"))
a = pheatmap::pheatmap(matrix1, 
                       display_numbers = matrix2, fontsize_number = 12, number_color = "black",
                       fontsize_color = 'black', fontsize_col = 10, fontsize_row = 11,
                       color = colorRampPalette(colors = c("white","#FA8072","red"))(10),
                       annotation_col = annotation_col, annotation_colors = ann_colors, 
                       fontfamily= "serif")
tree_col0 = trait_set2[ a[["tree_col"]][["order"]] ]
annotation_col = subset(annotation_col, row.names(annotation_col)%in% trait_set2)
PD = row.names(annotation_col)[annotation_col$TraitCategory1 == 'Psychiatric disorder']
ND = row.names(annotation_col)[annotation_col$TraitCategory1 == 'Neurological disorder']
BCP = row.names(annotation_col)[annotation_col$TraitCategory1 == 'Behavioral-cognitive phenotype']

Cognitive = row.names(annotation_col)[annotation_col$TraitCategory2 == 'Cognitive']
Personality = row.names(annotation_col)[annotation_col$TraitCategory2 == 'Personality']
Risky = row.names(annotation_col)[annotation_col$TraitCategory2 == 'Risky behavior']
Chronic = row.names(annotation_col)[annotation_col$TraitCategory2 == 'Chronic disorder']
Episodic = row.names(annotation_col)[annotation_col$TraitCategory2 == 'Episodic disorder']
Neurodegenerative = row.names(annotation_col)[annotation_col$TraitCategory2 == 'Neurodegenerative disease']
Compulsive = row.names(annotation_col)[annotation_col$TraitCategory2 == 'Compulsive disorder']
Internalizing = row.names(annotation_col)[annotation_col$TraitCategory2 == 'Internalizing disorder']
Neurodevelopmental = row.names(annotation_col)[annotation_col$TraitCategory2 == 'Neurodevelopmental disorder']
Psychotic = row.names(annotation_col)[annotation_col$TraitCategory2 == 'Psychotic disorder']


tree_col = c(tree_col0[tree_col0 %in% PD & tree_col0 %in% Psychotic], 
             tree_col0[tree_col0 %in% PD & tree_col0 %in% Neurodevelopmental],
             tree_col0[tree_col0 %in% PD & tree_col0 %in% Internalizing],
             tree_col0[tree_col0 %in% PD & tree_col0 %in% Compulsive],
             tree_col0[tree_col0 %in% ND & tree_col0 %in% Neurodegenerative],
             tree_col0[tree_col0 %in% ND & tree_col0 %in% Episodic],
             tree_col0[tree_col0 %in% ND & tree_col0 %in% Chronic],
             tree_col0[tree_col0 %in% BCP & tree_col0 %in% Risky],
             tree_col0[tree_col0 %in% BCP & tree_col0 %in% Personality],
             tree_col0[tree_col0 %in% BCP & tree_col0 %in% Cognitive])
tree_col[15:19] = c("SC","DPW","GRT","ASP","NEU")


matrix1 = matrix1[, match(tree_col, colnames(matrix1))]
matrix2 = matrix2[, match(tree_col, colnames(matrix2))]

pheatmap::pheatmap(matrix1, cluster_col = F,
                   display_numbers = matrix2, fontsize_number = 12, number_color = "black",
                   fontsize_color = 'black', fontsize_col = 10, fontsize_row = 11,
                   color = colorRampPalette(colors = c("white","#FA8072","red"))(10),
                   annotation_col = annotation_col, annotation_colors = ann_colors,
                   #gaps_col = c(9, 14), border_col = 'grey',
                   fontfamily= "serif")


trait_mapping = data.frame(trait_set1)
trait_mapping$`Trait Type` = 'Brain volume'
trait_mapping$`Trait Type`[grep('FA', trait_mapping$trait_set1)] = 'White matter microstructure'
rownames(trait_mapping) = trait_mapping$trait_set1
annotation_row = data.frame(trait_mapping$`Trait Type`)
rownames(annotation_row) = trait_mapping$trait_set1
names(annotation_row) = 'TraitCategory3'

BV = row.names(annotation_row)[which(annotation_row$`TraitCategory3` == 'Brain volume')] 
WM = row.names(annotation_row)[which(annotation_row$`TraitCategory3` == 'White matter microstructure')]
tree_row0 = trait_set1[ a[["tree_row"]][["order"]] ]

tree_row = c(tree_row0[which(tree_row0 %in% BV)], tree_row0[sort(which(tree_row0 %in% WM),decreasing=T)])

matrix1 = matrix1[match(tree_row, rownames(matrix1)), ]
matrix2 = matrix2[match(tree_row, rownames(matrix2)), ]

library(cols4all)
c4a_palettes("div")
mycol <- c4a("rainbow_wh_rd", 25)
c4a_plot(mycol)
set.seed(12)
mycol <- sample(mycol, 15)
ann_colors = list(`Trait_Type_col` = c(`Behavioral-cognitive phenotype` = "#E7298A", 
                                       `Psychiatric disorder` = "#66A61E",  `Neurological disorder` = "#D95F02"),
                  `Trait_Type_row` = c(`Brain volume` = "#FEB90E",  `White matter microstructure` ="#C58CA9"))
ann_colors = list(`TraitCategory1` = c(`Behavioral-cognitive phenotype` =mycol[1], 
                                       `Psychiatric disorder`=mycol[2],
                                       `Neurological disorder`=mycol[3]),
                  `TraitCategory2` = c(`Cognitive`=mycol[4],
                                       `Personality`=mycol[5],
                                       `Risky behavior`=mycol[6],
                                       `Chronic disorder`=mycol[7],
                                       `Episodic disorder`=mycol[8],
                                       `Neurodegenerative disease`=mycol[9],
                                       `Compulsive disorder`=mycol[10],
                                       `Internalizing disorder`=mycol[11],
                                       `Neurodevelopmental disorder`=mycol[12],
                                       `Psychotic disorder`=mycol[13]),    
                  `TraitCategory3` = c(`Brain volume`="#FEB90E",
                                       `White matter microstructure`="#C58CA9"))
pheatmap::pheatmap(matrix1, 
                   cluster_col = F, cluster_rows = F,
                   display_numbers = matrix2, fontsize_number = 12, number_color = "black",
                   fontsize_color = 'black', fontsize_col = 10, fontsize_row = 11,
                   color = colorRampPalette(colors = c("white","red","#a90101"))(20),
                   annotation_col = annotation_col, 
                   annotation_row = annotation_row, 
                   annotation_colors = ann_colors,
                   annotation_names_row = F,
                   border_col = 'black',
                   angle_col = 45,
                   gaps_col = c(9, 14),
                   fontfamily= "serif")

matrix1_t = t(matrix1)
matrix2_t = t(matrix2)

pheatmap::pheatmap(matrix1_t, 
                   cluster_col = F, cluster_rows = F,
                   display_numbers = matrix2_t, fontsize_number = 12, number_color = "black",
                   fontsize_color = 'black', fontsize_col = 10, fontsize_row = 11,
                   color = colorRampPalette(colors = c("white","red","#a90101"))(20),
                   annotation_col = annotation_row, 
                   annotation_row = annotation_col, 
                   annotation_colors = ann_colors,
                   annotation_names_row = F,
                   border_col = 'black',
                   angle_col = 90,
                   gaps_row = c(9, 14),
                   fontfamily= "serif")

