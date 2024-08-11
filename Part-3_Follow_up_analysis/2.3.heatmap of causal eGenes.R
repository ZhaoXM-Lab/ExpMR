#Heatmap for IDPs###############################################################################################
rm(list=ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)
library(tidyr)


calcu_Ngene = function(data){
  Ngene = data.frame(stringsAsFactors = F)
  
  pheno_all = unique(data$pheno)
  tissue_all = unique(data$tissue)
  for (pheno_idx in pheno_all){
    for (tissue_idx in tissue_all){
      subset_file = subset(data, data$tissue == tissue_idx & data$pheno == pheno_idx)
      N = nrow(subset_file)
      Ngene.tmp = data.frame(pheno = pheno_idx, tissue = tissue_idx, N = N)
      Ngene = rbind(Ngene.tmp, Ngene)
    }
  }
  return(Ngene)
}
f1 = fread('results/signif_genes.CellToIDP.fdr.txt', data.table = F)
f1 = f1[, c("pheno", "tissue","gene","gene_name")]



Ngene1 = calcu_Ngene(f1)
Ngene = Ngene1
Ngene = Ngene[-which((Ngene$N) == 0), ]
data = Ngene[, c("pheno","tissue","N")]


matrix0 = matrix(0,  length(unique(Ngene$pheno)), length(unique(Ngene$tissue))) #行-cell, 列-gwas文件
rownames(matrix0) = unique(Ngene$pheno)
colnames(matrix0) = unique(Ngene$tissue)
for (i in 1:nrow(Ngene)){
  pheno = Ngene$pheno[i]
  tissue = Ngene$tissue[i]
  matrix0[which(rownames(matrix0)==pheno), which(colnames(matrix0) == tissue)] = Ngene[i,'N']
}

row_sums = rowSums(matrix0) 
matrix = matrix0[row_sums >= 3, ]
matrix = t(matrix)



trait_mapping = data.frame(pheno = colnames(matrix))
trait_mapping$`TraitCategory` = 'Brain regional volume'
trait_mapping$`TraitCategory`[grep('FA', trait_mapping$pheno)] = 'White matter microstructure'
rownames(trait_mapping) = trait_mapping$pheno
annotation = data.frame(trait_mapping$`TraitCategory`)
rownames(annotation) = trait_mapping$pheno
names(annotation) = 'TraitCategory'


library(cols4all)
c4a_palettes("div")
mycol <- c4a("rainbow_wh_rd", 25)
c4a_plot(mycol)
set.seed(12)
mycol <- sample(mycol, 15)

ann_colors = list(TraitCategory = c(`Brain regional volume` = '#FEB90E', 
                                    `White matter microstructure` = '#CC5D5B'))

a = pheatmap::pheatmap(matrix, 
                   fontsize_number = 12, number_color = "black",
                   fontsize_color = 'black', fontsize_col = 10,fontsize_row = 11,
                   color = colorRampPalette(colors = c("white","#FA8072","red"))(20),
                   annotation_col = annotation, 
                   annotation_colors = ann_colors,
                   angle_col = 45,
                   fontfamily= "serif")

BV = row.names(annotation)[which(annotation$`TraitCategory` == 'Brain regional volume')]
WM = row.names(annotation)[which(annotation$`TraitCategory` == 'White matter microstructure')]
tree_col0 = colnames(matrix)[ a[["tree_col"]][["order"]] ]
tree_col = c(tree_col0[which(tree_col0 %in% BV)], tree_col0[sort(which(tree_col0 %in% WM),decreasing=T)])
matrix = matrix[, tree_col]

tree_row = rownames(matrix)[ a[["tree_row"]][["order"]] ]
matrix = matrix[tree_row,]

matrix2 = matrix 
matrix2[which(matrix2 == 0)] = ''
pheatmap::pheatmap(matrix, 
                   cluster_col = F, cluster_row = F,
                   #display_numbers = matrix2, 
                   fontsize_number = 10, number_color = "black", border_color = NA,
                   fontsize_color = 'black', fontsize_col = 10, fontsize_row = 11,
                   color = colorRampPalette(colors = c("#DCDDDD", "#000188"))(20),
                   cellwidth = 10.5, cellheight	= 12,
                   annotation_col = annotation, 
                   annotation_colors = ann_colors,
                   angle_col = 90,
                   gaps_col = c(62),
                   fontfamily= "serif")





#Heatmap for DBs###############################################################################################
rm(list=ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)
library(tidyr)


calcu_Ngene = function(data){
  Ngene = data.frame(stringsAsFactors = F)
  
  pheno_all = unique(data$pheno)
  tissue_all = unique(data$tissue)
  for (pheno_idx in pheno_all){
    for (tissue_idx in tissue_all){
      
      subset_file = subset(data, data$tissue == tissue_idx & data$pheno == pheno_idx)
      N = nrow(subset_file)
      Ngene.tmp = data.frame(pheno = pheno_idx, tissue = tissue_idx, N = N)
      Ngene = rbind(Ngene.tmp, Ngene)
    }
  }
  return(Ngene)
}
f2 = fread('results/signif_genes.CellToDisorder.fdr.txt', data.table = F)
f2 = f2[, c("pheno", "tissue","gene","gene_name")]


Ngene2 = calcu_Ngene(f2)
Ngene = Ngene2
Ngene = Ngene[-which((Ngene$N) == 0), ]
data = Ngene[, c("pheno","tissue","N")]

matrix0 = matrix(0,  length(unique(Ngene$pheno)), length(unique(Ngene$tissue)))
rownames(matrix0) = unique(Ngene$pheno)
colnames(matrix0) = unique(Ngene$tissue)
for (i in 1:nrow(Ngene)){
  pheno = Ngene$pheno[i]
  tissue = Ngene$tissue[i]
  matrix0[which(rownames(matrix0)==pheno), which(colnames(matrix0) == tissue)] = Ngene[i,'N']
}
matrix = matrix0
matrix = t(matrix)


trait_mapping = read.csv('gwas_meta_all.csv', header = T)
trait_mapping = trait_mapping[, c('TraitCategory1','TraitCategory2','TraitAbbr')]
annotation_col = data.frame(trait_mapping$TraitCategory2, trait_mapping$TraitCategory1)
rownames(annotation_col) = trait_mapping$TraitAbbr
names(annotation_col) = c('TraitCategory2', 'TraitCategory1')


ann_colors = list(TraitCategory1 = c(`Behavioral-cognitive phenotype` = "#E7298A", 
                                     `Psychiatric disorder` = "#66A61E",  `Neurological disorder` = "#D95F02"))
a = pheatmap::pheatmap(matrix, 
                       fontsize_number = 12, number_color = "black",
                       fontsize_color = 'black', fontsize_col = 10, fontsize_row = 11,
                       color = colorRampPalette(colors = c("white","#FA8072","red"))(10),
                       annotation_col = annotation_col, annotation_colors = ann_colors, 
                       fontfamily= "serif")

tree_col0 = unique(Ngene$pheno)[ a[["tree_col"]][["order"]] ]
annotation_col = subset(annotation_col, row.names(annotation_col)%in% unique(Ngene$pheno))
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


matrix = matrix[, match(tree_col, colnames(matrix))]


tree_row = unique(Ngene$tissue)[ a[["tree_row"]][["order"]] ]
matrix = matrix[c("Inhibitory.neurons","OPCs...COPs","Microglia","Endothelial.cells",
                    "Pericytes","Excitatory.neurons","Oligodendrocytes","Astrocytes"), ]


library(cols4all)
c4a_palettes("div")
mycol <- c4a("rainbow_wh_rd", 25)
c4a_plot(mycol)
set.seed(12)
mycol <- sample(mycol, 15)
ann_colors = list(`Trait_Type_col` = c(`Behavioral-cognitive phenotype` = "#E7298A", 
                                       `Psychiatric disorder` = "#66A61E",  `Neurological disorder` = "#D95F02"),
                  `Trait_Type_row` = c(`Brain volume` = "#7EAED8",  `White matter microstructure` ="#CEA07A"))
ann_colors = list(`TraitCategory1` = c(`Behavioral-cognitive phenotype` =mycol[1], 
                                       `Psychiatric disorder`=mycol[2],
                                       `Neurological disorder`=mycol[3]))
annotation_col = annotation_col[-which(names(annotation_col) == 'TraitCategory2')]

matrix2 = matrix
matrix2[which(matrix2 == 0)] = '' 
pheatmap::pheatmap(matrix, 
                   cluster_col = F, cluster_row = F,
                   fontsize_number = 10, number_color = "black", border_color = NA,
                   fontsize_color = 'black', fontsize_col = 11, fontsize_row = 11,
                   cellwidth = 12, cellheight	= 12,
                   color = colorRampPalette(colors = c("#DCDDDD", "#000188"))(20),
                   annotation_col = annotation_col, 
                   annotation_colors = ann_colors,
                   angle_col = 90,
                   gaps_col = c(10, 18), border_col = 'white',
                   fontfamily= "serif") 

