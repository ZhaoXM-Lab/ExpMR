rm(list=ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)
library(stats)
library(tidyr)

f1 = fread('results/signif_genes.CellToIDP.fdr.txt', data.table = F)
f1 = f1[, c("pheno", "tissue","gene","gene_name","gene_biotype")]
f1$trait_type = 'IDP'

f1 = f1[grepl('FA', f1$pheno), ]  #WM

f1$pheno = gsub('[.]',' ', f1$pheno)
f1$tissue = gsub('Inhibitory.neurons','Inhibitory neurons', f1$tissue)
f1$tissue = gsub('OPCs...COPs','OPCs', f1$tissue)
f1$tissue = gsub('Excitatory.neurons','Excitatory neurons', f1$tissue)
f1$tissue = gsub('Endothelial.cells','Endothelial cells', f1$tissue)

f2 = fread('results/signif_genes.CellToDisorder.fdr.txt', data.table = F)
f2 = f2[, c("pheno", "tissue","gene","gene_name","gene_biotype")]
trait_mapping = read.csv('gwas_meta_all.csv', header = T)
trait_mapping = trait_mapping[, c('TraitCategory1', 'TraitAbbr')]
f2$trait_type = trait_mapping[match(f2$pheno, trait_mapping$TraitAbbr), 'TraitCategory1']

f2$tissue = gsub('Inhibitory.neurons','Inhibitory neurons', f2$tissue)
f2$tissue = gsub('OPCs...COPs','OPCs', f2$tissue)
f2$tissue = gsub('Excitatory.neurons','Excitatory neurons', f2$tissue)
f2$tissue = gsub('Endothelial.cells','Endothelial cells', f2$tissue)

dat = rbind(f1, f2)
trait_all = unique(dat$pheno)
cell_all = unique(dat$tissue)


causal_genes = data.frame(stringsAsFactors = F)
for (trait in trait_all){
  for (cell in cell_all){
    sub_dat = subset(dat, dat$pheno == trait & dat$tissue == cell)
    if (nrow(sub_dat) == 0) next
    causal_gene_list = paste(sub_dat$gene_name, collapse=',')
    #cell-gene
    cell_causal_gene_list = paste(cell, sub_dat$gene_name, sep = ' - ')
    cell_causal_gene_list = paste(cell_causal_gene_list, collapse = ',')
    
    causal_genes.tmp = data.frame(trait = trait, cell = cell, causal_gene = causal_gene_list, 
                                  cell_causal_gene = cell_causal_gene_list, 
                                  trait_type = sub_dat$trait_type[1], stringsAsFactors = F)
    causal_genes = rbind(causal_genes, causal_genes.tmp)
    
    cat(trait, 'and', cell, 'is OK.', '\n')
  }
}


cell_gene_all = paste(causal_genes$cell_causal_gene, collapse = ',')
cell_gene_all = unlist(strsplit(cell_gene_all, split = ','))
cell_gene_all = unique(cell_gene_all)

table = matrix(data=0, nrow = length(cell_gene_all), ncol = length(trait_all))
table = data.frame(table, stringsAsFactors = F)
row.names(table) = cell_gene_all
names(table) = trait_all

for (i in 1:nrow(table)){
  cell_gene = row.names(table)[i]
  sub_causal_genes = causal_genes[grep(cell_gene, causal_genes$cell_causal_gene), ]
  traits = sub_causal_genes$trait
  table[i, match(traits, names(table))] = 1
}
table = cbind(table, Total = rowSums(table))
table = table[,c(ncol(table), 1:(ncol(table)-1))]

table = subset(table, table$Total > 1)


trait_mapping = dat[,c('pheno', 'trait_type')]
trait_mapping = trait_mapping[!(duplicated(trait_mapping$pheno)), ]


#IDP
IDP_all = subset(trait_mapping, trait_mapping$trait_type =='IDP') $pheno
table_IDP = table[, names(table) %in% IDP_all]
table_IDP = cbind(table_IDP, Total_IDP = rowSums(table_IDP))
table_IDP = table_IDP[,c(ncol(table_IDP), 1:(ncol(table_IDP)-1))]
#Behavioral-cognitive phenotype
BCP_all = subset(trait_mapping, trait_mapping$trait_type =='Behavioral-cognitive phenotype') $pheno
table_BCP = table[, names(table) %in% BCP_all]
table_BCP = cbind(table_BCP, Total_BCP = rowSums(table_BCP))
table_BCP = table_BCP[,c(ncol(table_BCP), 1:(ncol(table_BCP)-1))]
#Neurological disorder
ND_all = subset(trait_mapping, trait_mapping$trait_type =='Neurological disorder') $pheno
table_ND = table[, names(table) %in% ND_all]
table_ND = cbind(table_ND, Total_ND = rowSums(table_ND))
table_ND = table_ND[,c(ncol(table_ND), 1:(ncol(table_ND)-1))]
#Psychiatric disorder
PD_all = subset(trait_mapping, trait_mapping$trait_type =='Psychiatric disorder') $pheno
table_PD = table[, names(table) %in% PD_all]
table_PD = cbind(table_PD, Total_PD = rowSums(table_PD))
table_PD = table_PD[,c(ncol(table_PD), 1:(ncol(table_PD)-1))]


table_IDP_BCP_ND_PD = cbind(table_IDP, table_BCP)
table_IDP_BCP_ND_PD = cbind(table_IDP_BCP_ND_PD, table_ND)
table_IDP_BCP_ND_PD = cbind(table_IDP_BCP_ND_PD, table_PD)
col_total = which(names(table_IDP_BCP_ND_PD) %in% c('Total_IDP','Total_BCP','Total_ND','Total_PD'))
table_IDP_BCP_ND_PD = table_IDP_BCP_ND_PD[,c(col_total, setdiff(1:ncol(table_IDP_BCP_ND_PD), col_total))]


plot = subset(table_IDP_BCP_ND_PD, table_IDP_BCP_ND_PD$Total_IDP > 0 &  table_IDP_BCP_ND_PD$Total_PD > 0)
plot = plot[,names(plot) %in% c(names(table_IDP), names(table_PD))]
plot = rbind(plot, Total = colSums(plot))
plot = plot[,-which(plot['Total',]==0)]

plot = plot[, c(-1,-2)]
plot = plot[-nrow(plot), ]


for (i in 1:ncol(plot)){
  plot[,i] = gsub(1, i, plot[,i])
  plot[,i] = as.numeric(plot[,i])
  plot[,i][which(plot[,i] == 0 )] = NA
  plot[,i] = as.numeric(plot[,i])
}
plot$value1 = NA
plot$value2 = NA
plot = plot[,c(ncol(plot)-1, ncol(plot), 1:(ncol(plot)-2))]
for (i in 1:nrow(plot)){
  plot[i,'value1'] = min(plot[i,1:ncol(plot)], na.rm = T)
  plot[i,'value2'] = max(plot[i,1:ncol(plot)], na.rm = T)
}
plot$group = row.names(plot)
plot = plot[,c(ncol(plot), 1:(ncol(plot)-1))]
plot$group = factor(plot$group, levels = unique(plot$group))
for (i in 2:ncol(plot)) plot[,i] = as.numeric(plot[,i])


library(ggplot2)
p <- ggplot(plot) +
  geom_segment(aes(x=group, xend=group, y=value1, yend=value2), color='grey', linewidth=1) +
  theme_minimal() +
  xlab("Cell-specific Genes") +
  ylab("Traits") +
  scale_y_continuous(limits=c(1,ncol(plot)-3), breaks=seq(1,ncol(plot)-3,1),labels = c(names(plot)[4:ncol(plot)])) +
  theme(axis.text.x=element_text(angle=90, size=14, vjust = 0.25, hjust = 1,color = 'black'),
        axis.text.y=element_text(size=14, color = 'black'), 
        plot.title = element_text(size=18, color = 'black')) +
  labs(x = NULL, y = NULL) +
  ggtitle("Shared genes in specific cell types between white matter microstructure and psychiatric disorders") + 
  coord_flip()

for (i in c(4:ncol(plot))){
  if (!(names(plot)[i] %in% PD_all)) {
    p <- p + 
      eval(parse(text = paste("geom_point( aes(x=group, y= plot[,",i,"]),shape = 21,color='black',fill='#0069E9', size=4)"))) 
  }
  if (names(plot)[i] %in% PD_all){
    p <- p + 
      eval(parse(text = paste("geom_point( aes(x=group, y= plot[,",i,"]),shape = 21,color='black',fill='#FEB500', size=4)")))
  } 
}
p


color_cell_mapping = c(`Astrocytes` = '#F0027F', `Inhibitory neurons`='#FDC086', `Excitatory neurons`='#BEAED4', 
                       Oligodendrocytes = '#386CB0', OPCs = '#4DC97F', Microglia = '#FFFF99', Pericytes = '#BF5B17')
cell_set = separate(plot, col = group, into = c('cell', 'gene'), sep = ' - ')
cell_set = cell_set$cell
my_colors = color_cell_mapping[match(cell_set, names(color_cell_mapping))]
barplot(rep(1,times=length(my_colors)), col=my_colors, border=my_colors, axes=FALSE, horiz=T)



plot = subset(table_IDP_BCP_ND_PD, table_IDP_BCP_ND_PD$Total_IDP > 0 & table_IDP_BCP_ND_PD$Total_ND > 0)
plot = plot[,names(plot) %in% c(names(table_IDP), names(table_ND))]
plot = rbind(plot, Total = colSums(plot)) 
plot = plot[,-which(plot['Total',]==0)]
plot = plot[, c(-1,-2)]
plot = plot[-nrow(plot), ]


for (i in 1:ncol(plot)){
  plot[,i] = gsub(1, i, plot[,i])
  plot[,i] = as.numeric(plot[,i])
  plot[,i][which(plot[,i] == 0 )] = NA
  plot[,i] = as.numeric(plot[,i])
}

plot$value1 = NA
plot$value2 = NA
plot = plot[,c(ncol(plot)-1, ncol(plot), 1:(ncol(plot)-2))]
for (i in 1:nrow(plot)){
  plot[i,'value1'] = min(plot[i,1:ncol(plot)], na.rm = T)
  plot[i,'value2'] = max(plot[i,1:ncol(plot)], na.rm = T)
}
plot$group = row.names(plot)
plot = plot[,c(ncol(plot), 1:(ncol(plot)-1))]
plot$group = factor(plot$group, levels = unique(plot$group))
for (i in 2:ncol(plot)) plot[,i] = as.numeric(plot[,i])


library(ggplot2)
p <- ggplot(plot) +
  geom_segment(aes(x=group, xend=group, y=value1, yend=value2), color='grey',linewidth=1) +
  theme_minimal() +
  xlab("Cell-specific Genes") +
  ylab("Traits") +
  scale_y_continuous(limits=c(1,ncol(plot)-3), breaks=seq(1,ncol(plot)-3,1),labels = c(names(plot)[4:ncol(plot)])) +
  theme(axis.text.x=element_text(angle=90, size=14, vjust = 0.25, hjust = 1,color = 'black'),
        axis.text.y=element_text(size=14, color = 'black'), 
        plot.title = element_text(size=18, color = 'black')) +
  labs(x = NULL, y = NULL) + 
  ggtitle("Shared genes in specific cell types between white matter microstructure and neurological disorders") + 
  coord_flip()

for (i in c(4:ncol(plot))){
  if (!(names(plot)[i] %in% ND_all)) {
    p <- p + 
      eval(parse(text = paste("geom_point( aes(x=group, y= plot[,",i,"]),shape = 21,color='black',fill='#0069E9', size=4)"))) 
  }
  if (names(plot)[i] %in% ND_all){
    p <- p + 
      eval(parse(text = paste("geom_point( aes(x=group, y= plot[,",i,"]),shape = 21,color='black',fill='#E44D00', size=4)")))
  } 
}
p


color_cell_mapping = c(`Astrocytes` = '#F0027F', `Inhibitory neurons`='#FDC086', `Excitatory neurons`='#BEAED4', 
                       Oligodendrocytes = '#386CB0', OPCs = '#4DC97F', Microglia = '#FFFF99', Pericytes = '#BF5B17')
cell_set = separate(plot, col = group, into = c('cell', 'gene'), sep = ' - ')
cell_set = cell_set$cell
my_colors = color_cell_mapping[match(cell_set, names(color_cell_mapping))]
barplot(rep(1,times=length(my_colors)), col=my_colors, border=my_colors, axes=FALSE, horiz=T)



plot = subset(table_IDP_BCP_ND_PD, table_IDP_BCP_ND_PD$Total_IDP > 0 & table_IDP_BCP_ND_PD$Total_BCP > 0)
plot = plot[,names(plot) %in% c(names(table_IDP), names(table_BCP))]
plot = rbind(plot, Total = colSums(plot))
plot = plot[,-which(plot['Total',]==0)]

plot = plot[, c(-1,-2)]
plot = plot[-nrow(plot), ]


for (i in 1:ncol(plot)){
  plot[,i] = gsub(1, i, plot[,i])
  plot[,i] = as.numeric(plot[,i])
  plot[,i][which(plot[,i] == 0 )] = NA
  plot[,i] = as.numeric(plot[,i])
}

plot$value1 = NA
plot$value2 = NA
plot = plot[,c(ncol(plot)-1, ncol(plot), 1:(ncol(plot)-2))]
for (i in 1:nrow(plot)){
  plot[i,'value1'] = min(plot[i,1:ncol(plot)], na.rm = T)
  plot[i,'value2'] = max(plot[i,1:ncol(plot)], na.rm = T)
}
plot$group = row.names(plot)
plot = plot[,c(ncol(plot), 1:(ncol(plot)-1))]
plot$group = factor(plot$group, levels = unique(plot$group))
for (i in 2:ncol(plot)) plot[,i] = as.numeric(plot[,i])


library(ggplot2)
p <- ggplot(plot) +
  geom_segment(aes(x=group, xend=group, y=value1, yend=value2), color='grey',linewidth=1) +
  theme_minimal() +
  xlab("Cell-specific Genes") +
  ylab("Traits") +
  scale_y_continuous(limits=c(1,ncol(plot)-3), breaks=seq(1,ncol(plot)-3,1),labels = c(names(plot)[4:ncol(plot)])) +
  theme(axis.text.x=element_text(angle=90, size=14, vjust = 0.25, hjust = 1,color = 'black'),
        axis.text.y=element_text(size=14, color = 'black'), 
        plot.title = element_text(size=18, color = 'black')) +
  labs(x = NULL, y = NULL) +
  ggtitle("Shared genes in specific cell types between white matter microstructure and behavioral-cognitive phenotypes") + 
  coord_flip()

for (i in c(4:ncol(plot))){
  if (!(names(plot)[i] %in% BCP_all)) {

    p <- p + 
      eval(parse(text = paste("geom_point( aes(x=group, y= plot[,",i,"]),shape = 21,color='black',fill='#0069E9', size=4)"))) #WM
  }
  if (names(plot)[i] %in% BCP_all){
    p <- p + 
      eval(parse(text = paste("geom_point( aes(x=group, y= plot[,",i,"]),shape = 21,color='black',fill='#C58CA9', size=4)")))
  } 
}
p


color_cell_mapping = c(`Astrocytes` = '#F0027F', `Inhibitory neurons`='#FDC086', `Excitatory neurons`='#BEAED4', 
                       Oligodendrocytes = '#386CB0', OPCs = '#4DC97F', Microglia = '#FFFF99', Pericytes = '#BF5B17')
cell_set = separate(plot, col = group, into = c('cell', 'gene'), sep = ' - ')
cell_set = cell_set$cell
my_colors = color_cell_mapping[match(cell_set, names(color_cell_mapping))]
barplot(rep(1,times=length(my_colors)), col=my_colors, border=my_colors, axes=FALSE, horiz=T)
