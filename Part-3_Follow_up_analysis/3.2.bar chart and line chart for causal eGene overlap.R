#number of shared cell-specific causal eGenes between each pair of phenotypes####################################
rm(list=ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)
library(stats)

f1 = fread('results/signif_genes.CellToIDP.fdr.txt', data.table = F)
f1 = f1[, c("pheno", "tissue","gene","gene_name","gene_biotype")]
f1$trait_type = 'Brain regional volume'  #Brain regional volume
f1$trait_type[grepl('FA', f1$pheno)] = 'White matter microstructure'

f2 = fread('results/signif_genes.CellToDisorder.fdr.txt', data.table = F)
f2 = f2[, c("pheno", "tissue","gene","gene_name","gene_biotype")]
trait_mapping = read.csv('gwas_meta_all.csv', header = T)
trait_mapping = trait_mapping[, c('TraitCategory1', 'TraitAbbr')]
f2$trait_type = trait_mapping[match(f2$pheno, trait_mapping$TraitAbbr), 'TraitCategory1']

dat = rbind(f1, f2)
trait_all = unique(dat$pheno)
cell_all = unique(dat$tissue)


causal_genes = data.frame(stringsAsFactors = F)
for (trait in trait_all){
  for (cell in cell_all){
    sub_dat = subset(dat, dat$pheno == trait & dat$tissue == cell)
    if (nrow(sub_dat) == 0) next
    causal_gene_list = paste(sub_dat$gene, collapse=',')
    causal_genes.tmp = data.frame(trait = trait, cell = cell, causal_gene = causal_gene_list, 
                                  trait_type = sub_dat$trait_type[1], stringsAsFactors = F)
    causal_genes = rbind(causal_genes, causal_genes.tmp)
    
    cat(trait, 'and', cell, 'is OK.', '\n')
  }
}


result = data.frame(stringsAsFactors = F)
for (i in 1:nrow(causal_genes)){
  trait1 = causal_genes[i,'trait']
  trait1_type = causal_genes[i,'trait_type']
  cell1 = causal_genes[i,'cell']
  gene_list1 = unlist(strsplit(causal_genes[i,'causal_gene'], split = ","))
  n1 = length(gene_list1)

  for (j in 1:nrow(causal_genes)){
    trait2 = causal_genes[j,'trait']
    trait2_type = causal_genes[j,'trait_type']
    cell2 = causal_genes[j,'cell']
    
    if (trait1 == trait2 | cell1 != cell2) next
    gene_list2 = unlist(strsplit(causal_genes[j,'causal_gene'], split = ","))
    n2 = length(gene_list2)
    
    overlap_gene = intersect(gene_list1, gene_list2)
    
    result.tmp = data.frame(trait1 = trait1, trait1_type = trait1_type, trait2 = trait2, trait2_type = trait2_type, cell = cell1, 
                            gene_list1 = paste(gene_list1,collapse=','), n1 = n1,
                            gene_list2 = paste(gene_list2, collapse=','), n2 = n2, 
                            overlap_gene = paste(overlap_gene, collapse=','), n_overlap = length(overlap_gene))
    result = rbind(result, result.tmp)
    cat(i, j, 'is OK.','\n')
  }
}


result$n1_prop = result$n_overlap / result$n1  


table = data.frame(TraitCategory1 = 'White matter microstructure', TraitAbbr = unique(result[result$trait1_type=='White matter microstructure','trait1']))
table = rbind(trait_mapping, table)

table$`Brain regional volume` = NA

#types = c('Brain regional volume','White matter microstructure','Behavioral-cognitive phenotype',
#          'Neurological disorder','Psychiatric disorder')
types = c('Brain regional volume')
for (i in 1:nrow(table)){
  for (type in types){
    trait = table[i,'TraitAbbr']

    sub_result = subset(result, result$trait2 == trait & result$trait1_type == type)
    #mean of shared proportion
    table[i,type] = mean(sub_result$n1_prop) 
  }
}

table$causal_count = NA
for (i in 1:nrow(table)){
  trait_idx = table[i, 'TraitAbbr']
  sub_causal_genes = subset(causal_genes, causal_genes$trait == trait_idx)
  gene_list = paste( sub_causal_genes$causal_gene, collapse=',')
  gene_list = unlist(strsplit(gene_list, split = ","))
  count = length(gene_list)
  
  table[i, 'causal_count'] = count
}
table = table[, c("TraitCategory1","TraitAbbr","Brain regional volume","causal_count")]
write.table(table, '坚果云/shared perc between BRV and other traits.txt', sep='\t', col.names = T, row.names = F, quote = F)

a = subset(table, table$TraitCategory1 == 'White matter microstructure'); max(a$`Brain regional volume`); min(a$`Brain regional volume`)
a = subset(table, table$TraitCategory1 == 'Psychiatric disorder'); max(a$`Brain regional volume`); min(a$`Brain regional volume`)
a = subset(table, table$TraitCategory1 == 'Behavioral-cognitive phenotype'); max(a$`Brain regional volume`); min(a$`Brain regional volume`)



#bar chart & line chart (BRV & DB)#############
plot = subset(table, table$TraitCategory1 %in% c("Behavioral-cognitive phenotype","Neurological disorder","Psychiatric disorder"))
plot = plot[order(plot$TraitCategory1, plot$`Brain regional volume`, decreasing = T),]
plot$TraitAbbr = factor(plot$TraitAbbr, levels = unique(plot$TraitAbbr))

library(ggplot2)
values=c(rep('#4AB74C',10), rep('#E44D00',8), rep('#0069E9',8))
ggplot(data=plot, aes(x = TraitAbbr, y = `Brain regional volume`))+
  geom_bar(stat="identity", width=0.5, position='dodge',fill = values)+
  #scale_fill_manual(values=c(rep('#999999',8), rep('#E69F00',8), rep('blue',10)))+
  ylim(0, max(plot$`Brain regional volume`)+0.02)+
  geom_text(label=paste(round(plot$`Brain regional volume` *100,2),'%',sep=''), color="black", size=3.5, position=position_dodge(0.5),vjust=-0.5)+
  theme_minimal() +
  ggtitle("Percentage sharing of cell type-specific genes between brain regional volume group with disorders and behaviors") + 
  ylab("Percentage sharing") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(data=plot, aes(x = TraitAbbr, y = causal_count)) + 
  geom_point(size=4, shape=20) +
  geom_line(group = 1, color="blue", linewidth=1) +
  geom_text(label=plot$causal_count, color="black", size=3.5, position=position_dodge(0.5),vjust=-0.5)+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



##bar chart & line chart (BRV & WMMP)####
plot = subset(table, table$TraitCategory1 %in% c("White matter microstructure"))
plot = plot[order(plot$TraitCategory1, plot$`Brain regional volume`, decreasing = T),]
plot$TraitAbbr = factor(plot$TraitAbbr, levels = unique(plot$TraitAbbr))

library(ggplot2)
ggplot(data=plot, aes(x = TraitAbbr, y = `Brain regional volume`))+
  geom_bar(stat="identity", width=0.5, position='dodge', fill = '#C58CA9')+
  #scale_fill_manual(values=c(rep('#999999',8), rep('#E69F00',8), rep('blue',10)))+
  ylim(0, max(plot$`Brain regional volume`)+0.02)+
  geom_text(label=paste(round(plot$`Brain regional volume` *100,2),'%',sep=''), color="black", size=3.5, position=position_dodge(0.5),vjust=-0.5)+
  theme_minimal() +
  ggtitle("Percentage sharing of cell type-specific genes between brain regional volume group and white matter tracts") +
  ylab("Percentage sharing") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(data=plot, aes(x = TraitAbbr, y = causal_count)) + 
  geom_point(size=4, shape=20) +
  geom_line(group = 1, color="blue", linewidth=1) +
  geom_text(label=plot$causal_count, color="black", size=3.5, position=position_dodge(0.5),vjust=-0.5)+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
