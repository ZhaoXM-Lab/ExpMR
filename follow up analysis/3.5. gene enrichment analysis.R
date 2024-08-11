rm(list = ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)
library(tidyverse)
library(AnnotationHub)
library(org.Hs.eg.db) 
library(clusterProfiler)
library(dplyr)
library(ggplot2)


f1 = fread('results/signif_genes.CellToIDP.fdr.txt', data.table = F)
f2 = fread('results/signif_genes.CellToDisorder.fdr.txt', data.table = F)
f1$tissue_gene = paste(f1$tissue, f1$gene, sep='.')
f2$tissue_gene = paste(f2$tissue, f2$gene, sep='.')
trait_mapping = read.csv('gwas_meta_all.csv', header = T)
trait_mapping = trait_mapping[, c('TraitCategory1', 'TraitAbbr')]
PD = trait_mapping[trait_mapping$TraitCategory1 == 'Psychiatric disorder', 'TraitAbbr']
ND = trait_mapping[trait_mapping$TraitCategory1 == 'Neurological disorder', 'TraitAbbr']
BP = trait_mapping[trait_mapping$TraitCategory1 == 'Behavioral-cognitive phenotype', 'TraitAbbr']

set1 = f1[!grepl('FA', f1$pheno), ]  #Brain regional volume
set2 = f1[grepl('FA', f1$pheno), ]  #White matter microstructure
set3 = subset(f2, f2$pheno %in% PD)  #Psychiatric disorder
set4 = subset(f2, f2$pheno %in% ND)  #Neurological disorder
set5 = subset(f2, f2$pheno %in% BP)  #Behavioral-cognitive disorder


p = list()
p <- list(`set1 & set1` = names(table(set1$tissue_gene))[which(table(set1$tissue_gene) > 1)], #shared in at least two traits in BRV
          `set2 & set2` = names(table(set2$tissue_gene))[which(table(set2$tissue_gene) > 1)],
          `set3 & set3` = names(table(set3$tissue_gene))[which(table(set3$tissue_gene) > 1)],
          `set4 & set4` = names(table(set4$tissue_gene))[which(table(set4$tissue_gene) > 1)],
          `set5 & set5` = names(table(set5$tissue_gene))[which(table(set5$tissue_gene) > 1)],
          
          `set1 & set2` = Reduce(intersect, list(set1$tissue_gene, set2$tissue_gene)), 
          `set1 & set3` = Reduce(intersect, list(set1$tissue_gene, set3$tissue_gene)),
          `set1 & set4` = Reduce(intersect, list(set1$tissue_gene, set4$tissue_gene)),
          `set1 & set5` = Reduce(intersect, list(set1$tissue_gene, set5$tissue_gene)),
          `set2 & set3` = Reduce(intersect, list(set2$tissue_gene, set3$tissue_gene)),
          `set2 & set4` = Reduce(intersect, list(set2$tissue_gene, set4$tissue_gene)),
          `set2 & set5` = Reduce(intersect, list(set2$tissue_gene, set5$tissue_gene)),
          `set3 & set4` = Reduce(intersect, list(set3$tissue_gene, set4$tissue_gene)),
          `set3 & set5` = Reduce(intersect, list(set3$tissue_gene, set5$tissue_gene)),
          `set4 & set5` = Reduce(intersect, list(set4$tissue_gene, set5$tissue_gene)),
          `set1 & set2 & set3` = Reduce(intersect, list(set1$tissue_gene,set2$tissue_gene,set3$tissue_gene)),
          `set1 & set2 & set4` = Reduce(intersect, list(set1$tissue_gene,set2$tissue_gene,set4$tissue_gene)),
          `set1 & set2 & set5` = Reduce(intersect, list(set1$tissue_gene,set2$tissue_gene,set5$tissue_gene)),
          `set2 & set3 & set4` = Reduce(intersect, list(set2$tissue_gene,set3$tissue_gene,set4$tissue_gene)),
          `set2 & set3 & set5` = Reduce(intersect, list(set2$tissue_gene,set3$tissue_gene,set5$tissue_gene)),
          `set3 & set4 & set5` = Reduce(intersect, list(set3$tissue_gene,set4$tissue_gene,set5$tissue_gene)),
          `set1 & set2 & set3 & set4` = Reduce(intersect, list(set1$tissue_gene,set2$tissue_gene,set3$tissue_gene,set4$tissue_gene)),
          `set1 & set2 & set3 & set5` = Reduce(intersect, list(set1$tissue_gene,set2$tissue_gene,set3$tissue_gene,set5$tissue_gene)),
          `set2 & set3 & set4 & set5` = Reduce(intersect, list(set2$tissue_gene,set3$tissue_gene,set4$tissue_gene,set5$tissue_gene)),
          `set1 & set2 & set3 & set4 & set5` = Reduce(intersect, list(set1$tissue_gene,set2$tissue_gene,set3$tissue_gene,set4$tissue_gene,set5$tissue_gene))
)

mapping = rbind(f1[,c("gene",'tissue_gene')], f2[, c("gene",'tissue_gene')])
mapping = mapping[!(duplicated(mapping$tissue_gene)), ]

for (i in 1:length(p)){
  p[[i]] = subset(mapping, mapping$tissue_gene %in% p[[i]]) $ gene
}


gene_info = fread('F:/类脑/1_实操文件/genexpMR/data/gene_loc/Homo_sapiens.GRCh37.87.gtf.gz', data.table = F, skip = 5)
gene_info = separate(gene_info, col = V9, 
                     into = c("gene_id", "gene_version","gene_name","gene_source","gene_biotype"),sep = "; ")
gene_info = subset(gene_info, gene_info$V3 == 'gene')
gene_info = gene_info[, c('gene_name', 'V1', 'V4', 'V5','V7', 'gene_id','gene_biotype')]
gene_info$gene_id = gsub("gene_id ", "", gene_info$gene_id)
gene_info$gene_id = gsub('"', '', gene_info$gene_id)

EG2Ensembl = toTable(org.Hs.egENSEMBL)


result = data.frame(stringsAsFactors = F)
for (i in 1:length(p)){
  result.tmp = data.frame(stringsAsFactors = F)
  
  ensembl_id = unique(p[[i]])
  #if (length(ensembl_id) <= 10) next
  id = EG2Ensembl$gene_id[match(ensembl_id, EG2Ensembl$ensembl_id)]
  
  if (length(id) == 0) next
  ego = clusterProfiler::enrichGO(OrgDb="org.Hs.eg.db", gene = id, ont = "ALL",
                                  readable= TRUE, pvalueCutoff = 1, pAdjustMethod='none')
  if (length(ego) == 0) next
  gene_set_id = ego@result[["ID"]]
  gene_set_descrip = ego@result[["Description"]]
  gene_set_count = ego@result[["Count"]]
  gene_set_p = ego@result[["pvalue"]]
  
  
  if (length(gene_set_id) == 0) {
    result.tmp = t(data.frame(c(group = names(p)[i], gene_set_id='nan', gene_set_descrip='nan', gene_set_count='nan', gene_set_p='nan'), stringsAsFactors = F))
  } else {
    result.tmp = data.frame(group = names(p)[i], gene_set_id, gene_set_descrip, gene_set_count, gene_set_p)
  }
  
  result = rbind(result, result.tmp)
  cat(i, '\n')
}

result = result[!is.na(result$gene_set_id), ]

result$group = gsub('set1', 'BRV', result$group) #Brain regional volume
result$group = gsub('set2', 'WMM',result$group)  #White matter microstructure
result$group = gsub('set3', 'PD', result$group)  #Psychiatric disorder
result$group = gsub('set4', 'ND', result$group)  #Neurological disorder
result$group = gsub('set5', 'BCP', result$group) #Behavioral-cognitive disorder

write.table(result, 'results/gsea_ALL_results.2024-03-17.txt', sep="\t", row.names=F, col.names=T, quote=F)