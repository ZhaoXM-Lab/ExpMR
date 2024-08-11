rm(list=ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)
library(stats)

f1 = fread('results/signif_genes.CellToIDP.fdr.txt', data.table = F)
f2 = fread('results/signif_genes.CellToDisorder.fdr.txt', data.table = F)

IDP_all = unique(f1$pheno)
disorder_all = unique(f2$pheno)
cell_all = unique(c(f1$tissue, f2$tissue))

result = data.frame(stringsAsFactors = F)
for (IDP in IDP_all){
  for (disorder in disorder_all){
    for (cell in cell_all){
      sub_f1 = subset(f1, f1$pheno == IDP & f1$tissue == cell)
      sub_f2 = subset(f2, f2$pheno == disorder & f2$tissue == cell)
      
      gene_list1 = sub_f1$gene_name
      gene_list2 = sub_f2$gene_name
      inter_gene = intersect(gene_list1, gene_list2)
      if (length(inter_gene) == 0) next

      q = length(inter_gene)
      m = length(gene_list1)
      n = 26861 - m
      k = length(gene_list2)
      p_hyper = phyper(q-1, m, n, k, lower.tail=F)

      result.tmp = data.frame(IDP = IDP, Disorder = disorder, Cell = cell, Causal_genes1 = paste(gene_list1,collapse=','),
                              Causal_genes2 = paste(gene_list2, collapse=','),
                              Overlap_genes = paste(inter_gene,collapse=','), p_hyper = p_hyper, stringsAsFactors = F)
      result = rbind(result, result.tmp)
      cat(IDP, 'and', disorder, 'is OK.', '\n') 
    }
  }
}
length(unique(result$IDP))
length(unique(result$Disorder))
write.table(result, 'results/overlap_genes.hyper_test.fdr.txt', sep="\t", row.names=F, col.names=T, quote=F)

