rm(list=ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)
library(tidyr)


f1 = fread('results/signif_genes.CellToDisorder.fdr.txt', data.table = F)
f1 = f1[, c("pheno", "tissue","gene","gene_name","gene_biotype")]
f1 = subset(f1, f1$pheno %in% c('SCZ','BIP','ASD'))
f1$pheno_gene = paste(f1$pheno, f1$gene, sep = '_')


#PsychEncode DGE---
dge = fread('PsychEncode_DGE.txt', data.table = F)
dge = dge[, c('ensembl_gene_id', 'ASD.fdr', 'SCZ.fdr', 'BD.fdr')]

scz_gene = subset(dge, dge$SCZ.fdr < 0.05)
scz_gene$pheno = 'SCZ'
scz_gene = scz_gene[, c('pheno', 'ensembl_gene_id')]
bip_gene = subset(dge, dge$BD.fdr < 0.05)
bip_gene$pheno = 'BIP'
bip_gene = bip_gene[, c('pheno', 'ensembl_gene_id')]
asd_gene = subset(dge, dge$ASD.fdr < 0.05)
asd_gene$pheno = 'ASD'
asd_gene = asd_gene[, c('pheno', 'ensembl_gene_id')]

dge = Reduce(rbind, list(scz_gene, bip_gene, asd_gene))
dge$pheno_gene = paste(dge$pheno, dge$ensembl_gene_id, sep = '_')


#PsychEncode DTE---
dte = fread('PsychEncode_DTE.txt', data.table = F)
dte = dte[, c('ensembl_gene_id', 'ASD.fdr', 'SCZ.fdr', 'BD.fdr')]

scz_gene = subset(dte, dte$SCZ.fdr < 0.05)
scz_gene$pheno = 'SCZ'
scz_gene = scz_gene[, c('pheno', 'ensembl_gene_id')]
bip_gene = subset(dte, dte$BD.fdr < 0.05)
bip_gene$pheno = 'BIP'
bip_gene = bip_gene[, c('pheno', 'ensembl_gene_id')]
asd_gene = subset(dte, dte$ASD.fdr < 0.05)
asd_gene$pheno = 'ASD'
asd_gene = asd_gene[, c('pheno', 'ensembl_gene_id')]

dte = Reduce(rbind, list(scz_gene, bip_gene, asd_gene))
dte$pheno_gene = paste(dte$pheno, dte$ensembl_gene_id, sep = '_')


f1$If_repl = ifelse(f1$pheno_gene %in% c(dge$pheno_gene, dte$pheno_gene), 'Y','N')
write.table(f1, 'results/signif_genes.repPsychEncode.CellToDisorder.fdr.txt', sep="\t", row.names=F, col.names=T, quote=F)


###overlap
length(dge$pheno_gene)
length(dte$pheno_gene)
length(unique(c(dge$pheno_gene), dte$pheno_gene))
##dge + dte
length(unique(c(dte$pheno_gene, dge$pheno_gene))); length(unique(f1$pheno_gene[which(f1$pheno_gene %in% c(dge$pheno_gene, dte$pheno_gene))] ));length(unique(f1$pheno_gene)) 
q = length(unique(f1$pheno_gene[which(f1$pheno_gene %in% c(dge$pheno_gene, dte$pheno_gene))] ))
m = length(unique(f1$pheno_gene))
n = 26861*3 - m
k = length(unique(c(dte$pheno_gene, dge$pheno_gene)))
phyper(q-1, m, n, k, lower.tail=F)


#29.5%, 23 out of 78
length(unique(f1$pheno_gene[which(f1$pheno_gene %in% c(dge$pheno_gene, dte$pheno_gene))] )) / length(unique(f1$pheno_gene)) 

