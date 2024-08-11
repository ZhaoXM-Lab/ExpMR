rm(list =ls())
setwd('F:/类脑/1_实操文件/genexpMR/data/gene_expression/')
library(data.table)
library(tidyverse)
library(TwoSampleMR)


lookup_table = fread('Malhotra_lab/original_expdata/snp_pos.txt.gz', data.table = F)
lookup_table = lookup_table[,c('SNP','SNP_id_hg38')]
lookup_table = separate(lookup_table, SNP_id_hg38, into= c("chr","variant_pos"), sep= ":")


g1000_freq = fread('../gwas_summary/g1000_eur.frq', data.table = F)
lookup_table = subset(lookup_table, lookup_table$SNP %in% g1000_freq$SNP) 


gene_info = fread('F:/类脑/1_实操文件/genexpMR/data/gene_loc/Homo_sapiens.GRCh37.87.gtf.gz', data.table = F, skip = 5)
gene_info = separate(gene_info, col = V9, 
                     into = c("gene_id", "gene_version","gene_name","gene_source","gene_biotype"),sep = "; ")
gene_info = subset(gene_info, gene_info$V3 == 'gene')
gene_info = gene_info[, c('gene_name', 'V1', 'V4', 'V5','V7', 'gene_id','gene_biotype')]
gene_info$gene_id = gsub("gene_id ", "", gene_info$gene_id)
gene_info$gene_id = gsub('"', '', gene_info$gene_id)


cell_types = c('Astrocytes','Endothelial.cells','Excitatory.neurons','Inhibitory.neurons',
               'Microglia','Oligodendrocytes','OPCs...COPs','Pericytes')


for (i in 1:length(cell_types)) {
  cell_type = cell_types[i]
  
  data = data.frame(stringsAsFactors = F)
  for (chr in 1:22){
    file_name = paste('Malhotra_lab/original_expdata/', cell_type, '.', chr, '.gz', sep='')
    data.tmp = fread(file_name, data.table = F)
    data = rbind(data, data.tmp)
    cat(i,', chr', chr, '\n')
  }
  names(data) = c('Gene_id', 'SNP_id','Distance to TSS', 'Nominal p-value', 'Beta')
  
  data = data[, c("SNP_id","Gene_id","Nominal p-value","Beta")]
  data = subset(data, data$`Nominal p-value` < 5e-8)
  data$Zscore = sign(data$Beta) * abs(qnorm(0.5 * data$`Nominal p-value`))
  data$se = data$Beta / data$Zscore
  data = data[, c("SNP_id","Gene_id","Beta","se","Nominal p-value")]
  
  lookup_table_subset = subset(lookup_table, lookup_table$SNP %in% data$SNP_id)
  
  data_mr = merge(data, lookup_table_subset, by.x = 'SNP_id',by.y = 'SNP', all = F)
  data_mr = data_mr[, c('SNP_id','Gene_id','Beta','se','Nominal p-value','chr')]
  

  g1000_freq_subset = subset(g1000_freq, g1000_freq$SNP %in% data_mr$SNP_id)
  data_mr = merge(data_mr, g1000_freq_subset[,2:5], by.x = 'SNP_id', by.y = 'SNP', all.x = T)  
  

  data_mr = data_mr[, c('SNP_id','Gene_id','A1','A2','MAF','Beta','se','Nominal p-value')]
  names(data_mr) = c('SNP','gene','effect_allele','other_allele','eaf','beta','se','pval')
  
  data_mr = separate(data_mr, col = SNP, into = c("SNP1", "SNP2"),sep = ",")
  data_mr = data_mr[, c('SNP1','gene','effect_allele','other_allele','eaf','beta','se','pval')]
  names(data_mr)[1] = 'SNP'
  
  data_mr$id = paste(data_mr$SNP, data_mr$gene)
  data_mr = data_mr[!(duplicated(data_mr$id)), ]
  data_mr = data_mr[,-which(names(data_mr)=='id')]
  

  data_mr =  separate(data_mr, gene, into= c("gene_name","gene_id"), sep= "_", remove = F)
  data_mr = subset(data_mr,data_mr$gene_id %in% gene_info$gene_id)   
  data_mr = data_mr[, -which(names(data_mr) %in% c("gene_name","gene_id"))]
  
  data_mr$samplesize = 192
  
  outname = paste('Malhotra_lab/mr_format/TwoSampleMR/', cell_type,'.signif_pairs_5e-8.txt',sep='')
  write.table(data_mr, outname, sep=' ', row.names=F,quote=F)
  
  cat(i, '\n')
}


