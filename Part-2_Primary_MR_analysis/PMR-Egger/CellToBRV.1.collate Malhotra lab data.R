rm(list =ls())
setwd('F:/类脑/1_实操文件/genexpMR/data/gene_expression/Malhotra_lab/original_expdata')
library(data.table)
library(tidyr)


lookup_table = fread('snp_pos.txt.gz', data.table = F)
lookup_table = lookup_table[,c('SNP','SNP_id_hg19')]
lookup_table = separate(lookup_table, SNP_id_hg19, into= c("chr","variant_pos"), sep= ":")


g1000_freq = fread('../../../gwas_summary/g1000_eur.frq', data.table = F)
lookup_table = subset(lookup_table, lookup_table$SNP %in% g1000_freq$SNP) 


gene_info = fread('F:/类脑/1_实操文件/genexpMR/data/gene_loc/Homo_sapiens.GRCh37.87.gtf.gz', data.table = F, skip = 5)
gene_info = separate(gene_info, col = V9, 
                     into = c("gene_id", "gene_version","gene_name","gene_source","gene_biotype"),sep = "; ")
gene_info = subset(gene_info, gene_info$V3 == 'gene')
gene_info = gene_info[, c('gene_name', 'V1', 'V4', 'V5','V7', 'gene_id','gene_biotype')]
gene_info$gene_id = gsub("gene_id ", "", gene_info$gene_id)
gene_info$gene_id = gsub('"', '', gene_info$gene_id)

#cell types
cell_types = c('Astrocytes','Endothelial.cells','Excitatory.neurons','Inhibitory.neurons',
               'Microglia','Oligodendrocytes','OPCs...COPs','Pericytes')


for (i in 1:length(cell_types)) {
  cell_type = cell_types[i]
  
  data = data.frame(stringsAsFactors = F)
  for (chr in 1:22){
    file_name = paste(cell_type, '.', chr, '.gz', sep='')
    data.tmp = fread(file_name, data.table = F)
    data = rbind(data, data.tmp)
    cat(i,', chr', chr, '\n')
  }
  names(data) = c('Gene_id', 'SNP_id','Distance to TSS', 'Nominal p-value', 'Beta')
  
  data = data[, c("SNP_id","Gene_id","Nominal p-value","Beta")]
  data = subset(data, data$`Nominal p-value` < 5e-8)
  data$Zscore = sign(data$Beta) * abs(qnorm(0.5* data$`Nominal p-value`))
  
  lookup_table_subset = subset(lookup_table, lookup_table$SNP %in% data$SNP_id) 
  
  data_mr = merge(data, lookup_table_subset, by.x = 'SNP_id',by.y = 'SNP', all = F)
  data_mr = data_mr[, c('Gene_id','SNP_id','chr','variant_pos','Zscore')]
  names(data_mr) = c('Gene','SNP','chr','bp','Zscore')
  data_mr = data_mr[order(data_mr$chr),]
  

  data_mr =  separate(data_mr, Gene, into= c("gene_name","gene_id"), sep= "_", remove = F)
  data_mr = subset(data_mr,data_mr$gene_id %in% gene_info$gene_id)  
  data_mr = data_mr[,-which(names(data_mr) %in% c("gene_name","gene_id"))]
  

  outname = paste('../mr_format/PMR-Egger/', cell_type,'.signif_variant_gene_pairs_5e-8.txt',sep='')
  write.table(data_mr, outname, sep='\t', row.names=F,quote=F)
  

  selected_snp = unique(data_mr$SNP)
  outname = paste('../allele/', cell_type, '.allele', sep='')
  write.table(selected_snp, outname, sep='\t', row.names=F, col.names=F, quote=F)
}
