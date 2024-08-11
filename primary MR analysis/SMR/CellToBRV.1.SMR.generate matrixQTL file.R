#make matrixQTL file##########################################################
rm(list =ls())
setwd('F:/类脑/1_实操文件/genexpMR/data/gene_expression/')
library(data.table)
library(tidyverse)


lookup_table = fread('Malhotra_lab/original_expdata/snp_pos.txt.gz', data.table = F)
lookup_table = lookup_table[,c('SNP','SNP_id_hg19')]
lookup_table = separate(lookup_table, SNP_id_hg19, into= c("chr","variant_pos"), sep= ":")


g1000_freq = fread('../gwas_summary/g1000_eur.frq', data.table = F)
lookup_table = subset(lookup_table, lookup_table$SNP %in% g1000_freq$SNP)

cell_types = c('Astrocytes','Endothelial.cells','Excitatory.neurons','Inhibitory.neurons',
               'Microglia','Oligodendrocytes','OPCs...COPs','Pericytes')


gene_info = fread('F:/类脑/1_实操文件/genexpMR/data/gene_loc/Homo_sapiens.GRCh37.87.gtf.gz', data.table = F, skip = 5)
gene_info = separate(gene_info, col = V9, 
                     into = c("gene_id", "gene_version","gene_name","gene_source","gene_biotype"),sep = "; ")
gene_info = subset(gene_info, gene_info$V3 == 'gene')
gene_info = gene_info[, c('gene_name', 'V1', 'V4', 'V5','V7', 'gene_id','gene_biotype')]
gene_info$gene_id = gsub("gene_id ", "", gene_info$gene_id)
gene_info$gene_id = gsub('"', '', gene_info$gene_id)


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
  names(data)[2] = 'gene_id'
  data = subset(data, data$`Nominal p-value` < 5e-8)
  data$Zscore = sign(data$Beta) * abs(qnorm(0.5 * data$`Nominal p-value`))
  
  lookup_table_subset = subset(lookup_table, lookup_table$SNP %in% data$SNP_id)
  
  data_mr = merge(data, lookup_table_subset, by.x = 'SNP_id',by.y = 'SNP', all = F)
  data_mr = data_mr[, c('SNP_id','gene_id','Beta','Zscore','Nominal p-value','chr')]
  
  data_mr$FDR = 0
  data_mr = data_mr[,-which(names(data_mr)=='chr')]
  names(data_mr) = c('SNP','gene','beta','t-stat','p-value','FDR')
  

  data_mr = separate(data_mr, col = SNP, into = c("SNP1", "SNP2"),sep = ",")
  data_mr = data_mr[, c('SNP1','gene','beta','t-stat','p-value','FDR')]
  names(data_mr)[1] = 'SNP'
  

  data_mr$id = paste(data_mr$SNP, data_mr$gene)
  data_mr = data_mr[!(duplicated(data_mr$id)), ]
  data_mr = data_mr[,-which(names(data_mr)=='id')]
  

  data_mr =  separate(data_mr, gene, into= c("gene_name","gene_id"), sep= "_", remove = F)
  data_mr = subset(data_mr,data_mr$gene_id %in% gene_info$gene_id)   
  data_mr = data_mr[,-which(names(data_mr) %in% c("gene_name","gene_id"))]
  

  outname = paste('Malhotra_lab/mr_format/SMR/', cell_type,'.signif_pairs_5e-8.mateQTL.txt',sep='')
  write.table(data_mr, outname, sep=' ', row.names=F,quote=F)
  
  cat(i, '\n')
}



#Make a BESD file from Matrix eQTL output######################################################
rm(list =ls())
setwd('/public/home/yanganyi/Desktop/genexpMR/data/gene_expression/')
library(data.table)

tissues = c('Astrocytes','Endothelial.cells','Excitatory.neurons','Inhibitory.neurons',
            'Microglia','Oligodendrocytes','OPCs...COPs','Pericytes')

for (i in 1:length(tissues)){
  tissue = tissues[i]
  command = paste('/public/home/yanganyi/Desktop/Software/smr-1.3.1-linux-x86_64/smr-1.3.1', 
                  ' --eqtl-summary /public/home/yanganyi/Desktop/genexpMR/data/gene_expression/Malhotra_lab/mr_format/SMR/',tissue,'.signif_pairs_5e-8.mateQTL.txt',
                  ' --matrix-eqtl-format --make-besd',
                  ' --out /public/home/yanganyi/Desktop/genexpMR/data/gene_expression/Malhotra_lab/mr_format/SMR/',tissue, sep = '')
  system(command, intern=F)
}



#update .esi, .epi file######################################################
#update .esi file
rm(list =ls())
setwd('F:/类脑/1_实操文件/genexpMR/data/gene_expression/')
library(data.table)

tissues = c('Astrocytes','Endothelial.cells','Excitatory.neurons','Inhibitory.neurons',
            'Microglia','Oligodendrocytes','OPCs...COPs','Pericytes')


g1000_freq = fread('../gwas_summary/g1000_eur.frq', data.table = F)
snp_loc = fread('Malhotra_lab/original_expdata/snp_pos.txt.gz', data.table = F)
snp_loc = snp_loc[,c('SNP','SNP_id_hg19')]
snp_loc = separate(snp_loc, SNP_id_hg19, into= c("chr","BP"), sep= ":")

snp_loc = subset(snp_loc, snp_loc$SNP %in% g1000_freq$SNP)
g1000_freq = subset(g1000_freq, g1000_freq$SNP %in%snp_loc$SNP)

snp_loc = snp_loc[match(g1000_freq$SNP, snp_loc$SNP),]
g1000_freq$BP = snp_loc$BP
head(g1000_freq)

for (i in 1:length(tissues)){
  tissue = tissues[i]
  data = fread(paste('Malhotra_lab/mr_format/SMR/',tissue,'.esi',sep=''), data.table = F)
  head(data)
  g1000_freq_sub = g1000_freq[match(data$V2, g1000_freq$SNP),]
  #1    rs1001  0   744055  A   G   0.23
  #1    rs1002  0   765522  C   G   0.06
  data$V1 = g1000_freq_sub$CHR
  data$V4 = g1000_freq_sub$BP
  data$V5 = g1000_freq_sub$A1
  data$V6 = g1000_freq_sub$A2
  data$V7 = g1000_freq_sub$MAF
  
  outname = paste('Malhotra_lab/mr_format/SMR/',tissue,'.esi',sep='')
  write.table(data, outname, sep='\t', row.names=F, col.names = F, quote=F)
  cat(i, '\n')
}


#update .epi file---------------------------------
rm(list =ls())
setwd('F:/类脑/1_实操文件/genexpMR/data/gene_expression/')
library(data.table)
library(tidyverse)

tissues = c('Astrocytes','Endothelial.cells','Excitatory.neurons','Inhibitory.neurons',
            'Microglia','Oligodendrocytes','OPCs...COPs','Pericytes')


gene_info = fread('F:/类脑/1_实操文件/genexpMR/data/gene_loc/Homo_sapiens.GRCh37.87.gtf.gz', data.table = F, skip = 5)
gene_info = separate(gene_info, col = V9, 
                     into = c("gene_id", "gene_version","gene_name","gene_source","gene_biotype"),sep = "; ")
gene_info = subset(gene_info, gene_info$V3 == 'gene')
gene_info = gene_info[, c('gene_name', 'V1', 'V4', 'V5','V7', 'gene_id','gene_biotype')]
gene_info$gene_id = gsub("gene_id ", "", gene_info$gene_id)
gene_info$gene_id = gsub('"', '', gene_info$gene_id)


for (i in 1:length(tissues)){
  tissue = tissues[i]
  data = fread(paste('Malhotra_lab/mr_format/SMR/',tissue,'.epi',sep=''), data.table = F)
  data =  separate(data, V2, into= c("prefix","gene_id"), sep= "_", remove = F)
  data = data[,-which(names(data)=='prefix')]
    
  #1    probe1001   0   924243  Gene01  +
  #1    probe1002   0   939564  Gene02  -
  gene_info_sub = gene_info[match(data$gene_id, gene_info$gene_id),]
  data$V1 = gene_info_sub$V1
  data$V4 = gene_info_sub$V4
  data$V5 = data$V2
  data$V6 = gene_info_sub$V7
  data = data[,-which(names(data)=='gene_id')] #还用原来的带有前缀的gene_id
  
  outname = paste('Malhotra_lab/mr_format/SMR/',tissue,'.epi',sep='')
  write.table(data, outname, sep='\t', row.names=F, col.names = F, quote=F)
  cat(i, '\n')
}
