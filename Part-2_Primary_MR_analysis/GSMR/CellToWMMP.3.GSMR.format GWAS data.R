rm(list = ls())
setwd("/public/home/yanganyi/Desktop/genexpMR/data/gwas_summary")
library(data.table)
library(tidyr)

path = '/public/home/yanganyi/Desktop/Data/gwas_sumstat/zhuhongtu/'
file_names = list.files('/public/home/yanganyi/Desktop/Data/gwas_sumstat/zhuhongtu/')
file_names = file_names[grep('.fastGWA', file_names)]

g1000_freq = fread('g1000_eur.frq', data.table = F)
names(g1000_freq)[5] = 'freq'


for (i in 1:length(file_names)){
  file = data.table::fread(paste(path, file_names[i], sep = ''), header = T, data.table = F, sep='\t')
  
  file = file[, c('SNP','A1','A2','CHR','BETA','SE','P','N')]
  names(file) = c('SNP','A1','A2','freq','b','se','p','n')
  
  #用g1000的A1,A2和freq
  file = subset(file, file$SNP %in% g1000_freq$SNP)
  g1000_freq_sub = g1000_freq[match(file$SNP, g1000_freq$SNP),]
  file$A1 = g1000_freq_sub$A1
  file$A2 = g1000_freq_sub$A2
  file$freq = g1000_freq_sub$freq
  file = file[, c('SNP','A1','A2','freq','b','se','p','n')]
  
  names(file) = c('SNP','a1','a2','a1_freq','bzy','bzy_se','bzy_pval','bzy_n')
  
  out_name = paste('/public/home/yanganyi/Desktop/genexpMR/data/gwas_summary/mr_format/GSMR/IDP/', file_names[i], '.format.txt', sep = '')
  write.table(file, out_name, row.names = F, col.names= T, sep="\t", quote = F)
}

