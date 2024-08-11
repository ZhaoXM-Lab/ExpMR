rm(list = ls())
setwd("/public/home/yanganyi/Desktop/genexpMR/data/gwas_summary")
library(data.table)
library(tidyr)

file_names = list.files('/public/home/yanganyi/Desktop/Data/gwas_sumstat/zhuhongtu/')
file_names = file_names[grep('.fastGWA', file_names)]


for (i in 1:length(file_names)){
  file = fread(paste('/public/home/yanganyi/Desktop/Data/gwas_sumstat/zhuhongtu/', file_names[i], sep = ''), header = T, data.table = F, sep='\t')

  file$Zscore = file$BETA / file$SE
  file = file[,c('SNP', 'Zscore')]

  out_name = paste('/public/home/yanganyi/Desktop/genexpMR/data/gwas_summary/mr_format/PMR-Egger/IDP/', 
                   file_names[i], '.format.txt', sep = '')
  write.table(file, out_name, row.names = F, col.names= T, sep="\t", quote = F)
  cat(i, '\t', file_names[i])
}
