rm(list = ls())
setwd("/public/home/yanganyi/Desktop/genexpMR/data/gwas_summary")
library(data.table)
library(tidyr)

file_names = list.files('/public/home/yanganyi/Desktop/Data/gwas_sumstat/yyc_postgwas/ukb/BrainVolume/')
file_names = file_names[grep('ukb', file_names)]

format_gwas = function(i){
  file = data.table::fread(paste('/public/home/yanganyi/Desktop/Data/gwas_sumstat/yyc_postgwas/ukb/BrainVolume/', file_names[i], sep = ''), header = F, data.table = F, sep=' ')
  file = tidyr::separate(file, V9, into= c("V10","V11"), sep= "\t")
  names(file) = c('CHR','SNP','BP','A1','TEST', 'NMISS','BETA','STAT','P','A2')

  file$Zscore = file$STAT
  file = file[,c('SNP', 'Zscore')]
  out_name = paste('/public/home/yanganyi/Desktop/genexpMR/data/gwas_summary/mr_format/PMR-Egger/IDP/', 
                   gsub('.txt', '', file_names[i]), '.format.txt', sep = '')
  write.table(file, out_name, row.names = F, col.names= T, sep="\t", quote = F)
  return(0)
}

library(doParallel)
library(foreach)

cl = makeCluster(5) 
registerDoParallel(cl)
result = foreach(i = 1:length(file_names), .combine = 'rbind') %dopar% format_gwas(i)

stopCluster(cl)
