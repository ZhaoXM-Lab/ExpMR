rm(list = ls())
setwd("/public/home/yanganyi/Desktop/genexpMR/data/gwas_summary")
library(data.table)
library(tidyr)

path = '/public/home/yanganyi/Desktop/Data/gwas_sumstat/yyc_postgwas/ukb/BrainVolume/'
file_names = list.files(path)
file_names = file_names[grep('ukb', file_names)]

g1000_freq = fread('g1000_eur.frq', data.table = F)
names(g1000_freq)[5] = 'freq'


format_gwas = function(i){
  file = data.table::fread(paste(path, file_names[i], sep = ''), header = F, data.table = F, sep=' ')
  file = tidyr::separate(file, V9, into= c("V10","V11"), sep= "\t")
  names(file) = c('CHR','SNP','BP','A1','TEST', 'NMISS','BETA','STAT','P','A2')
  
  file$se = file$BETA / file$STAT
  file$n = 19629
  
  file = file[,c('SNP','A1','A2','TEST','BETA','se','P','n')]
  names(file) = c('SNP','A1','A2','freq','b','se','p','n')
  
  file = subset(file, file$SNP %in% g1000_freq$SNP)
  g1000_freq_sub = g1000_freq[match(file$SNP, g1000_freq$SNP),]
  file$A1 = g1000_freq_sub$A1
  file$A2 = g1000_freq_sub$A2
  file$freq = g1000_freq_sub$freq
  file = file[, c('SNP','A1','A2','freq','b','se','p','n')]
  
  names(file) = c('SNP','a1','a2','a1_freq','bzy','bzy_se','bzy_pval','bzy_n')
  
  out_name = paste('/public/home/yanganyi/Desktop/genexpMR/data/gwas_summary/mr_format/GSMR/IDP/', 
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


