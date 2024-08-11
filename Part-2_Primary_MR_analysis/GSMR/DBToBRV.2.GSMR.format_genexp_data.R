rm(list = ls())
setwd('/public/home/yanganyi/Desktop/genexpMR/')
library(data.table)
library(tidyr)


trait_info = read.csv('data/Nsample.csv', header = T, stringsAsFactors = F)
exposure_all = trait_info[which(trait_info$Category == 'GWAS'), 'Trait']

for (exposure in exposure_all){
  file = fread(paste('/public/home/yanganyi/Desktop/genexpMR/data/gwas_summary/mr_format/GSMR/brain_disorder/', exposure, '.format.txt', sep=''), data.table = F)
  file = subset(file, file$bzy_pval < 5e-8)
  
  if (nrow(file) == 0) next
  file$gene = exposure
  file = file[, c('SNP','gene', 'a1', 'a2', 'a1_freq', 'bzy', 'bzy_se', 'bzy_pval', 'bzy_n')]
  names(file) = c('SNP','gene','a1','a2','a1_freq','bzx','bzx_se','bzx_pval', 'bzx_n')
  
  outname = paste('/public/home/yanganyi/Desktop/genexpMR/data/gwas_summary/mr_format/GSMR/brain_disorder/',exposure, '.EUR.signif_pairs_5e-8.txt',sep='')
  write.table(file, outname, sep=' ', row.names=F,quote=F)
  
  cat(exposure, '\n')
}
