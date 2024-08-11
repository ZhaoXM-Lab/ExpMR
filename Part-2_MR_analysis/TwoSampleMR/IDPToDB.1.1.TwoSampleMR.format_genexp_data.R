rm(list =ls())
setwd('/public/home/yanganyi/Desktop/genexpMR/')
library(data.table)
library(tidyr)


exposure_all = list.files('/public/home/yanganyi/Desktop/Data/gwas_sumstat/yyc_postgwas/ukb/BrainVolume/')
exposure_all = gsub('.txt', '', exposure_all)


for (exposure in exposure_all){
  file = fread(paste('/public/home/yanganyi/Desktop/genexpMR/data/gwas_summary/mr_format/SMR/IDP/', exposure, '.format.txt', sep=''), data.table = F)
  file = subset(file, file$p < 5e-8)
  if (nrow(file) == 0) next
  file$gene = exposure
  file = file[, c('SNP','gene', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'n')]
  names(file) = c('SNP','gene','effect_allele','other_allele','eaf','beta','se','pval', 'N')
  
  outname = paste('/public/home/yanganyi/Desktop/genexpMR/data/gwas_summary/mr_format/SMR/IDP/',exposure, '.EUR.signif_pairs_5e-8.txt',sep='')
  write.table(file, outname, sep=' ', row.names=F,quote=F)
  
  cat(exposure, '\n')
}