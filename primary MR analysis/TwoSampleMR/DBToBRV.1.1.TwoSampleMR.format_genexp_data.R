rm(list = ls())
setwd("F:/类脑/1_实操文件/genexpMR/data/gwas_summary/pheno/brain_disorder")
library(data.table)
library(tidyr)


trait_info = read.csv('../../../Nsample.csv', header = T, stringsAsFactors = F)
exposure_all = trait_info[which(trait_info$Category == 'GWAS'), 'Trait']

for (exposure in exposure_all){
  file = fread(paste('../../mr_format/SMR/', exposure, '.format.txt', sep=''), data.table = F)
  file = subset(file, file$p < 5e-8)
  if (nrow(file) == 0) next
  file$gene = exposure
  file = file[, c('SNP','gene', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'n')]
  names(file) = c('SNP','gene','effect_allele','other_allele','eaf','beta','se','pval', 'N')
  
  outname = paste('../../mr_format/SMR/',exposure, '.EUR.signif_pairs_5e-8.txt',sep='')
  write.table(file, outname, sep=' ', row.names=F,quote=F)
  
  cat(exposure, '\n')
}