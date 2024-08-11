rm(list = ls())
setwd('/public/home/yanganyi/Desktop/genexpMR/')
library(data.table)
library(tidyr)


exposure_all = list.files('/public/home/yanganyi/Desktop/Data/gwas_sumstat/yyc_postgwas/ukb/BrainVolume/')
exposure_all = gsub('.txt', '', exposure_all)

for (exposure in exposure_all){
  file = fread(paste('/public/home/yanganyi/Desktop/genexpMR/data/gwas_summary/mr_format/PMR-Egger/IDP/', exposure, '.format.txt', sep=''), data.table = F)
  file$gene = exposure
  file$pval = ifelse(file$Zscore >= 0, 2-2*pnorm(file$Zscore), 2*2*pnorm(file$Zscore))
  
  file = subset(file, file$pval < 5e-8)
  if (nrow(file) == 0) next
  
  file = file[,c('gene', 'SNP', 'Zscore', 'pval')]
  names(file) = c('Gene','SNP','Zscore','pval')
  
  outname = paste('/public/home/yanganyi/Desktop/genexpMR/data/gwas_summary/mr_format/PMR-Egger/IDP/',exposure, '.EUR.signif_pairs_5e-8.txt',sep='')
  write.table(file, outname, sep=' ', row.names=F,quote=F)
  
  selected_snp = unique(file$SNP)
  outname = paste('/public/home/yanganyi/Desktop/genexpMR/data/gwas_summary/mr_format/PMR-Egger/IDP/allele/', exposure,'.EUR.signif_pairs_5e-8.allele',sep='')
  write.table(selected_snp, outname, sep='\t', row.names=F, col.names=F, quote=F)
  
  
  cat(exposure, '\n')
}

