rm(list = ls())
setwd('/public/home/yanganyi/Desktop/genexpMR/')
library(data.table)
library(tidyr)


exposure_all = list.files('/public/home/yanganyi/Desktop/Data/gwas_sumstat/zhuhongtu/')
exposure_all = exposure_all[grep('.fastGWA', exposure_all)]


for (exposure in exposure_all){
  file = fread(paste('/public/home/yanganyi/Desktop/Data/gwas_sumstat/zhuhongtu/', exposure, sep=''), data.table = F)
  file$gene = exposure
  file$Zscore = file$BETA / file$SE
  file = subset(file, file$P < 5e-8)
  if (nrow(file) == 0) next
  
  file = file[,c('gene', 'SNP', 'Zscore', 'P')]
  names(file) = c('Gene','SNP','Zscore','pval')
  
  outname = paste('/public/home/yanganyi/Desktop/genexpMR/data/gwas_summary/mr_format/PMR-Egger/IDP/',exposure, '.EUR.signif_pairs_5e-8.txt',sep='')
  write.table(file, outname, sep=' ', row.names=F,quote=F)
  
  selected_snp = unique(file$SNP)
  outname = paste('/public/home/yanganyi/Desktop/genexpMR/data/gwas_summary/mr_format/PMR-Egger/IDP/allele/', exposure,'.EUR.signif_pairs_5e-8.allele',sep='')
  write.table(selected_snp, outname, sep='\t', row.names=F, col.names=F, quote=F)
  
  
  cat(exposure, '\n')
}