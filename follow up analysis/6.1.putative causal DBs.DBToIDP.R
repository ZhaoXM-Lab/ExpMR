rm(list=ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)
library(tidyr)
library(ggplot2)


pheno_info1 = read.table('data/brain_volume_ref.txt', stringsAsFactors = F, header = T)
pheno_info1$file_name = paste('ukb_roi_volume_may12_2019_phase1and2_', pheno_info1$pheno_d,'_allchr_withA2', sep = '')
pheno_info1 = pheno_info1[,c('file_name', 'pheno_name')]
pheno_info1$pheno_name = gsub('[.]', ' ', pheno_info1$pheno_name)
capitalize_first <- function(string) {
  if (nchar(string) > 0) {
    string = paste(toupper(substr(string, 1, 1)), substring(string, 2), sep = "")
  }
  return(string)
}
pheno_info1$pheno_name <- as.character(lapply(pheno_info1$pheno_name, capitalize_first))

pheno_info2 = read.csv('data/110_mean.csv')
pheno_info2 = pheno_info2[grep('FA', pheno_info2$phenocode),]
pheno_info2$file_name = gsub('ukbiobank_mean_','ukb_phase1to3_dti441_dec21_2019_', pheno_info2$assoc_files)
pheno_info2$file_name = gsub('_aug2020.zip','.fastGWA', pheno_info2$file_name)
pheno_info2 = pheno_info2[,c('file_name', 'phenocode')]
names(pheno_info2) = c('file_name', 'pheno_name')

pheno_info = rbind(pheno_info1, pheno_info2)


#PMR-Egger--
file11 = fread('results/PMR-Egger/PMR-Egger.DisorderToIDP.all.txt', data.table = F)
file12 = fread('results/PMR-Egger/PMR-Egger.DisorderToWM.all.txt', data.table = F)
file1 = rbind(file11, file12)
file1 = file1[,c('pheno2','tissue','causal_effect','causal_pvalue','pleiotropy_pvalue')]
names(file1) = c('pheno','tissue','b.pmregger','p.pmregger','p_pleiotropy.pmregger')
file1 = merge(file1, pheno_info, by.x = 'pheno', by.y = 'file_name', all.x = T)
file1 = file1[,c('pheno_name','tissue','b.pmregger','p.pmregger','p_pleiotropy.pmregger')]
names(file1)[1] = 'pheno'

#GSMR--
file21 = fread('results/GSMR/GSMR.DisorderToIDP.all.txt', data.table = F)
file22 = fread('results/GSMR/GSMR.DisorderToWM.all.txt', data.table = F)
file2 = rbind(file21, file22)
file2 = file2[,c('outcome','tissue','bxy','bxy_pval')]
names(file2) = c('pheno','tissue','b.gsmr','p.gsmr')
file2 = merge(file2, pheno_info, by.x = 'pheno', by.y = 'file_name', all.x = T)
file2 = file2[,c('pheno_name','tissue','b.gsmr','p.gsmr')]
names(file2)[1] = 'pheno'

#TwoSampleMR--
file31 = fread('results/TwoSampleMR/TwoSampleMR.DisorderToIDP.all.txt', data.table = F)
file32 = fread('results/TwoSampleMR/TwoSampleMR.DisorderToWM.all.txt', data.table = F) 
file3 = rbind(file31, file32)
file3 = file3[,c('outcome','tissue','method', 'b','pval','Q_pval','pval_egger_intercept')]
names(file3) = c('pheno','tissue','method','b.wr_ivw','p.wr_ivw','Q_pval','pval_egger_intercept')
file3 = merge(file3, pheno_info, by.x = 'pheno', by.y = 'file_name', all.x = T) 
file3 = file3[,c('pheno_name','tissue','method','b.wr_ivw','p.wr_ivw','Q_pval','pval_egger_intercept')]
names(file3)[1] = 'pheno'
file31 = subset(file3, file3$method == "MR Egger")
file32 = subset(file3, file3$method == "Weighted median")
file33 = subset(file3, file3$method == "Inverse variance weighted")
file34 = subset(file3, file3$method == "Simple mode")
file35 = subset(file3, file3$method == "Weighted mode")
file36 = subset(file3, file3$method == "Wald ratio")

names(file31) = c('pheno','tissue','method','b.mr_egger','p.mr_egger','Q_pval.mr_egger','pval_egger_intercept.mr_egger')
names(file32) = c('pheno','tissue','method','b.weighted_median','p.weighted_median','Q_pval.weighted_median','pval_egger_intercept.weighted_median')
names(file33) = c('pheno','tissue','method','b.ivw','p.ivw','Q_pval.ivw','pval_egger_intercept.ivw')
names(file34) = c('pheno','tissue','method','b.simple_mode','p.simple_mode','Q_pval.simple_mode','pval_egger_intercept.simple_mode')
names(file35) = c('pheno','tissue','method','b.weighted_mode','p.weighted_mode','Q_pval.weighted_mode','pval_egger_intercept.weighted_mode')
names(file36) = c('pheno','tissue','method','b.wald_ratio','p.wald_ratio','Q_pval.wald_ratio','pval_egger_intercept.wald_ratio')

file31 = file31[,-which(names(file31)=='method')]
file32 = file32[,-which(names(file32)=='method')]
file33 = file33[,-which(names(file33)=='method')]
file34 = file34[,-which(names(file34)=='method')]
file35 = file35[,-which(names(file35)=='method')]
file36 = file36[,-which(names(file36)=='method')]


#merge
summary = merge(file1, file2, by = c('pheno','tissue'), all = T)
summary = merge(summary, file31, by = c('pheno','tissue'), all = T)
summary = merge(summary, file32, by = c('pheno','tissue'), all = T)
summary = merge(summary, file33, by = c('pheno','tissue'), all = T)
summary = merge(summary, file34, by = c('pheno','tissue'), all = T)
summary = merge(summary, file35, by = c('pheno','tissue'), all = T)
summary = merge(summary, file36, by = c('pheno','tissue'), all = T)

DB_selec = c('CI','EPI','ICH','IS','ALS','MS','OCD','TS',
             'ADHD','AUD','ASD','PTSD','SCZ')
summary = subset(summary, summary$tissue %in% DB_selec)

#significant DBs------------------------
#pmregger--
data1 = summary[,c('pheno','tissue','b.pmregger','p.pmregger')]
data1 = subset(data1, !(is.na(data1$p.pmregger)))
data1$fdr.pmregger = p.adjust(data1$p.pmregger, method = 'fdr', n = nrow(data1))
#gsmr--
data2 = summary[,c('pheno','tissue','b.gsmr','p.gsmr')]
data2 = subset(data2, !(is.na(data2$p.gsmr)))
data2$fdr.gsmr = p.adjust(data2$p.gsmr, method = 'fdr', n = nrow(data2))
#TwoSampleMR--
delete = file3[, c('pheno','tissue','Q_pval','pval_egger_intercept')]
delete = delete[which(delete$Q_pval < 0.05 | delete$pval_egger_intercept < 0.05), ]
delete = delete[!duplicated(paste(delete$pheno, delete$tissue)), ]
delete$pheno_tissue = paste(delete$pheno, delete$tissue, sep = '-')

data31 = summary[,c('pheno','tissue','b.mr_egger','p.mr_egger')]
data31 = subset(data31, !(is.na(data31$p.mr_egger)))
data31$pheno_tissue = paste(data31$pheno, data31$tissue, sep = '-')
data31 = subset(data31, !(data31$pheno_tissue %in% delete$pheno_tissue))
data31 = data31[, -which(names(data31) == 'pheno_tissue')]
data31$fdr.mr_egger = p.adjust(data31$p.mr_egger, method = 'fdr', n = nrow(data31))

data32 = summary[,c('pheno','tissue','b.weighted_median','p.weighted_median')]
data32 = subset(data32, !(is.na(data32$p.weighted_median)))
data32$pheno_tissue = paste(data32$pheno, data32$tissue, sep = '-')
data32 = subset(data32, !(data32$pheno_tissue %in% delete$pheno_tissue))
data32 = data32[, -which(names(data32) == 'pheno_tissue')]
data32$fdr.weighted_median = p.adjust(data32$p.weighted_median, method = 'fdr', n = nrow(data32))

data33 = summary[,c('pheno','tissue','b.ivw','p.ivw')]
data33 = subset(data33, !(is.na(data33$p.ivw)))
data33$pheno_tissue = paste(data33$pheno, data33$tissue, sep = '-')
data33 = subset(data33, !(data33$pheno_tissue %in% delete$pheno_tissue))
data33 = data33[, -which(names(data33) == 'pheno_tissue')]
data33$fdr.ivw = p.adjust(data33$p.ivw, method = 'fdr', n = nrow(data33))

data34 = summary[,c('pheno','tissue','b.simple_mode','p.simple_mode')]
data34 = subset(data34, !(is.na(data34$p.simple_mode)))
data34$pheno_tissue = paste(data34$pheno, data34$tissue, sep = '-')
data34 = subset(data34, !(data34$pheno_tissue %in% delete$pheno_tissue))
data34 = data34[, -which(names(data34) == 'pheno_tissue')]
data34$fdr.simple_mode = p.adjust(data34$p.simple_mode, method = 'fdr', n = nrow(data34))

data35 = summary[,c('pheno','tissue','b.weighted_mode','p.weighted_mode')]
data35 = subset(data35, !(is.na(data35$p.weighted_mode)))
data35$pheno_tissue = paste(data35$pheno, data35$tissue, sep = '-')
data35 = subset(data35, !(data35$pheno_tissue %in% delete$pheno_tissue))
data35 = data35[, -which(names(data35) == 'pheno_tissue')]
data35$fdr.weighted_mode = p.adjust(data35$p.weighted_mode, method = 'fdr', n = nrow(data35))

data36 = summary[,c('pheno','tissue','b.wald_ratio','p.wald_ratio')]
data36 = subset(data36, !(is.na(data36$p.wald_ratio)))
data36$pheno_tissue = paste(data36$pheno, data36$tissue, sep = '-')
data36 = subset(data36, !(data36$pheno_tissue %in% delete$pheno_tissue))
data36 = data36[, -which(names(data36) == 'pheno_tissue')]
data36$fdr.wald_ratio = p.adjust(data36$p.wald_ratio, method = 'fdr', n = nrow(data36))


data1$Ifsigif.pmregger = ifelse(data1$fdr.pmregger < 0.05, 1, 0)
data2$Ifsigif.gsmr = ifelse(data2$fdr.gsmr < 0.05, 1, 0)
data31$Ifsigif.mr_egger = ifelse(data31$fdr.mr_egger < 0.05, 1, 0)
data32$Ifsigif.weighted_median = ifelse(data32$fdr.weighted_median < 0.05, 1, 0)
data33$Ifsigif.ivw = ifelse(data33$fdr.ivw < 0.05, 1, 0)
data34$Ifsigif.simple_mode = ifelse(data34$fdr.simple_mode < 0.05, 1, 0)
data35$Ifsigif.weighted_mode = ifelse(data35$fdr.weighted_mode < 0.05, 1, 0)
data36$Ifsigif.wald_ratio = ifelse(data36$fdr.wald_ratio < 0.05, 1, 0)


summary2 = merge(data1, data2, by = c('pheno','tissue'), all = T)
summary2 = merge(summary2, data31, by = c('pheno','tissue'), all = T)
summary2 = merge(summary2, data32, by = c('pheno','tissue'), all = T)
summary2 = merge(summary2, data33, by = c('pheno','tissue'), all = T)
summary2 = merge(summary2, data34, by = c('pheno','tissue'), all = T)
summary2 = merge(summary2, data35, by = c('pheno','tissue'), all = T)


index = data.frame(summary2$Ifsigif.pmregger, summary2$Ifsigif.gsmr, summary2$Ifsigif.mr_egger,
                   summary2$Ifsigif.weighted_median, summary2$Ifsigif.ivw, summary2$Ifsigif.simple_mode,
                   summary2$Ifsigif.weighted_mode)
summary2$N_Ifsigif = rowSums(index, na.rm = TRUE)
signif_gene = subset(summary2, summary2$N_Ifsigif >= 2)



library(dplyr)
signif_gene$direction = 'Inconsistent'
for (i in 1:nrow(signif_gene)){
  sign_vector = c(sign(signif_gene$b.pmregger[i]),
                  sign(signif_gene$b.gsmr[i]),
                  sign(signif_gene$b.mr_egger[i]),
                  sign(signif_gene$b.weighted_median[i]),
                  sign(signif_gene$b.ivw[i]),
                  sign(signif_gene$b.simple_mode[i]),
                  sign(signif_gene$b.weighted_mode[i])
                  )
  sign_vector = sign_vector[!(is.na(sign_vector))]
  
  if (abs(sum(sign_vector)) == length(sign_vector)){
    signif_gene$direction[i] = 'Consistent'
  }
  cat(length(sign_vector), abs(sum(sign_vector)), '\n')
}
signif_gene = subset(signif_gene, signif_gene$direction == 'Consistent')
signif_gene[is.na(signif_gene)] = '-'
signif_gene$DisorderToIDP = paste(signif_gene$tissue, signif_gene$pheno, sep = '-')

write.table(signif_gene, paste('results/','signif_DisorderToIDP.fdr.new.txt',sep=''), sep="\t", row.names=F, col.names=T, quote=F)
