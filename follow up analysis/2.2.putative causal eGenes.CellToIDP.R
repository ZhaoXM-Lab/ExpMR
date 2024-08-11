rm(list=ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)
library(tidyr)
library(ggplot2)
library(tidyverse)


#name mapping
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

#TWAS--
file11 = fread('results/TWAS/TWAS.14BrainTissuesToIDP.allmosetsignif.txt', data.table = F) 
file12 = fread('results/TWAS/TWAS.14BrainTissuesToWM.allmosetsignif.txt', data.table = F) 
file1 = rbind(file11, file12)
file1 = file1[,c('pheno','tissue','gene','b.twas','p.twas')]
file1 = merge(file1, pheno_info, by.x = 'pheno', by.y = 'file_name', all.x = T)
file1 = file1[,c('pheno_name','tissue','gene','b.twas','p.twas')]; names(file1)[1] = 'pheno'
file1 = file1[,-which(names(file1) == 'tissue')]
#PMR-Egger--
file21 = fread('results/PMR-Egger/PMR-Egger.CellToIDP.all.txt', data.table = F)  #volume
file22 = fread('results/PMR-Egger/PMR-Egger.CellToWM.all.txt', data.table = F)  #white matter
file2 = rbind(file21, file22)
file2 = file2[,c('pheno2','tissue','pheno1','causal_effect','causal_pvalue','pleiotropy_pvalue')]
names(file2) = c('pheno','tissue','gene','b.pmregger','p.pmregger','p_pleiotropy.pmregger')
file2 = separate(file2, 'gene', into= c("gene_name","gene"),sep= "_") 
file2 = file2[,-which(names(file2) == 'gene_name')]
file2 = merge(file2, pheno_info, by.x = 'pheno', by.y = 'file_name', all.x = T)
file2 = file2[,c('pheno_name','tissue','gene','b.pmregger','p.pmregger','p_pleiotropy.pmregger')]
names(file2)[1] = 'pheno'
#SMR--
file31 = fread('results/SMR/SMR.CellToIDP.all.txt', data.table = F)
file32 = fread('results/SMR/SMR.CellToWM.all.txt', data.table = F)
file3 = rbind(file31, file32)
file3 = file3[,c('pheno','tissue','Gene','b_SMR','p_SMR','p_HEIDI')]
names(file3) = c('pheno','tissue','gene','b.smr','p.smr','p_HEIDI.smr')
file3 = separate(file3, 'gene', into= c("gene_name","gene"),sep= "_") 
file3 = file3[,-which(names(file3) == 'gene_name')]
file3 = merge(file3, pheno_info, by.x = 'pheno', by.y = 'file_name', all.x = T)
file3 = file3[,c('pheno_name','tissue','gene','b.smr','p.smr','p_HEIDI.smr')]
names(file3)[1] = 'pheno'
#TwoSampleMR--
file41 = fread('results/TwoSampleMR/TwoSampleMR.CellToIDP.all.txt', data.table = F)
file42 = fread('results/TwoSampleMR/TwoSampleMR.CellToWM.all.txt', data.table = F)
file4 = rbind(file41, file42)
file4 = subset(file4, file4$method %in% c('Wald ratio','Inverse variance weighted'))
file4 = file4[,c('outcome','tissue','exposure','b','pval','Q_pval','pval_egger_intercept')]
names(file4) = c('pheno','tissue','gene','b.wr_ivw','p.wr_ivw', 'Q_pval.wr_ivw', 'pval_egger_intercept.wr_ivw')
file4 = separate(file4, 'gene', into= c("gene_name","gene"),sep= "_") 
file4 = file4[,-which(names(file4) == 'gene_name')]
file4 = merge(file4, pheno_info, by.x = 'pheno', by.y = 'file_name', all.x = T)
file4 = file4[,c('pheno_name','tissue','gene','b.wr_ivw','p.wr_ivw', 'Q_pval.wr_ivw', 'pval_egger_intercept.wr_ivw')]
names(file4)[1] = 'pheno'
#GSMR--
file51 = fread('results/GSMR/GSMR.CellToIDP.all.txt', data.table = F)
file52 = fread('results/GSMR/GSMR.CellToWM.all.txt', data.table = F)
file5 = rbind(file51, file52)
file5 = file5[,c('outcome','tissue','exposure','bxy','bxy_pval')]
names(file5) = c('pheno','tissue','gene','b.gsmr','p.gsmr')
file5 = separate(file5, 'gene', into= c("gene_name","gene"),sep= "_")
file5 = file5[,-which(names(file5) == 'gene_name')]
file5 = merge(file5, pheno_info, by.x = 'pheno', by.y = 'file_name', all.x = T) 
file5 = file5[,c('pheno_name','tissue','gene','b.gsmr','p.gsmr')]
names(file5)[1] = 'pheno'


summary = merge(file1, file2, by = c('pheno','gene'), all = T)
summary = merge(summary, file3, by = c('pheno','tissue','gene'), all = T)
summary = merge(summary, file4, by = c('pheno','tissue','gene'), all = T)
summary = merge(summary, file5, by = c('pheno','tissue','gene'), all = T)
 

nrow(file1[!is.na(file1$p.twas), ]) 
nrow(file3[!is.na(file3$p.smr), ])
nrow(file4[!is.na(file4$p.wr_ivw), ])
nrow(file2[!is.na(file2$p.pmregger), ])
nrow(file5[!is.na(file5$p.gsmr), ]) 



summary1 = subset(summary, summary$p.twas < 0.05); nrow(summary1)
length(unique(paste(summary1$pheno, summary1$gene)))
nrow(summary1)
nrow(summary1[!is.na(summary1$p.smr), ])
nrow(summary1[!is.na(summary1$p.wr_ivw), ])
nrow(summary1[!is.na(summary1$p.pmregger), ])
nrow(summary1[!is.na(summary1$p.gsmr), ])



calcu_Ngene = function(data){
  Ngene = data.frame(stringsAsFactors = F)
  for (pheno_idx in pheno_all){
    subset_file = subset(data, data$pheno == pheno_idx)
    N = nrow(subset_file)
    Ngene.tmp = data.frame(pheno = pheno_idx, N = N)
    Ngene = rbind(Ngene.tmp, Ngene)
  }
  return(Ngene)
}
pheno_all = Reduce(intersect, list(file1$pheno,file2$pheno,file3$pheno,file4$pheno,file5$pheno))
Ngene = calcu_Ngene(summary1)

Ngene = Ngene[order(Ngene$N, decreasing = T), ]
Ngene$pheno = factor(Ngene$pheno, levels = unique(Ngene$pheno))

ggplot(Ngene, aes(x = pheno, y = N)) +
  geom_bar(stat = "identity", colour = "black", fill = "#FED9A6") +
  theme(axis.text.x = element_text(size=15, angle=90, hjust=1, vjust=1), 
        axis.text.y = element_text(size=8), 
        legend.title = element_text(), legend.position =  c(0.85, 0.1)) +
  ylab('Number of signature genes') +
  xlab('')



pheno_all = Reduce(intersect, list(file1$pheno,file2$pheno,file3$pheno,file4$pheno,file5$pheno))
tissue_all = Reduce(intersect, list(file2$tissue,file3$tissue,file4$tissue,file5$tissue))

calcu_Ngene = function(data){
  Ngene = data.frame(stringsAsFactors = F)
  for (pheno_idx in pheno_all){
    for (tissue_idx in tissue_all){
      subset_file = subset(data, data$tissue == tissue_idx & data$pheno == pheno_idx)
      N = nrow(subset_file)
      Ngene.tmp = data.frame(pheno = pheno_idx, tissue = tissue_idx, N = N)
      Ngene = rbind(Ngene.tmp, Ngene)
    }
  }
  return(Ngene)
}

#smr
data1 = summary1[,c('pheno','tissue','gene', 'p.smr', 'p_HEIDI.smr')]
data1 = subset(data1, !(is.na(data1$p.smr)) & data1$p_HEIDI.smr >0.05)
data1$fdr.smr = p.adjust(data1$p.smr, method = 'fdr', n = nrow(data1))
data1 = subset(data1, data1$fdr.smr < 0.05)
Ngene1 = calcu_Ngene(data1); Ngene1$method = 'smr'
#wr_ivw
data2 = summary1[,c('pheno','tissue','gene', 'p.wr_ivw', 'Q_pval.wr_ivw', 'pval_egger_intercept.wr_ivw')]
data2 = subset(data2, !(is.na(data2$p.wr_ivw)))
data2 = subset(data2, (data2$Q_pval.wr_ivw>=0.05) | is.na(data2$Q_pval.wr_ivw))
data2 = subset(data2, (data2$pval_egger_intercept.wr_ivw>=0.05) | is.na(data2$pval_egger_intercept.wr_ivw))
data2$fdr.wr_ivw = p.adjust(data2$p.wr_ivw, method = 'fdr', n = nrow(data2))
data2 = subset(data2, data2$fdr.wr_ivw < 0.05)
Ngene2 = calcu_Ngene(data2); Ngene2$method = 'wr_ivw'
#pmregger
data3 = summary1[,c('pheno','tissue','gene', 'p.pmregger')]
data3 = subset(data3, !(is.na(data3$p.pmregger)))
data3$fdr.pmregger = p.adjust(data3$p.pmregger, method = 'fdr', n = nrow(data3))
data3 = subset(data3, data3$fdr.pmregger < 0.05)
Ngene3 = calcu_Ngene(data3); Ngene3$method = 'pmregger'
#gsmr
data4 = summary1[,c('pheno','tissue','gene', 'p.gsmr')]
data4 = subset(data4, !(is.na(data4$p.gsmr)))
data4$fdr.gsmr = p.adjust(data4$p.gsmr, method = 'fdr', n = nrow(data4))
data4 = subset(data4, data4$fdr.gsmr < 0.05)
Ngene4 = calcu_Ngene(data4); Ngene4$method = 'gsmr'

names(Ngene1) = names(Ngene2) = names(Ngene3) = names(Ngene4)
Ngene = Reduce(rbind, list(Ngene1, Ngene2, Ngene3, Ngene4))

Ngene = Ngene[order(Ngene$N, decreasing = T), ]
Ngene$pheno = factor(Ngene$pheno, levels = unique(Ngene$pheno))
Ngene$method = factor(Ngene$method, levels = unique(Ngene$method))
ggplot(Ngene) +
  geom_boxplot(aes(pheno, N, fill = method), varwidth = F,colour = 'black', outlier.colour = 'red') +
  theme(axis.text.x = element_text(size = 13, angle=90, hjust=1, vjust=1), 
        axis.text.y = element_text(size = 13)) +
  ylab('Number of signif genes') +
  xlab('')



#signif causal eGenes------------------------
#smr
data1 = summary1[,c('pheno','tissue','gene','b.smr', 'p.smr', 'p_HEIDI.smr')]
data1 = subset(data1, !(is.na(data1$p.smr)) & data1$p_HEIDI.smr >0.05)
data1$fdr.smr = p.adjust(data1$p.smr, method = 'fdr', n = nrow(data1))
#wr_ivw
data2 = summary1[,c('pheno','tissue','gene','b.wr_ivw','p.wr_ivw', 'Q_pval.wr_ivw', 'pval_egger_intercept.wr_ivw')]
data2 = subset(data2, !(is.na(data2$p.wr_ivw)))
data2 = subset(data2, (data2$Q_pval.wr_ivw>=0.05) | is.na(data2$Q_pval.wr_ivw))
data2 = subset(data2, (data2$pval_egger_intercept.wr_ivw>=0.05) | is.na(data2$pval_egger_intercept.wr_ivw)) 
data2$fdr.wr_ivw = p.adjust(data2$p.wr_ivw, method = 'fdr', n = nrow(data2))
#pmregger
data3 = summary1[,c('pheno','tissue','gene','b.pmregger','p.pmregger')]
data3 = subset(data3, !(is.na(data3$p.pmregger)))
data3$fdr.pmregger = p.adjust(data3$p.pmregger, method = 'fdr', n = nrow(data3))
#gsmr
data4 = summary1[,c('pheno','tissue','gene','b.gsmr','p.gsmr')]
data4 = subset(data4, !(is.na(data4$p.gsmr)))
data4$fdr.gsmr = p.adjust(data4$p.gsmr, method = 'fdr', n = nrow(data4))
data1$Ifsigif.smr = ifelse(data1$fdr.smr < 0.05, T, F)
data2$Ifsigif.wr_ivw = ifelse(data2$fdr.wr_ivw < 0.05, T, F)
data3$Ifsigif.pmregger = ifelse(data3$fdr.pmregger < 0.05, T, F)
data4$Ifsigif.gsmr = ifelse(data4$fdr.gsmr < 0.05, T, F)

summary2 = merge(data1, data2, by = c('pheno','tissue','gene'), all = T)
summary2 = merge(summary2, data3, by = c('pheno','tissue','gene'), all = T)
summary2 = merge(summary2, data4, by = c('pheno','tissue','gene'), all = T)


signif_gene = subset(summary2, (summary2$fdr.smr < 0.05 | summary2$fdr.wr_ivw < 0.05) & 
                                    (summary2$fdr.pmregger < 0.05 | summary2$fdr.gsmr < 0.05))

library(dplyr)
signif_gene$direction = 'Inconsistent'
for (i in 1:nrow(signif_gene)){
  sign_vector = c(sign(signif_gene$b.smr[i]),sign(signif_gene$b.wr_ivw[i]),sign(signif_gene$b.pmregger[i]),sign(signif_gene$b.gsmr[i]))
  sign_vector = sign_vector[!(is.na(sign_vector))]
  
  if (abs(sum(sign_vector)) == length(sign_vector)){
    signif_gene$direction[i] = 'Consistent'
  }
  cat(abs(sum(sign_vector)), '\n')
}
signif_gene = subset(signif_gene, signif_gene$direction == 'Consistent')
signif_gene[is.na(signif_gene)] = '-'


gene_info = fread('F:/类脑/1_实操文件/genexpMR/data/gene_loc/Homo_sapiens.GRCh37.87.gtf.gz', data.table = F, skip = 5)
gene_info = separate(gene_info, col = V9, 
                     into = c("gene_id", "gene_version","gene_name","gene_source","gene_biotype"),sep = "; ")
gene_info = subset(gene_info, gene_info$V3 == 'gene')
gene_info = gene_info[, c('gene_name', 'V1', 'V4', 'V5','V7', 'gene_id','gene_biotype')]
gene_info$gene_id = gsub("gene_id ", "", gene_info$gene_id)
gene_info$gene_id = gsub('"', '', gene_info$gene_id)
gene_info$gene_name = gsub("gene_name ", "", gene_info$gene_name)
gene_info$gene_name = gsub('"', '', gene_info$gene_name)
gene_info$gene_biotype = gsub("gene_biotype ", "", gene_info$gene_biotype)
gene_info$gene_biotype = gsub('"', '', gene_info$gene_biotype)
gene_info$gene_biotype = gsub(';','', gene_info$gene_biotype)

signif_gene = merge(signif_gene, gene_info[,c('gene_name', 'gene_id','gene_biotype')], by.x = 'gene', by.y = 'gene_id', all.x = T)
signif_gene = signif_gene[,c(2,3,1,ncol(signif_gene)-1,ncol(signif_gene),4:(ncol(signif_gene)-2))]
length(unique(signif_gene$gene))
length(unique(signif_gene$pheno))
length(unique(signif_gene$tissue))


write.table(signif_gene, paste('results/','signif_genes.CellToIDP.fdr.txt',sep=''), sep="\t", row.names=F, col.names=T, quote=F)
