#CellToDB##############################################
rm(list=ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)


signif_gene = fread('results/signif_genes.CellToDisorder.fdr.txt',data.table = F)
signif_gene$gene_pheno = paste(signif_gene$gene, signif_gene$pheno, sep = '-')
find_gene_pheno = unique(signif_gene$gene_pheno)


#SMR
rfile1 = fread('results/SMR/SMR.BrainMetaToDisorder.all.txt', data.table = F)
rfile1 = rfile1[,c('pheno','tissue','Gene','b_SMR','p_SMR','p_HEIDI')]
names(rfile1) = c('pheno','tissue','gene','b.smr','p.smr','p_HEIDI.smr')
rfile1 = tidyr::separate(rfile1, gene, c('gene_id', 'gene_version'), sep = '[.]')
rfile1$gene_pheno = paste(rfile1$gene_id, rfile1$pheno, sep = '-')
#wr_ivw
rfile2 = fread('results/TwoSampleMR/TwoSampleMR.BrainMetaToDisorder.all.txt', data.table = F)
rfile2 = subset(rfile2, rfile2$method %in% c('Wald ratio','Inverse variance weighted'))
rfile2 = rfile2[,c('outcome','tissue','exposure','b','pval','Q_pval','pval_egger_intercept')]
names(rfile2) = c('pheno','tissue','gene','b.wr_ivw','p.wr_ivw', 'Q_pval.wr_ivw', 'pval_egger_intercept.wr_ivw')
rfile2 = tidyr::separate(rfile2, gene, c('gene_id', 'gene_version'), sep = '[.]')
rfile2$gene_pheno = paste(rfile2$gene_id, rfile2$pheno, sep = '-')
#PMR-Egger
rfile3 = fread('results/PMR-Egger/PMR-Egger.BrainMetaToDisorder.all.txt', data.table = F)
rfile3 = rfile3[,c('pheno2','tissue','pheno1','causal_effect','causal_pvalue','pleiotropy_pvalue')]
names(rfile3) = c('pheno','tissue','gene','b.pmregger','p.pmregger','p_pleiotropy.pmregger')
rfile3= tidyr::separate(rfile3, gene, c('gene_id', 'gene_version'), sep = '[.]')
rfile3$gene_pheno = paste(rfile3$gene_id, rfile3$pheno, sep = '-')
#GSMR
rfile4 = fread('results/GSMR/GSMR.BrainMetaToDisorder.all.txt', data.table = F)
rfile4 = rfile4[,c('outcome','tissue','exposure','bxy','bxy_pval')]
names(rfile4) = c('pheno','tissue','gene','b.gsmr','p.gsmr')
rfile4 = tidyr::separate(rfile4, col=gene, into=c('gene_id', 'gene_version'), sep = '[.]')
rfile4$gene_pheno = paste(rfile4$gene_id, rfile4$pheno, sep = '-')


#replication rate---------------------
cutoff = 0.05
rdata1 = subset(rfile1, rfile1$p.smr < cutoff& rfile1$p_HEIDI.smr > 0.05)
rdata2 = subset(rfile2, rfile2$p.wr_ivw < cutoff)
rdata2 = subset(rdata2, (rdata2$Q_pval.wr_ivw>=0.05) | is.na(rdata2$Q_pval.wr_ivw)) 
rdata2 = subset(rdata2, (rdata2$pval_egger_intercept.wr_ivw>=0.05) | is.na(rdata2$pval_egger_intercept.wr_ivw))
rdata3 = subset(rfile3, rfile3$p.pmregger < cutoff)
rdata4 = subset(rfile4, rfile4$p.gsmr < cutoff)

rep_gene_pheno = Reduce(union, list(rdata1$gene_pheno, rdata2$gene_pheno, rdata3$gene_pheno, rdata4$gene_pheno))
rep_gene_pheno = unique(rep_gene_pheno)

replicated_pairs = intersect(rep_gene_pheno, find_gene_pheno)
length(replicated_pairs) / length(find_gene_pheno) 

signif_gene$If_repl = ifelse(signif_gene$gene_pheno %in% replicated_pairs, 'Y','N')
signif_gene = signif_gene[,-which(names(signif_gene) == 'gene_pheno')]
write.table(signif_gene, 'results/signif_genes.repBrainMeta.CellToDisorder.fdr.txt', sep="\t", row.names=F, col.names=T, quote=F)



length(rep_gene_pheno); length(replicated_pairs); length(find_gene_pheno)
q = length(replicated_pairs)
m = length(find_gene_pheno)
n = 26861*length(unique(signif_gene$pheno)) - m
k = length(rep_gene_pheno)
phyper(q-1, m, n, k, lower.tail=F)
m-q; k-q



#CellToIDP##############################################
rm(list=ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)


signif_gene = fread('results/signif_genes.CellToIDP.fdr.txt',data.table = F)
signif_gene$gene_pheno = paste(signif_gene$gene, signif_gene$pheno, sep = '-')
find_gene_pheno = unique(signif_gene$gene_pheno)


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


#SMR
rfile1 = fread('results/SMR/SMR.BrainMetaToIDP.all.txt', data.table = F)
rfile1 = rfile1[,c('pheno','tissue','Gene','b_SMR','p_SMR','p_HEIDI')]
names(rfile1) = c('pheno','tissue','gene','b.smr','p.smr','p_HEIDI.smr')
rfile1 = tidyr::separate(rfile1, gene, c('gene_id', 'gene_version'), sep = '[.]')
rfile1 = merge(rfile1, pheno_info, by.x = 'pheno', by.y = 'file_name', all.x = T)
rfile1$gene_pheno = paste(rfile1$gene_id, rfile1$pheno_name, sep = '-')
#wr_ivw
rfile2 = fread('results/TwoSampleMR/TwoSampleMR.BrainMetaToIDP.all.txt', data.table = F)
rfile2 = subset(rfile2, rfile2$method %in% c('Wald ratio','Inverse variance weighted'))
rfile2 = rfile2[,c('outcome','tissue','exposure','b','pval','Q_pval','pval_egger_intercept')]
names(rfile2) = c('pheno','tissue','gene','b.wr_ivw','p.wr_ivw', 'Q_pval.wr_ivw', 'pval_egger_intercept.wr_ivw')
rfile2 = tidyr::separate(rfile2, gene, c('gene_id', 'gene_version'), sep = '[.]')
rfile2 = merge(rfile2, pheno_info, by.x = 'pheno', by.y = 'file_name', all.x = T)
rfile2$gene_pheno = paste(rfile2$gene_id, rfile2$pheno_name, sep = '-')
#PMR-Egger
rfile3 = fread('results/PMR-Egger/PMR-Egger.BrainMetaToIDP.all.txt', data.table = F)
rfile3 = rfile3[,c('pheno2','tissue','pheno1','causal_effect','causal_pvalue','pleiotropy_pvalue')]
names(rfile3) = c('pheno','tissue','gene','b.pmregger','p.pmregger','p_pleiotropy.pmregger')
rfile3 = tidyr::separate(rfile3, gene, c('gene_id', 'gene_version'), sep = '[.]')
rfile3 = merge(rfile3, pheno_info, by.x = 'pheno', by.y = 'file_name', all.x = T)
rfile3$gene_pheno = paste(rfile3$gene_id, rfile3$pheno_name, sep = '-')
#GSMR
rfile4 = fread('results/GSMR/GSMR.BrainMetaToIDP.all.txt', data.table = F)
rfile4 = rfile4[,c('outcome','tissue','exposure','bxy','bxy_pval')]
names(rfile4) = c('pheno','tissue','gene','b.gsmr','p.gsmr')
rfile4 = tidyr::separate(rfile4, col=gene, into=c('gene_id', 'gene_version'), sep = '[.]')
rfile4 = merge(rfile4, pheno_info, by.x = 'pheno', by.y = 'file_name', all.x = T)
rfile4$gene_pheno = paste(rfile4$gene_id, rfile4$pheno_name, sep = '-')


#replication rate---------------------
cutoff = 0.05
rdata1 = subset(rfile1, rfile1$p.smr < cutoff& rfile1$p_HEIDI.smr > 0.05)
rdata2 = subset(rfile2, rfile2$p.wr_ivw < cutoff)
rdata2 = subset(rdata2, (rdata2$Q_pval.wr_ivw>=0.05) | is.na(rdata2$Q_pval.wr_ivw))
rdata2 = subset(rdata2, (rdata2$pval_egger_intercept.wr_ivw>=0.05) | is.na(rdata2$pval_egger_intercept.wr_ivw))
rdata3 = subset(rfile3, rfile3$p.pmregger < cutoff)
rdata4 = subset(rfile4, rfile4$p.gsmr < cutoff)
#Replication rate = number of replication associations / number of both replicated and not replicated associations
rep_gene_pheno = Reduce(union, list(rdata1$gene_pheno, rdata2$gene_pheno, rdata3$gene_pheno, rdata4$gene_pheno))
rep_gene_pheno = unique(rep_gene_pheno)

replicated_pairs = intersect(rep_gene_pheno, find_gene_pheno)
length(replicated_pairs) / length(find_gene_pheno) 

signif_gene$If_repl = ifelse(signif_gene$gene_pheno %in% replicated_pairs, 'Y','N')
signif_gene = signif_gene[,-which(names(signif_gene) == 'gene_pheno')]
write.table(signif_gene, 'results/signif_genes.repBrainMeta.CellToIDP.fdr.txt', sep="\t", row.names=F, col.names=T, quote=F)


length(rep_gene_pheno); length(replicated_pairs); length(find_gene_pheno)
q = length(replicated_pairs)
m = length(find_gene_pheno)
n = 26861*length(unique(signif_gene$pheno)) - m
k = length(rep_gene_pheno)
phyper(q-1, m, n, k, lower.tail=F)
m-q; k-q

