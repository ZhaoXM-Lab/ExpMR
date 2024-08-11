#Merge TWAS results+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rm(list =ls())
setwd('F:/类脑/1_实操文件/genexpMR')
library(data.table)

trait_info = read.csv('data/Nsample.csv', header = T, stringsAsFactors = F)
tissue_all = trait_info[which(trait_info$Category == 'GTEx'), 'Trait']
tissue_all = c(tissue_all, 'DLPFC')
pheno_all0 = read.csv('data/110_mean.csv')
pheno_all0 = pheno_all0[grep('FA', pheno_all0$phenocode),]
pheno_all = pheno_all0$assoc_files
pheno_all = gsub('ukbiobank_mean_','ukb_phase1to3_dti441_dec21_2019_', pheno_all)
pheno_all = gsub('_aug2020.zip','.fastGWA', pheno_all)

merge.data = data.frame(stringsAsFactors = F)
for (pheno in pheno_all){
  for (tissue in tissue_all){
    file_name = paste('results/TWAS/22WM_14GTEx/', pheno,'.',tissue, '.MetaXcan.csv', sep='')
    
    new.data = read.table(file_name, header = T, sep = ",", stringsAsFactors = F)
    new.data$pheno = pheno
    new.data$tissue = tissue
    new.data = new.data[,c(ncol(new.data)-1,ncol(new.data),1:(ncol(new.data)-2))]
    
    merge.data = rbind(merge.data, new.data)
  }
  cat(pheno, 'is OK.\n')
}
write.table(merge.data, 'results/TWAS/TWAS.14BrainTissuesToWM.all.txt', sep="\t", row.names=F, col.names=T, quote=F)



#Pre-select gene-phenotype pairs+++++++++++++++++++++++++++++++++++++++++++++++++++++++
rm(list=ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)
library(tidyr)

file = fread('results/TWAS/TWAS.14BrainTissuesToWM.all.txt', data.table = F)
file = file[,c('pheno','tissue','gene','effect_size','pvalue')]
names(file) = c('pheno','tissue','gene','b.twas','p.twas')

file = tidyr::separate(file, gene, into= c("gene_id","version"),sep= "[.]")
file = file[,-which(names(file) == 'version')]
names(file)[which(names(file)=='gene_id')] = 'gene'

trait_info = read.csv('data/Nsample.csv', header = T, stringsAsFactors = F)
tissue_all = trait_info[which(trait_info$Category == 'GTEx'), 'Trait']
tissue_all = c(tissue_all, 'DLPFC')
pheno_all0 = read.csv('data/110_mean.csv')
pheno_all0 = pheno_all0[grep('FA', pheno_all0$phenocode),]
pheno_all = pheno_all0$assoc_files
pheno_all = gsub('ukbiobank_mean_','ukb_phase1to3_dti441_dec21_2019_', pheno_all)
pheno_all = gsub('_aug2020.zip','.fastGWA', pheno_all)


run_each = function(pheno){
  sub_file = file[(file$pheno == pheno), ]
  result = data.frame(stringsAsFactors = F)
  for (gene in unique(sub_file$gene)){
    sub_sub_file = sub_file[sub_file$gene == gene, ]
    result.tmp = sub_sub_file[which.min(sub_sub_file$p.twas), ]
    result = rbind(result, result.tmp)
  }
  return(result)
}

library(doParallel)
library(foreach)

cl = makeCluster(8) 
registerDoParallel(cl)
result = foreach(pheno = pheno_all, .combine = 'rbind') %dopar% run_each(pheno)
stopCluster(cl)
write.table(result, 'results/TWAS/TWAS.14BrainTissuesToWM.allmosetsignif.txt', sep="\t", row.names=F, col.names=T, quote=F)




#Merge PMR-Egger results++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rm(list =ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)

tissue_all = c('Astrocytes','Endothelial.cells','Excitatory.neurons','Inhibitory.neurons',
               'Microglia','Oligodendrocytes','OPCs...COPs','Pericytes')
pheno_all0 = read.csv('data/110_mean.csv')
pheno_all0 = pheno_all0[grep('FA', pheno_all0$phenocode),]
pheno_all = pheno_all0$assoc_files
pheno_all = gsub('ukbiobank_mean_','ukb_phase1to3_dti441_dec21_2019_', pheno_all)
pheno_all = gsub('_aug2020.zip','.fastGWA', pheno_all)


merge.data = data.frame(stringsAsFactors = F)
for (pheno in pheno_all){
  for (tissue in tissue_all){
    file_name = paste('results/PMR-Egger/CellToIDP/CellToWM/', pheno,'.',tissue, '.PMR-Egger.txt', sep='')
    new.data = read.table(file_name, header = T, sep = "\t", stringsAsFactors = F)
    merge.data = rbind(merge.data, new.data)
  }
  cat(pheno, '\n')
}
write.table(merge.data, "results/PMR-Egger/PMR-Egger.CellToWM.all.txt", sep="\t", row.names=F, col.names=T, quote=F)



#Merge SMR results++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rm(list =ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)

tissue_all = c('Astrocytes','Endothelial.cells','Excitatory.neurons','Inhibitory.neurons',
               'Microglia','Oligodendrocytes','OPCs...COPs','Pericytes')
pheno_all0 = read.csv('data/110_mean.csv')
pheno_all0 = pheno_all0[grep('FA', pheno_all0$phenocode),]
pheno_all = pheno_all0$assoc_files
pheno_all = gsub('ukbiobank_mean_','ukb_phase1to3_dti441_dec21_2019_', pheno_all)
pheno_all = gsub('_aug2020.zip','.fastGWA', pheno_all)

merge.data = data.frame(stringsAsFactors = F)
for (pheno in pheno_all){
  for (tissue in tissue_all){
    file_name = paste('results/SMR/CellToIDP/CellToWM/', pheno,'.',tissue,'.smr', sep='')
    new.data = fread(file_name, header = T, sep = "\t", stringsAsFactors = F, data.table = F)
    if (nrow(new.data) == 0) next

    new.data$pheno = pheno
    new.data$tissue = tissue
    new.data = new.data[,c(ncol(new.data)-1, ncol(new.data), 1:(ncol(new.data)-2))]
    merge.data = rbind(merge.data, new.data)
  }
  cat(pheno, '\n')
}
write.table(merge.data, "results/SMR/SMR.CellToWM.all.txt", sep="\t", row.names=F, col.names=T, quote=F)



#Merge TwoSampleMR results++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rm(list =ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)

tissue_all = c('Astrocytes','Endothelial.cells','Excitatory.neurons','Inhibitory.neurons',
               'Microglia','Oligodendrocytes','OPCs...COPs','Pericytes')
pheno_all0 = read.csv('data/110_mean.csv')
pheno_all0 = pheno_all0[grep('FA', pheno_all0$phenocode),]
pheno_all = pheno_all0$assoc_files
pheno_all = gsub('ukbiobank_mean_','ukb_phase1to3_dti441_dec21_2019_', pheno_all)
pheno_all = gsub('_aug2020.zip','.fastGWA', pheno_all)


merge.data = data.frame(stringsAsFactors = F)
for (pheno in pheno_all){
  for (tissue in tissue_all){
    file_name = paste('results/TwoSampleMR/CellToIDP/CellToWM/', pheno,'.',tissue,'.TwoSampleMR.txt', sep='')
    new.data = read.table(file_name, header = T, sep = "\t", stringsAsFactors = F)
    
    merge.data = rbind(merge.data, new.data)
  }
  cat(pheno, '\n')
}
write.table(merge.data, "results/TwoSampleMR/TwoSampleMR.CellToWM.all.txt", sep="\t", row.names=F, col.names=T, quote=F)



#Merge GSMR results++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rm(list =ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)

tissue_all = c('Astrocytes','Endothelial.cells','Excitatory.neurons','Inhibitory.neurons',
               'Microglia','Oligodendrocytes','OPCs...COPs','Pericytes')
pheno_all0 = read.csv('data/110_mean.csv')
pheno_all0 = pheno_all0[grep('FA', pheno_all0$phenocode),]
pheno_all = pheno_all0$assoc_files
pheno_all = gsub('ukbiobank_mean_','ukb_phase1to3_dti441_dec21_2019_', pheno_all)
pheno_all = gsub('_aug2020.zip','.fastGWA', pheno_all)

merge.data = data.frame(stringsAsFactors = F)
for (pheno in pheno_all){
  for (tissue in tissue_all){
    file_name = paste('results/GSMR/CellToIDP/CellToWM/', pheno,'.',tissue, '.GSMR.txt', sep='')
    new.data = read.table(file_name, header = T, sep = "\t", stringsAsFactors = F)
    merge.data = rbind(merge.data, new.data)
  }
  cat(pheno, 'is OK.\n')
}
write.table(merge.data, paste('results/GSMR/GSMR.CellToWM.all.txt',sep=''), sep="\t", row.names=F, col.names=T, quote=F)
