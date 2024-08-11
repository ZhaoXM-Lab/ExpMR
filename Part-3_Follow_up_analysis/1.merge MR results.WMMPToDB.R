#Merge PMR-Egger results++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rm(list =ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)

tissue_all0 = read.csv('data/110_mean.csv')
tissue_all0 = tissue_all0[grep('FA', tissue_all0$phenocode),]
tissue_all = tissue_all0$assoc_files
tissue_all = gsub('ukbiobank_mean_','ukb_phase1to3_dti441_dec21_2019_', tissue_all)
tissue_all = gsub('_aug2020.zip','.fastGWA', tissue_all)
trait_info = read.csv('data/Nsample.csv', header = T, stringsAsFactors = F)
pheno_all = trait_info[which(trait_info$Category == 'GWAS'), 'Trait']

merge.data = data.frame(stringsAsFactors = F)
for (pheno in pheno_all){
  for (tissue in tissue_all){
    file_name = paste('results/PMR-Egger/IDPToDisorder/WMToDisorder/', pheno,'.',tissue, '.PMR-Egger.txt', sep='')
    if (!file.exists(file_name)) next
    
    try_error <- try({
      new.data = fread(file_name, header = T, sep = "\t", stringsAsFactors = F)
      merge.data = rbind(merge.data, new.data)
    })
    if ('try-error' %in% class(try_error)) next
  }
  cat(pheno, 'is OK.\n')
}
write.table(merge.data, "results/PMR-Egger/PMR-Egger.WMToDisorder.all.txt", sep="\t", row.names=F, col.names=T, quote=F)



#Merge TwoSmapleMR results++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rm(list =ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)

tissue_all0 = read.csv('data/110_mean.csv')
tissue_all0 = tissue_all0[grep('FA', tissue_all0$phenocode),]
tissue_all = tissue_all0$assoc_files
tissue_all = gsub('ukbiobank_mean_','ukb_phase1to3_dti441_dec21_2019_', tissue_all)
tissue_all = gsub('_aug2020.zip','.fastGWA', tissue_all)
trait_info = read.csv('data/Nsample.csv', header = T, stringsAsFactors = F)
pheno_all = trait_info[which(trait_info$Category == 'GWAS'), 'Trait']

merge.data = data.frame(stringsAsFactors = F)
for (pheno in pheno_all){
  for (tissue in tissue_all){
    file_name = paste('results/TwoSampleMR/IDPToDisorder/WMToDisorder/', pheno,'.',tissue,'.TwoSampleMR.txt', sep='')
    if (!file.exists(file_name)) next
    
    new.data = read.table(file_name, header = T, sep = "\t", stringsAsFactors = F)
    merge.data = rbind(merge.data, new.data)
  }
  cat(pheno, 'is OK.\n')
}
merge.data = subset(merge.data, !is.na(merge.data$outcome))
write.table(merge.data, "results/TwoSampleMR/TwoSampleMR.WMToDisorder.all.txt", sep="\t", row.names=F, col.names=T, quote=F)



#Merge GSMR results++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rm(list =ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)

tissue_all0 = read.csv('data/110_mean.csv')
tissue_all0 = tissue_all0[grep('FA', tissue_all0$phenocode),]
tissue_all = tissue_all0$assoc_files
tissue_all = gsub('ukbiobank_mean_','ukb_phase1to3_dti441_dec21_2019_', tissue_all)
tissue_all = gsub('_aug2020.zip','.fastGWA', tissue_all)
trait_info = read.csv('data/Nsample.csv', header = T, stringsAsFactors = F)
pheno_all = trait_info[which(trait_info$Category == 'GWAS'), 'Trait']

merge.data = data.frame(stringsAsFactors = F)
for (pheno in pheno_all){
  for (tissue in tissue_all){
    file_name = paste('results/GSMR/IDPToDisorder/WMToDisorder/', pheno,'.',tissue, '.GSMR.txt', sep='')
    if (!file.exists(file_name)) next
    
    new.data = read.table(file_name, header = T, sep = "\t", stringsAsFactors = F)
    merge.data = rbind(merge.data, new.data)
  }
  cat(pheno, 'is OK.\n')
}
merge.data = subset(merge.data, !is.na(merge.data$outcome))
write.table(merge.data, paste('results/GSMR/GSMR.WMToDisorder.all.txt',sep=''), sep="\t", row.names=F, col.names=T, quote=F)
