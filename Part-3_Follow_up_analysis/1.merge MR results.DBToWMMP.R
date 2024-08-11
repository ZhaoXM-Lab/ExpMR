#Merge PMR-Egger results++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rm(list =ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)

trait_info = read.csv('data/Nsample.csv', header = T, stringsAsFactors = F)
tissue_all = trait_info[which(trait_info$Category == 'GWAS'), 'Trait']
pheno_all0 = read.csv('data/110_mean.csv')
pheno_all0 = pheno_all0[grep('FA', pheno_all0$phenocode),]
pheno_all = pheno_all0$assoc_files
pheno_all = gsub('ukbiobank_mean_','ukb_phase1to3_dti441_dec21_2019_', pheno_all)
pheno_all = gsub('_aug2020.zip','.fastGWA', pheno_all)

merge.data = data.frame(stringsAsFactors = F)
for (pheno in pheno_all){
  for (tissue in tissue_all){
    file_name = paste('results/PMR-Egger/DisorderToIDP/DisorderToWM/', pheno,'.',tissue, '.PMR-Egger.txt', sep='')
    if (!file.exists(file_name)) next
    
    try_error <- try({ 
      new.data = fread(file_name, header = T, sep = "\t", stringsAsFactors = F)
      merge.data = rbind(merge.data, new.data)
    })
    if ('try-error' %in% class(try_error)) {
      next
    }
  }
  cat(pheno, 'is OK.\n')
}
write.table(merge.data, "results/PMR-Egger/PMR-Egger.DisorderToWM.all.txt", sep="\t", row.names=F, col.names=T, quote=F)




#Merge TwoSmapleMR results++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rm(list =ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)

trait_info = read.csv('data/Nsample.csv', header = T, stringsAsFactors = F)
tissue_all = trait_info[which(trait_info$Category == 'GWAS'), 'Trait']
pheno_all0 = read.csv('data/110_mean.csv')
pheno_all0 = pheno_all0[grep('FA', pheno_all0$phenocode),]
pheno_all = pheno_all0$assoc_files
pheno_all = gsub('ukbiobank_mean_','ukb_phase1to3_dti441_dec21_2019_', pheno_all)
pheno_all = gsub('_aug2020.zip','.fastGWA', pheno_all)


merge.data = data.frame(stringsAsFactors = F)
for (pheno in pheno_all){
  for (tissue in tissue_all){
    file_name = paste('results/TwoSampleMR/DisorderToIDP/DisorderToWM/', pheno,'.',tissue,'.TwoSampleMR.txt', sep='')
    if (!file.exists(file_name)) next
    
    try_error <- try({
      new.data = fread(file_name, header = T, sep = "\t", stringsAsFactors = F)
      merge.data = rbind(merge.data, new.data)
    })
    if ('try-error' %in% class(try_error)) {
      next
    }
  }
  cat(pheno, 'is OK.\n')
}
merge.data = subset(merge.data, !is.na(merge.data$outcome))
write.table(merge.data, "results/TwoSampleMR/TwoSampleMR.DisorderToWM.all.txt", sep="\t", row.names=F, col.names=T, quote=F)



#Merge GSMR results++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rm(list =ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)

trait_info = read.csv('data/Nsample.csv', header = T, stringsAsFactors = F)
tissue_all = trait_info[which(trait_info$Category == 'GWAS'), 'Trait']
pheno_all0 = read.csv('data/110_mean.csv')
pheno_all0 = pheno_all0[grep('FA', pheno_all0$phenocode),]
pheno_all = pheno_all0$assoc_files
pheno_all = gsub('ukbiobank_mean_','ukb_phase1to3_dti441_dec21_2019_', pheno_all)
pheno_all = gsub('_aug2020.zip','.fastGWA', pheno_all)

merge.data = data.frame(stringsAsFactors = F)
for (pheno in pheno_all){
  for (tissue in tissue_all){
    file_name = paste('results/GSMR/DisorderToIDP/DisorderToWM/', pheno,'.',tissue, '.GSMR.txt', sep='')
    if (!file.exists(file_name)) next
    
    try_error <- try({
      new.data = fread(file_name, header = T, sep = "\t", stringsAsFactors = F)
      merge.data = rbind(merge.data, new.data)
    })
    if ('try-error' %in% class(try_error)) next
  }
  cat(pheno, 'is OK.\n')
}
merge.data = subset(merge.data, !is.na(merge.data$outcome))
write.table(merge.data, paste('results/GSMR/GSMR.DisorderToWM.all.txt',sep=''), sep="\t", row.names=F, col.names=T, quote=F)


