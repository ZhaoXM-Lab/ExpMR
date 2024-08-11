#Merge PMR-Egger results++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rm(list =ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)

trait_info = read.csv('data/Nsample.csv', header = T, stringsAsFactors = F)
tissue_all = trait_info[which(trait_info$Category == 'GWAS'), 'Trait']
pheno_all = paste('ukb_roi_volume_may12_2019_phase1and2_pheno', 1:101,'_allchr_withA2', sep = '')

merge.data = data.frame(stringsAsFactors = F)
for (pheno in pheno_all){
  for (tissue in tissue_all){
    file_name = paste('results/PMR-Egger/DisorderToIDP/', pheno,'.',tissue, '.PMR-Egger.txt', sep='')
    if (!file.exists(file_name)) next
        
    try_error <- try({
      new.data = fread(file_name, header = T, sep = "\t", stringsAsFactors = F)
      merge.data = rbind(merge.data, new.data)
      })
    if ('try-error' %in% class(try_error)) {
      next
    }
  }
}
write.table(merge.data, "results/PMR-Egger/PMR-Egger.DisorderToIDP.all.txt", sep="\t", row.names=F, col.names=T, quote=F)



#Merge TwoSampleMR results++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rm(list =ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)

trait_info = read.csv('data/Nsample.csv', header = T, stringsAsFactors = F)
tissue_all = trait_info[which(trait_info$Category == 'GWAS'), 'Trait']
pheno_all = paste('ukb_roi_volume_may12_2019_phase1and2_pheno', 1:101,'_allchr_withA2', sep = '')


merge.data = data.frame(stringsAsFactors = F)
for (i in 1:length(pheno_all)){
  for (j in 1:length(tissue_all)){
    pheno = pheno_all[i]
    tissue = tissue_all[j]
    file_name = paste('results/TwoSampleMR/DisorderToIDP/DisorderToVolume/', pheno,'.',tissue,'.TwoSampleMR.txt', sep='')
    if (!file.exists(file_name)) next
    
    try_error <- try({
      new.data = fread(file_name, header = T, sep = "\t", stringsAsFactors = F)
      merge.data = rbind(merge.data, new.data)
    })
    if ('try-error' %in% class(try_error)) {
      next
    }
  }
  cat(i, '\n')
}
write.table(merge.data, "results/TwoSampleMR/TwoSampleMR.DisorderToIDP.all.txt", sep="\t", row.names=F, col.names=T, quote=F)



#Merge GSMR results++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rm(list =ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)

trait_info = read.csv('data/Nsample.csv', header = T, stringsAsFactors = F)
tissue_all = trait_info[which(trait_info$Category == 'GWAS'), 'Trait']
pheno_all = paste('ukb_roi_volume_may12_2019_phase1and2_pheno', 1:101,'_allchr_withA2', sep = '')


merge.data = data.frame(stringsAsFactors = F)
for (pheno in pheno_all){
  for (tissue in tissue_all){
    file_name = paste('results/GSMR/DisorderToIDP/', pheno,'.',tissue, '.GSMR.txt', sep='')
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
write.table(merge.data, paste('results/GSMR/GSMR.DisorderToIDP.all.txt',sep=''), sep="\t", row.names=F, col.names=T, quote=F)


