rm(list =ls())
setwd('/public/home/yanganyi/Desktop/genexpMR/')
library(data.table)


trait_info = read.csv('data/Nsample.csv', header = T, stringsAsFactors = F)
tissue_all = trait_info[which(trait_info$Category == 'GWAS'), 'Trait']

outcome_all = list.files('/public/home/yanganyi/Desktop/Data/gwas_sumstat/zhuhongtu/')
outcome_all = outcome_all[grep('.fastGWA', outcome_all)]


for (outcome in outcome_all){
  filename = paste('/public/home/yanganyi/Desktop/genexpMR/data/gwas_summary/mr_format/PMR-Egger/IDP/',outcome, '.format.txt', sep = '')
  file1 =  fread(filename, data.table = F) 
  
  for (tissue in tissue_all){
    filename = paste('data/gwas_summary/mr_format/PMR-Egger/brain_disorder/', tissue,'.EUR.signif_pairs_5e-8.txt',sep='')
    file2 = fread(filename, data.table = F)
    
    file2 = subset(file2, file2$SNP %in% file1$SNP)
    
    outfilename = paste('data/gene_expression/GTEx/mr_format/PMR-Egger/SNP_in_outcome/',tissue,'.',outcome,'.SNP_interc.txt',sep='')
    write.table(file2, outfilename, row.names = F, col.names= T, sep="\t", quote = F)
    cat(tissue,' and ', outcome, 'is OK.\n')
  }
}



##########################################################################################
rm(list =ls())
setwd('/public/home/yanganyi/Desktop/genexpMR/')

trait_info = read.csv('/public/home/yanganyi/Desktop/genexpMR/data/Nsample.csv', header = T, stringsAsFactors = F)
tissue_all = trait_info[which(trait_info$Category == 'GWAS'), 'Trait']

outcome_all = list.files('/public/home/yanganyi/Desktop/Data/gwas_sumstat/zhuhongtu/')
outcome_all = outcome_all[grep('.fastGWA', outcome_all)]

for (outcome in outcome_all){
  for (tissue in tissue_all){
    command = paste('/public/home/yanganyi/Desktop/Software/plink1/plink ',
                    ' --bfile /public/home/yanganyi/Desktop/MAGMA_Annotation_file/g1000_eur_ref_file/g1000_eur',
                    ' --clump-p1 1',
                    ' --clump-r2 0.6',
                    ' --clump-kb 250',
                    ' --clump /public/home/yanganyi/Desktop/genexpMR/data/gene_expression/GTEx/mr_format/PMR-Egger/SNP_in_outcome/',tissue,'.',outcome,'.SNP_interc.txt',
                    ' --clump-snp-field SNP',
                    ' --clump-field pval',
                    ' --out /public/home/yanganyi/Desktop/genexpMR/data/gene_expression/GTEx/mr_format/PMR-Egger/SNP_in_outcome_clumped/',tissue,'.',outcome,
                    sep='')
    system(command, intern=F)
  }
}


##############################################
rm(list =ls())
setwd('/public/home/yanganyi/Desktop/genexpMR/data/gene_expression/GTEx/mr_format/PMR-Egger/')
library(data.table)
library(tidyverse)

trait_info = read.csv('/public/home/yanganyi/Desktop/genexpMR/data/Nsample.csv', header = T, stringsAsFactors = F)
tissue_all = trait_info[which(trait_info$Category == 'GWAS'), 'Trait']

outcome_all = list.files('/public/home/yanganyi/Desktop/Data/gwas_sumstat/zhuhongtu/')
outcome_all = outcome_all[grep('.fastGWA', outcome_all)]


for (outcome in outcome_all){
  for (tissue in tissue_all){
    file = fread(paste('SNP_in_outcome/', tissue,'.',outcome,'.SNP_interc.txt',sep=''), data.table = F) 
    clumped = fread(paste('SNP_in_outcome_clumped/', tissue,'.',outcome, '.clumped',sep=''), data.table = F)
    
    file  = subset(file, file$SNP %in% clumped$SNP)
    
    file = file[,-which(names(file)=='pval')]
    
    outname = paste('SNP_in_outcome_clumped_txt/', tissue,'.',outcome,'.clumped.txt',sep='')
    write.table(file, outname, sep=' ', row.names=F,quote=F)
    
    cat(tissue,' and ', outcome, 'is OK.\n')
  }
}
