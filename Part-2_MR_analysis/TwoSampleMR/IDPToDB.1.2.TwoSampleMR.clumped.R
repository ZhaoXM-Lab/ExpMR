rm(list =ls())
setwd('/public/home/yanganyi/Desktop/genexpMR/')
library(data.table)


tissue_all = list.files('data/gwas_summary/mr_format/SMR/IDP/')
tissue_all = tissue_all[grep('.EUR.signif_pairs_5e-8.txt', tissue_all)]
tissue_all = gsub('.EUR.signif_pairs_5e-8.txt', '', tissue_all)

trait_info = read.csv('/public/home/yanganyi/Desktop/genexpMR/data/Nsample.csv', header = T, stringsAsFactors = F)
outcome_all = trait_info[which(trait_info$Category == 'GWAS'), 'Trait']


for (outcome in outcome_all){

  filename = paste('data/gwas_summary/mr_format/SMR/brain_disorder/',outcome,'.format.txt', sep = '')
  file1 =  fread(filename, data.table = F) 
  
  for (tissue in tissue_all){
    filename = paste('data/gwas_summary/mr_format/SMR/IDP/', tissue,'.EUR.signif_pairs_5e-8.txt',sep='')
    file2 = fread(filename, data.table = F)
    
    file2 = subset(file2, file2$SNP %in% file1$SNP)
    
    #if (nrow(file2) == 0) next
    outfilename = paste('data/gene_expression/GTEx/mr_format/TwoSampleMR/SNP_in_outcome/',tissue,'.',outcome,'.SNP_interc.txt',sep='')
    write.table(file2, outfilename, row.names = F, col.names= T, sep="\t", quote = F)
    cat(tissue,' and ', outcome, 'is OK.\n')
  }
}



####################################################
rm(list =ls())
setwd('/public/home/yanganyi/Desktop/genexpMR/')

tissue_all = list.files('data/gwas_summary/mr_format/SMR/IDP/')
tissue_all = tissue_all[grep('.EUR.signif_pairs_5e-8.txt', tissue_all)]
tissue_all = gsub('.EUR.signif_pairs_5e-8.txt', '', tissue_all)
trait_info = read.csv('/public/home/yanganyi/Desktop/genexpMR/data/Nsample.csv', header = T, stringsAsFactors = F)
outcome_all = trait_info[which(trait_info$Category == 'GWAS'), 'Trait']


tissue = 'ukb_roi_volume_may12_2019_phase1and2_pheno12_allchr_withA2'; outcome = 'INT'
for (outcome in outcome_all){
  for (tissue in tissue_all[]){
    command = paste('/public/home/yanganyi/Desktop/Software/plink1/plink ',
                    ' --bfile /public/home/yanganyi/Desktop/MAGMA_Annotation_file/g1000_eur_ref_file/g1000_eur',
                    ' --clump-p1 1',
                    ' --clump-r2 0.1',
                    ' --clump-kb 250',
                    ' --clump /public/home/yanganyi/Desktop/genexpMR/data/gene_expression/GTEx/mr_format/TwoSampleMR/SNP_in_outcome/',tissue,'.',outcome,'.SNP_interc.txt',
                    ' --clump-snp-field SNP',
                    ' --clump-field pval',
                    ' --out /public/home/yanganyi/Desktop/genexpMR/data/gene_expression/GTEx/mr_format/TwoSampleMR/SNP_in_outcome_clumped/',tissue,'.',outcome,
                    sep='')
    system(command, intern=F)
  }
}



###################################################
rm(list =ls())
setwd('/public/home/yanganyi/Desktop/genexpMR/')
library(data.table)
library(tidyverse)

tissue_all = list.files('data/gwas_summary/mr_format/SMR/IDP/')
tissue_all = tissue_all[grep('.EUR.signif_pairs_5e-8.txt', tissue_all)]
tissue_all = gsub('.EUR.signif_pairs_5e-8.txt', '', tissue_all)
trait_info = read.csv('/public/home/yanganyi/Desktop/genexpMR/data/Nsample.csv', header = T, stringsAsFactors = F)
outcome_all = trait_info[which(trait_info$Category == 'GWAS'), 'Trait']


for (outcome in outcome_all){
  for (tissue in tissue_all){
    if (!(file.exists(paste('data/gene_expression/GTEx/mr_format/TwoSampleMR/SNP_in_outcome/', tissue,'.',outcome,'.SNP_interc.txt',sep='')))) next
    file = fread(paste('data/gene_expression/GTEx/mr_format/TwoSampleMR/SNP_in_outcome/', tissue,'.',outcome,'.SNP_interc.txt',sep=''), data.table = F) 
    
    if (!(file.exists(paste('data/gene_expression/GTEx/mr_format/TwoSampleMR/SNP_in_outcome_clumped/', tissue,'.',outcome, '.clumped',sep='')))) next
    clumped = fread(paste('data/gene_expression/GTEx/mr_format/TwoSampleMR/SNP_in_outcome_clumped/', tissue,'.',outcome, '.clumped',sep=''), data.table = F) #原来的文件样子
    
    file  = subset(file, file$SNP %in% clumped$SNP)
    outname = paste('data/gene_expression/GTEx/mr_format/TwoSampleMR/SNP_in_outcome_clumped_txt/', tissue,'.',outcome,'.clumped.txt',sep='')
    write.table(file, outname, sep=' ', row.names=F,quote=F)
    
    cat(tissue,' and ', outcome, 'is OK.\n')
  }
}
