# Extract the genotype data from a Disorder using GCTA
rm(list =ls())
setwd('/public/home/yanganyi/Desktop/genexpMR/data/gene_expression/')
library(data.table)


trait_info = read.csv('/public/home/yanganyi/Desktop/genexpMR/data/Nsample.csv', header = T, stringsAsFactors = F)
tissues = trait_info[which(trait_info$Category == 'GWAS'), 'Trait']

for (i in 1:length(tissues)){
  tissue = tissues[i]
  command = paste('/public/home/yanganyi/Desktop/Software/gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static', 
                  ' --bfile /public/home/yanganyi/Desktop/MAGMA_Annotation_file/g1000_eur_ref_file/g1000_eur',
                  ' --extract /public/home/yanganyi/Desktop/genexpMR/data/gwas_summary/mr_format/PMR-Egger/brain_disorder/allele/',tissue,'.EUR.signif_pairs_5e-8.allele',
                  ' --recode --out /public/home/yanganyi/Desktop/genexpMR/data/gwas_summary/mr_format/PMR-Egger/brain_disorder/xmat/',tissue, sep = '')
  system(command, intern=F)
}