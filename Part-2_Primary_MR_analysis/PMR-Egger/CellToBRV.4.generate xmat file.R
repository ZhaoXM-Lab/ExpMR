rm(list =ls())
setwd('/public/home/yanganyi/Desktop/genexpMR/data/gene_expression/Malhotra_lab')
library(data.table)

file_name = list.files('mr_format/PMR-Egger/')
file_name = file_name[grep('.signif_variant_gene_pairs_5e-8.txt', file_name)]
tissues = c('Astrocytes','Endothelial.cells','Excitatory.neurons','Inhibitory.neurons',
            'Microglia','Oligodendrocytes','OPCs...COPs','Pericytes')
  
for (i in 1:length(tissues)){
  tissue = tissues[i]
  command = paste('/public/home/yanganyi/Desktop/Software/gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static', 
                  ' --bfile /public/home/yanganyi/Desktop/MAGMA_Annotation_file/g1000_eur_ref_file/g1000_eur',
                  ' --extract /public/home/yanganyi/Desktop/genexpMR/data/gene_expression/Malhotra_lab/allele/',tissue,'.allele',
                  ' --recode --out /public/home/yanganyi/Desktop/genexpMR/data/gene_expression/Malhotra_lab/xmat/',tissue, sep = '')
  system(command, intern=F)
}
