rm(list =ls())
setwd('/public/home/yanganyi/Desktop/genexpMR/')
library(data.table)

tissue_all = c('Astrocytes','Endothelial.cells','Excitatory.neurons','Inhibitory.neurons',
               'Microglia','Oligodendrocytes','OPCs...COPs','Pericytes')
pheno_all = list.files('/share/home1/yanganyi/Desktop/Data/gwas_sumstat/yyc_postgwas/ukb/BrainVolume/')
pheno_all = gsub('.txt', '', pheno_all)


indiv = function(j){
  pheno = pheno_all[i]
  tissue = tissue_all[j]
  command = paste('/public/home/yanganyi/Desktop/Software/smr-1.3.1-linux-x86_64/smr-1.3.1', 
                  ' --bfile /public/home/yanganyi/Desktop/MAGMA_Annotation_file/g1000_eur_ref_file/g1000_eur',
                  ' --gwas-summary /public/home/yanganyi/Desktop/genexpMR/data/gwas_summary/mr_format/SMR/IDP/',pheno,'.format.txt',
                  ' --beqtl-summary /public/home/yanganyi/Desktop/genexpMR/data/gene_expression/Malhotra_lab/mr_format/SMR/',tissue,
                  ' --out /public/home/yanganyi/Desktop/genexpMR/results/SMR/CellToIDP/', pheno,'.',tissue, 
                  ' --thread-num 20', sep = '')
  system(command, intern=F)
  return(0)
}

library(doParallel)
library(foreach)

for (i in 1:length(pheno_all)){
  pheno = pheno_all[i]

  cl = makeCluster(2) 
  registerDoParallel(cl)
  clusterExport(cl, 'i') 
  result = foreach(j = 1:length(tissue_all), .combine = 'rbind') %dopar% indiv(j)
  stopCluster(cl)
  
  cat(i, 'th file has finished')
}
