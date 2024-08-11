rm(list =ls())
setwd('/public/home/yanganyi/Desktop/genexpMR/')
library(data.table)

tissue_all = c('Astrocytes','Endothelial.cells','Excitatory.neurons','Inhibitory.neurons',
               'Microglia','Oligodendrocytes','OPCs...COPs','Pericytes')
file_names = list.files('/public/home/yanganyi/Desktop/Data/gwas_sumstat/zhuhongtu/')
outcome_all = file_names[grep('.fastGWA', file_names)]


outcome = outcome_all[1]; tissue = tissue_all[1]
for (outcome in outcome_all){

  filename = paste('data/gwas_summary/mr_format/PMR-Egger/IDP/',outcome,'.format.txt', sep = '')
  file1 =  fread(filename, data.table = F) 
  
  for (tissue in tissue_all){
    filename = paste('data/gene_expression/Malhotra_lab/mr_format/PMR-Egger/', tissue,'.signif_variant_gene_pairs_5e-8.txt',sep='')
    file2 = fread(filename, data.table = F)
    
    file2 = subset(file2, file2$SNP %in% file1$SNP) 
    file2$pval = ifelse(file2$Zscore < 0, 2*pnorm(file2$Zscore), 2-2*pnorm(file2$Zscore))
      
    outfilename = paste('data/gene_expression/Malhotra_lab/mr_format/PMR-Egger/SNP_in_outcome/',tissue,'.',outcome,'.SNP_interc.txt',sep='')
    write.table(file2, outfilename, row.names = F, col.names= T, sep="\t", quote = F)
    cat(tissue,' and ', outcome, 'is OK.\n')
  }
}



############################################################################################
rm(list =ls())
setwd('/public/home/yanganyi/Desktop/genexpMR/')

tissue_all = c('Astrocytes','Endothelial.cells','Excitatory.neurons','Inhibitory.neurons',
               'Microglia','Oligodendrocytes','OPCs...COPs','Pericytes')
file_names = list.files('/public/home/yanganyi/Desktop/Data/gwas_sumstat/zhuhongtu/')
outcome_all = file_names[grep('.fastGWA', file_names)]

for (outcome in outcome_all){ 
  for (tissue in tissue_all){
    command = paste('/public/home/yanganyi/Desktop/Software/plink1/plink ',
                    ' --bfile /public/home/yanganyi/Desktop/MAGMA_Annotation_file/g1000_eur_ref_file/g1000_eur',
                    ' --clump-p1 1',
                    ' --clump-r2 0.9',
                    ' --clump-kb 250',
                    ' --clump /public/home/yanganyi/Desktop/genexpMR/data/gene_expression/Malhotra_lab/mr_format/PMR-Egger/SNP_in_outcome/',tissue,'.',outcome,'.SNP_interc.txt',
                    ' --clump-snp-field SNP',
                    ' --clump-field pval',
                    ' --out /public/home/yanganyi/Desktop/genexpMR/data/gene_expression/Malhotra_lab/mr_format/PMR-Egger/SNP_in_outcome_clumped/',tissue,'.',outcome,
                    sep='')
    system(command, intern=F)
  }
}



#############################################################################
rm(list =ls())
setwd('/public/home/yanganyi/Desktop/genexpMR/data/gene_expression/Malhotra_lab/mr_format/PMR-Egger/')
library(data.table)

tissue_all = c('Astrocytes','Endothelial.cells','Excitatory.neurons','Inhibitory.neurons',
               'Microglia','Oligodendrocytes','OPCs...COPs','Pericytes')
file_names = list.files('/public/home/yanganyi/Desktop/Data/gwas_sumstat/zhuhongtu/')
outcome_all = file_names[grep('.fastGWA', file_names)]

for (outcome in outcome_all){ 
  for (tissue in tissue_all){
    file = fread(paste('SNP_in_outcome/', tissue,'.',outcome, '.SNP_interc.txt', sep=''), data.table = F) 
    clumped = fread(paste('SNP_in_outcome_clumped/', tissue,'.',outcome, '.clumped',sep=''), data.table = F) 
    
    file  = subset(file, file$SNP %in% clumped$SNP)
    
    file = file[,-which(names(file)=='pval')]
    
    outname = paste('SNP_in_outcome_clumped_txt/', tissue,'.',outcome, '.clumped.txt', sep='')
    write.table(file, outname, sep=' ', row.names=F,quote=F)
    
    cat(tissue,' and ', outcome, 'is OK.\n')
  }
}


