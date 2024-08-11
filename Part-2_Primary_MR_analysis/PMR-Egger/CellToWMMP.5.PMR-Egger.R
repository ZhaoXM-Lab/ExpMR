rm(list =ls())
setwd('/public/home/yanganyi/Desktop/genexpMR/')

library(data.table)
library(PMR)
library(doParallel)
library(foreach)


run_mr = function(mr_file10, mr_file20, xmat0, n1, n2, j){
  library(PMR)

  pheno1 = unique(mr_file10$Gene)[j]
  iv = mr_file10[mr_file10$Gene == pheno1, 'SNP']
    

  mr_file1 = subset(mr_file10, mr_file10$Gene == pheno1 & mr_file10$SNP %in% iv)  
  mr_file2 = subset(mr_file20, mr_file20$SNP %in% iv) 
  xmat = xmat0[, which(colnames(xmat0) %in% iv)]
    
  
  try_error <- try({

    mr_file1 = mr_file1[match(iv, mr_file1$SNP), ]  
    mr_file2 = mr_file2[match(iv, mr_file2$SNP), ]
    xmat = xmat[, match(iv, colnames(xmat))]
    
    Zscore_1 = mr_file1$Zscore
    Zscore_2 = mr_file2$Zscore
    LDmatrix = cor(xmat)^2
    
    mr.result = PMR_summary_Egger(Zscore_1, Zscore_2, LDmatrix, LDmatrix, samplen1=n1, samplen2=n2, 
                                  lambda=0.15, max_iterin = 5000, epsin=1e-5, Heritability_geneexpression_threshold=0)
    mr.result = as.data.frame(mr.result)
    iv = paste(colnames(LDmatrix), collapse=",") 
    result.tmp = data.frame(tissue=tissue, pheno1=pheno1, pheno2=pheno2,
                            causal_effect=mr.result$causal_effect, causal_pvalue=mr.result$causal_pvalue,
                            pleiotropy_effect=mr.result$pleiotropy_effect, pleiotropy_pvalue=mr.result$pleiotropy_pvalue,
                            iv=iv, n1=n1, n2=n2)
  }, silent = T)
  if ('try-error' %in% class(try_error)) {
    result.tmp = data.frame(stringsAsFactors = F)
  }
  return(result.tmp)
}



#main--------------------------------------------------------------------------------
tissue_all = c('Astrocytes','Endothelial.cells','Excitatory.neurons','Inhibitory.neurons',
               'Microglia','Oligodendrocytes','OPCs...COPs','Pericytes')
file_names = list.files('/public/home/yanganyi/Desktop/Data/gwas_sumstat/zhuhongtu/')
pheno2_all = file_names[grep('.fastGWA', file_names)]

for (pheno2 in pheno2_all){
  for (tissue in tissue_all){

    file_name1 = paste('data/gene_expression/Malhotra_lab/mr_format/PMR-Egger/SNP_in_outcome_clumped_txt/',tissue,'.',pheno2,'.clumped.txt',sep='')
    file10 = fread(file_name1, data.table = F)
    
    xmat0 = fread(paste('data/gene_expression/Malhotra_lab/xmat/',tissue,'.xmat.gz',sep=''), data.table=F, header=T)
    xmat0 = xmat0[-1, ]
    xmat0 = xmat0[, -which(names(xmat0) %in% c('FID','IID'))] 
    xmat0 = apply(xmat0, 2, as.numeric)
    
    file_name2 = paste('data/gwas_summary/mr_format/PMR-Egger/IDP/', pheno2,'.format.txt',sep='')
    file20 = fread(file_name2, data.table = F)

    file20 = subset(file20, !(is.na(file20$Zscore)))


    snp_all = Reduce(intersect, list(colnames(xmat0), file10$SNP, file20$SNP))
    mr_file10 = subset(file10, file10$SNP %in% snp_all)
    mr_file20 = subset(file20, file20$SNP %in% snp_all)
    xmat0 = xmat0[, which(colnames(xmat0) %in% snp_all)]
    
    n1 = 192  
    n2 = 33292  
    
    
    cl = makeCluster(10) 
    registerDoParallel(cl)
    result = foreach(j = 1:length(unique(mr_file10$Gene)), .combine = 'rbind') %dopar% run_mr(mr_file10, mr_file20, xmat0, n1, n2, j)
    stopCluster(cl) #-----
    
    write.table(result, paste('results/PMR-Egger/CellToIDP/',pheno2,'.',tissue,'.PMR-Egger.txt',sep=''), row.names = F, col.names= T, sep="\t", quote = F)
    cat(tissue,' and ', pheno2, 'is OK.\n')
  }
}