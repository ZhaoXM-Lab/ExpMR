########################################################################################
rm(list=ls())
setwd('/public/home/yanganyi/Desktop/genexpMR/')
library("gsmr")
library(data.table)
library(doParallel)
library(foreach)


gsmr_analysis = function(gsmr_data, xmat, pheno2, tissue, pheno1){
  snp_coeff_id = colnames(xmat)
  ldrho = cor(xmat) # Calculate the LD correlation matrix
  colnames(ldrho) = rownames(ldrho) = colnames(xmat)
  
  #Standardization
  snpfreq = gsmr_data$a1_freq             # allele frequencies of the SNPs
  bzx = gsmr_data$bzx     # effects of the instruments on risk factor
  bzx_se = gsmr_data$bzx_se       # standard errors of bzx
  bzx_n = gsmr_data$bzx_n          # GWAS sample size for the risk factor
  std_zx = std_effect(snpfreq, bzx, bzx_se, bzx_n)    # perform standardisation
  gsmr_data$std_bzx = std_zx$b    # standardized bzx
  gsmr_data$std_bzx_se = std_zx$se    # standardized bzx_se
  
  #GSMR analysis
  bzx = gsmr_data$bzx    # SNP effects on the risk factor
  bzx_se = gsmr_data$bzx_se    # standard errors of bzx
  bzx_pval = gsmr_data$bzx_pval   # p-values for bzx
  bzy = gsmr_data$bzy    # SNP effects on the disease
  bzy_se = gsmr_data$bzy_se    # standard errors of bzy
  bzy_pval = gsmr_data$bzy_pval    # p-values for bzy
  n_ref = 503    # Sample size of the reference sample
  gwas_thresh = 1    # GWAS threshold to select SNPs as the instruments for the GSMR analysis
  single_snp_heidi_thresh = 0.05    # p-value threshold for single-SNP-based HEIDI-outlier analysis
  multi_snp_heidi_thresh = 0.05    # p-value threshold for multi-SNP-based HEIDI-outlier analysis
  nsnps_thresh = 2   # the minimum number of instruments required for the GSMR analysis
  heidi_outlier_flag = T    # flag for HEIDI-outlier analysis
  ld_r2_thresh = 0.9    # LD r2 threshold to remove SNPs in high LD
  ld_fdr_thresh = 1   # FDR threshold to remove the chance correlations between the SNP instruments
  gsmr2_beta = 0     # 0 - the original HEIDI-outlier method; 1 - the new HEIDI-outlier method that is currently under development
  
  # GSMR analysis
  try_error <- try({
    gsmr_results = gsmr(bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval, ldrho, snp_coeff_id, n_ref, heidi_outlier_flag, gwas_thresh, 
                        single_snp_heidi_thresh, multi_snp_heidi_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr2_beta)
  })
  
  if (!('try-error' %in% class(try_error))){
    tmp.result = data.frame(outcome = pheno2, tissue = tissue, exposure = pheno1, bxy = gsmr_results$bxy,
                            bxy_se = gsmr_results$bxy_se, bxy_pval = gsmr_results$bxy_pval, 
                            n_iv = length(gsmr_results$used_index), n_outliers = length(gsmr_results$pleio_snps), 
                            iv = paste(gsmr_data[gsmr_results$used_index,'SNP'],collapse = ','),
                            stringsAsFactors = F)
  }else{
    tmp.result = data.frame(outcome = NA, tissue = NA, exposure = NA, bxy = NA, bxy_se = NA, 
                            bxy_pval = NA, n_iv = NA, n_outliers = NA, iv=NA, stringsAsFactors = F)
  }
  return(tmp.result)
}


indiv = function(j){
  library("gsmr")
  
  pheno1 = unique(gsmr_data0$gene)[j]

  iv = unique(gsmr_data0[gsmr_data0$gene == pheno1, 'SNP'])
  
  #提取iv对应的信息
  gsmr_data = subset(gsmr_data0, gsmr_data0$gene == pheno1 & gsmr_data0$SNP %in% iv)  
  gsmr_data = gsmr_data[,-which(names(gsmr_data)=='gene')]
  if(length(iv) > 1){
    xmat = xmat0[, which(colnames(xmat0) %in% iv)]

    gsmr_data = gsmr_data[match(iv, gsmr_data$SNP), ]
    xmat = xmat[, match(iv, colnames(xmat))]
    
    tmp.result = gsmr_analysis(gsmr_data, xmat, pheno2, tissue, pheno1)
    
  } else {
    tmp.result = data.frame(outcome = NA, tissue = NA, exposure = NA, bxy = NA, bxy_se = NA, 
                            bxy_pval = NA, n_iv = NA, n_outliers = NA, iv=NA, stringsAsFactors = F)
  }
  return(tmp.result)
}



#main--------------------------------------------------------------------------------
tissue_all = c('Astrocytes','Endothelial.cells','Excitatory.neurons','Inhibitory.neurons',
               'Microglia','Oligodendrocytes','OPCs...COPs','Pericytes')
pheno2_all = list.files('/public/home/yanganyi/Desktop/Data/gwas_sumstat/zhuhongtu/')
pheno2_all = pheno2_all[grep('.fastGWA', pheno2_all)]


for (pheno2 in pheno2_all){
  for (tissue in tissue_all){
    file_name1 = paste('/public/home/yanganyi/Desktop/genexpMR/data/gene_expression/Malhotra_lab/mr_format/GSMR/',tissue,
                       '.signif_pairs_5e-8.txt',sep='')
    file10 = fread(file_name1, data.table = F)
    
    
    xmat0 = fread(paste('data/gene_expression/Malhotra_lab/xmat/',tissue,'.xmat.gz',sep=''), data.table=F, header=T)
    xmat0 = xmat0[-1, ]
    xmat0 = xmat0[, -which(names(xmat0) %in% c('FID','IID'))]
    xmat0 = apply(xmat0, 2, as.numeric)
    
    file_name2 = paste('data/gwas_summary/mr_format/GSMR/IDP/', pheno2,'.format.txt',sep='')
    file20 = fread(file_name2, data.table = F)
    
    snp_all = Reduce(intersect, list(colnames(xmat0), file10$SNP, file20$SNP))
    mr_file10 = subset(file10, file10$SNP %in% snp_all)
    mr_file20 = subset(file20, file20$SNP %in% snp_all)
    xmat0 = xmat0[, which(colnames(xmat0) %in% snp_all)]
    gsmr_data0 = merge(mr_file10, mr_file20, by = c('SNP','a1','a2','a1_freq'))
    

    cl = makeCluster(15) 
    registerDoParallel(cl)
    clusterExport(cl, c('xmat0','gsmr_data0','pheno2','tissue')) 
    result.all = foreach(j = 1:length(unique(gsmr_data0$gene)), .combine = 'rbind') %dopar% indiv(j)
    stopCluster(cl)
    
    cat(tissue,' and ', pheno2, 'is OK.\n')
    result.all = result.all[-which(is.na(result.all$outcome)), ]
    write.table(result.all, paste('results/GSMR/CellToIDP/',pheno2,'.',tissue,'.GSMR.txt',sep=''), row.names = F, col.names= T, sep="\t", quote = F)
  }
}
