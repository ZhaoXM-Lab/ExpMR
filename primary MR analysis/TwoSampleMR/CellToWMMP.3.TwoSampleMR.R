rm(list =ls())
setwd('/public/home/yanganyi/Desktop/genexpMR/')
library(data.table)
library(TwoSampleMR)
library(doParallel)
library(foreach)


run_mr = function(exposure_data.allGE, outcome_data, i){
  library(TwoSampleMR)
  GE = unique(exposure_data.allGE$exposure)[i]
  exposure_data = subset(exposure_data.allGE, exposure_data.allGE$exposure == GE)

  dat = harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)
  
  #MR analysis
  res.tmp = mr(dat)
  if (nrow(res.tmp) == 0){
    res.tmp = data.frame(matrix(rep(NA, 15), nrow = 1))
    names(res.tmp) = c("id.exposure","id.outcome","outcome","exposure","method","nsnp","b","se","pval",'Q','Q_df',
                       'Q_pval','egger_intercept','se_egger_intercept','pval_egger_intercept')
  } else if (nrow(res.tmp) == 1 & res.tmp$nsnp[1] == 1){
    res.tmp[, c('Q','Q_df','Q_pval','egger_intercept','se_egger_intercept','pval_egger_intercept')] = NA
  } else{
    heterogeneity = mr_heterogeneity(dat)
    pleiotropy = mr_pleiotropy_test(dat)
    
    res.tmp = merge(res.tmp, heterogeneity[,c('method','Q','Q_df','Q_pval')], by = 'method', all.x = T)
    res.tmp$egger_intercept = pleiotropy$egger_intercept
    res.tmp$se_egger_intercept = pleiotropy$se
    res.tmp$pval_egger_intercept = pleiotropy$pval
  }
  
  if (nrow(res.tmp) != 0) res.tmp$iv = paste(dat$SNP[dat$mr_keep == TRUE],collapse=',') 
  
  return(res.tmp)
}



#Main-------------------------------------------------------------------
tissue_all = c('Astrocytes','Endothelial.cells','Excitatory.neurons','Inhibitory.neurons',
               'Microglia','Oligodendrocytes','OPCs...COPs','Pericytes')
file_names = list.files('/public/home/yanganyi/Desktop/Data/gwas_sumstat/zhuhongtu/')
outcome_all = file_names[grep('.fastGWA', file_names)]


for (outcome in outcome_all){
  filename = paste('data/gwas_summary/mr_format/SMR/IDP/',outcome,'.format.txt', sep = '')
  outcome_data = read_outcome_data(filename, snps = NULL, sep = "\t",
                                   snp_col = "SNP", beta_col = "b",
                                   se_col = "se", eaf_col = "freq",
                                   effect_allele_col = "A1",
                                   other_allele_col = "A2", pval_col = "p",
                                   samplesize_col = "n")

  for (tissue in tissue_all){
    cat(tissue,' and ', outcome, 'starts.\n')
    
    filename = paste('data/gene_expression/Malhotra_lab/mr_format/TwoSampleMR/SNP_in_outcome_clumped_txt/',tissue,'.',outcome,'.clumped.txt',sep='')
    exposure_data.allGE = read_exposure_data(filename, clump = FALSE, sep = " ",
                                       phenotype_col = "gene", snp_col = "SNP", beta_col = "beta",
                                       se_col = "se", eaf_col = "eaf",
                                       effect_allele_col = "effect_allele",
                                       other_allele_col = "other_allele", pval_col = "pval", 
                                       samplesize_col = "samplesize")
    exposure_data.allGE = subset(exposure_data.allGE, exposure_data.allGE$SNP %in% outcome_data$SNP)
    
    
    cl = makeCluster(10) 
    registerDoParallel(cl)
    result = foreach(i = 1:length(unique(exposure_data.allGE$exposure)), .combine = 'rbind') %dopar% run_mr(exposure_data.allGE, outcome_data, i)
    stopCluster(cl)
    
    result$tissue = tissue
    result$outcome = outcome
    result = result[, c('tissue','exposure','outcome','method','nsnp','b','se','pval','Q','Q_df',
                        'Q_pval','egger_intercept','se_egger_intercept','pval_egger_intercept','iv')]
    
    outfilename = paste('results/TwoSampleMR/CellToIDP/CellToWM/', outcome, '.', tissue,'.TwoSampleMR.txt',sep='')
    write.table(result, outfilename, row.names = F, col.names= T, sep="\t", quote = F)
    
    cat(tissue,' and ', outcome, 'is OK.\n')
  }
}

