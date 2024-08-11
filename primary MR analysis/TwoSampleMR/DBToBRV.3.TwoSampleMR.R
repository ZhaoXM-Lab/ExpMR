rm(list =ls())
setwd('/public/home/yanganyi/Desktop/genexpMR/')
library(data.table)
library(TwoSampleMR)
library(doParallel)
library(foreach)



run_tissue = function(tissue, outcome_data, outcome){
  library(TwoSampleMR)
  filename = paste('data/gene_expression/GTEx/mr_format/TwoSampleMR/SNP_in_outcome_clumped_txt/',tissue,
                   '.ukb_roi_volume_may12_2019_phase1and2_pheno100_allchr_withA2.clumped.txt',sep='')
  exposure_data.allGE = read_exposure_data(filename, clump = FALSE, sep = " ",
                                           phenotype_col = "gene", snp_col = "SNP", beta_col = "beta",
                                           se_col = "se", eaf_col = "eaf",
                                           effect_allele_col = "effect_allele",
                                           other_allele_col = "other_allele", pval_col = "pval", 
                                           samplesize_col = "samplesize")

  exposure_data.allGE = subset(exposure_data.allGE, exposure_data.allGE$SNP %in% outcome_data$SNP)
  
  
  GE = unique(exposure_data.allGE$exposure)[1]
  exposure_data = subset(exposure_data.allGE, exposure_data.allGE$exposure == GE)

  dat = harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)
  

  result = mr(dat)
  if (nrow(result) == 1 & result$nsnp[1] == 1){
    result[, c('Q','Q_df','Q_pval','egger_intercept','se_egger_intercept','pval_egger_intercept')] = NA
  } else{
    heterogeneity = mr_heterogeneity(dat)
    pleiotropy = mr_pleiotropy_test(dat)
    
    result = merge(result, heterogeneity[,c('method','Q','Q_df','Q_pval')], by = 'method', all.x = T)
    result$egger_intercept = pleiotropy$egger_intercept
    result$se_egger_intercept = pleiotropy$se
    result$pval_egger_intercept = pleiotropy$pval
  }

  if (nrow(result) != 0) result$iv = paste(dat$SNP[dat$mr_keep == TRUE],collapse=',') 
  result$tissue = tissue
  result$outcome = outcome
  result = result[, c('tissue','exposure','outcome','method','nsnp','b','se','pval','iv',
                      'Q','Q_df','Q_pval','egger_intercept','se_egger_intercept','pval_egger_intercept')]
  
  outfilename = paste('results/TwoSampleMR/DisorderToIDP/',outcome,'.',tissue,'.TwoSampleMR.txt',sep='')
  write.table(result, outfilename, row.names = F, col.names= T, sep="\t", quote = F)
  
  return(0)
}


#Main-------------------------------------------------------------------
trait_info = read.csv('data/Nsample.csv', header = T, stringsAsFactors = F)
tissue_all = trait_info[which(trait_info$Category == 'GWAS'), 'Trait']

outcome_all = list.files('/public/home/yanganyi/Desktop/Data/gwas_sumstat/yyc_postgwas/ukb/BrainVolume/')
outcome_all = gsub('.txt', '', outcome_all)

for (outcome in outcome_all){
  filename = paste('data/gwas_summary/mr_format/SMR/IDP/',outcome,'.format.txt', sep = '')
  outcome_data = read_outcome_data(filename, snps = NULL, sep = "\t",
                                   snp_col = "SNP", beta_col = "b",
                                   se_col = "se", eaf_col = "freq",
                                   effect_allele_col = "A1",
                                   other_allele_col = "A2", pval_col = "p",
                                   samplesize_col = "n")
  
  for (tissue in tissue_all){
    result_all = run_tissue(tissue, outcome_data, outcome)
    cat(outcome,'and', tissue, 'is OK.', '\n')
  }
}

