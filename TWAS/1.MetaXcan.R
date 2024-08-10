cd /share/home1/yanganyi/Desktop/Software/MetaXcan/software/
conda activate imlabtools
cd /share/inspurStorage/home1/yanganyi/Desktop/

  
#101 ukb brain volumes & DLPFC-------------------------------------------------------
rm(list = ls())
setwd("/share/inspurStorage/home1/yanganyi/Desktop/genexpMR/")

#gwas summary file
path_gwas = "/share/inspurStorage/home1/yanganyi/Desktop/Data/gwas_sumstat/yyc_postgwas/ukb/BrainVolume/"
pheno_all = list.files(path_gwas) 

for (i in 1:length(pheno_all)){
  pheno = pheno_all[i]
  
  command = paste('python /share/inspurStorage/home1/yanganyi/Desktop/Software/MetaXcan/software/SPrediXcan.py', 
                  ' --model_db_path /share/inspurStorage/home1/yanganyi/Desktop/Data/TWAS/DLPFC/DLPFC_newMetax.db',
                  ' --covariance /share/inspurStorage/home1/yanganyi/Desktop/Data/TWAS/DLPFC/DLPFC.cov.txt.gz',
                  ' --gwas_folder /share/inspurStorage/home1/yanganyi/Desktop/Data/gwas_sumstat/yyc_postgwas/ukb/BrainVolume',
                  ' --gwas_file_pattern "', pheno, '"',
                  ' --snp_column SNP',
                  ' --effect_allele_column A1',
                  ' --non_effect_allele_column A2',
                  ' --beta_column BETA',
                  ' --pvalue_column P',
                  ' --output_file /share/inspurStorage/home1/yanganyi/Desktop/genexpMR/results/TWAS/',pheno,'.DLPFC.MetaXcan.csv',
                  sep ='')
  system(command, intern=F)
}


#22 WMMP & DLPFC-------------------------------------------------------
rm(list = ls())
setwd("/public/home/yanganyi/Desktop/")

#gwas summary file
pheno_all = list.files('/public/home/yanganyi/Desktop/Data/gwas_sumstat/zhuhongtu/')
pheno_all = pheno_all[grep('.fastGWA', pheno_all)]

for (i in 1:length(pheno_all)){
  pheno = pheno_all[i]
  
  command = paste('python /public/home/yanganyi/Desktop/Software/MetaXcan/software/SPrediXcan.py', 
                  ' --model_db_path /public/home/yanganyi/Desktop/Data/TWAS/DLPFC/DLPFC_newMetax.db',
                  ' --covariance /public/home/yanganyi/Desktop/Data/TWAS/DLPFC/DLPFC.cov.txt.gz',
                  ' --gwas_folder /public/home/yanganyi/Desktop/Data/gwas_sumstat/zhuhongtu',
                  ' --gwas_file_pattern "', pheno, '"',
                  ' --snp_column SNP',
                  ' --effect_allele_column A1',
                  ' --non_effect_allele_column A2',
                  ' --beta_column BETA',
                  ' --pvalue_column P',
                  ' --output_file /public/home/yanganyi/Desktop/genexpMR/results/TWAS/',pheno,'.DLPFC.MetaXcan.csv',
                  sep ='')
  system(command, intern=F)
}



#26 phenos & GTEx_brain_tissues-------------------------------------------------------
rm(list = ls())
setwd("/share/inspurStorage/home1/yanganyi/Desktop/genexpMR/")

trait_info = read.csv('data/Nsample.csv', header = T, stringsAsFactors = F)
tissue_all = trait_info[which(trait_info$Category == 'GTEx'), 'Trait']
pheno_all = trait_info[which(trait_info$Category == 'GWAS'), 'Trait']


for (i in 1:length(pheno_all)){
  for (j in 1:length(tissue_all)){
    pheno = pheno_all[i]
    tissue = tissue_all[j]
    
    command = paste('python /share/inspurStorage/home1/yanganyi/Desktop/Software/MetaXcan/software/SPrediXcan.py', 
                    ' --model_db_path /share/inspurStorage/home1/yanganyi/Desktop/Data/TWAS/GTEx-V8_elastic_net_eqtls_models/en_',tissue,'.db',
                    ' --covariance /share/inspurStorage/home1/yanganyi/Desktop/Data/TWAS/GTEx-V8_elastic_net_eqtls_models/en_',tissue,'.txt.gz',
                    ' --gwas_folder /share/inspurStorage/home1/yanganyi/Desktop/genexpMR/data/gwas_summary/mr_format/SMR',
                    ' --gwas_file_pattern "',pheno,'.format.txt"',
                    ' --snp_column SNP',
                    ' --effect_allele_column A1',
                    ' --non_effect_allele_column A2',
                    ' --beta_column b',
                    ' --pvalue_column p',
                    ' --output_file /share/inspurStorage/home1/yanganyi/Desktop/genexpMR/results/TWAS/',pheno,'.',tissue,'.MetaXcan.csv',
                    sep ='')
    system(command, intern=F)
    cat(pheno, '-', tissue, '\n')
  }
}


#26 phenos & GTEx_blood-------------------------------------------------------
cd /public/home/yanganyi/Desktop/Software
git clone https://github.com/hakyimlab/MetaXcan
cd /public/home/yanganyi/Desktop/Software/MetaXcan/software
conda env create -f conda_env.yaml
conda activate imlabtools

ssh bnode01 -X
cd /public/home/yanganyi/Desktop/Software/MetaXcan/software
conda activate imlabtools


rm(list = ls())
setwd("/public/home/yanganyi/Desktop/genexpMR/")

trait_info = read.csv('data/Nsample.csv', header = T, stringsAsFactors = F)
tissue_all = c('Whole_Blood')
pheno_all = trait_info[which(trait_info$Category == 'GWAS'), 'Trait']

for (i in 1:length(pheno_all)){
  for (j in 1:length(tissue_all)){
    pheno = pheno_all[i]
    tissue = tissue_all[j]
    
    command = paste('python /public/home/yanganyi/Desktop/Software/MetaXcan/software/SPrediXcan.py', 
                    ' --model_db_path /public/home/yanganyi/Desktop/Data/TWAS/GTEx-V8_elastic_net_eqtls_models/en_',tissue,'.db',
                    ' --covariance /public/home/yanganyi/Desktop/Data/TWAS/GTEx-V8_elastic_net_eqtls_models/en_',tissue,'.txt.gz',
                    ' --gwas_folder /public/home/yanganyi/Desktop/genexpMR/data/gwas_summary/mr_format/SMR/brain_disorder',
                    ' --gwas_file_pattern "',pheno,'.format.txt"',
                    ' --snp_column SNP',
                    ' --effect_allele_column A1',
                    ' --non_effect_allele_column A2',
                    ' --beta_column b',
                    ' --pvalue_column p',
                    ' --output_file /public/home/yanganyi/Desktop/genexpMR/results/TWAS/',pheno,'.',tissue,'.MetaXcan.csv',
                    sep ='')
    system(command, intern=F)
    cat(pheno, '-', tissue, '\n')
  }
}



#101 phenos & 13 GTEx_tissues-------------------------------------------------------
ssh bnode01 -X
cd /public/home/yanganyi/Desktop/Software/MetaXcan/software
conda activate imlabtools


rm(list = ls())
setwd("/public/home/yanganyi/Desktop/genexpMR/")

trait_info = read.csv('data/Nsample.csv', header = T, stringsAsFactors = F)
tissue_all = trait_info[which(trait_info$Category == 'GTEx'), 'Trait']
path_gwas = "/public/home/yanganyi/Desktop/Data/gwas_sumstat/yyc_postgwas/ukb/BrainVolume/" #绝对路径最前面要加'/'
pheno_all = list.files(path_gwas) 

for (i in 1:length(pheno_all)){
  for (j in 1:length(tissue_all)){
    pheno = pheno_all[i]
    tissue = tissue_all[j]
    
    command = paste('python /public/home/yanganyi/Desktop/Software/MetaXcan/software/SPrediXcan.py', 
                    ' --model_db_path /public/home/yanganyi/Desktop/Data/TWAS/GTEx-V8_elastic_net_eqtls_models/en_',tissue,'.db',
                    ' --covariance /public/home/yanganyi/Desktop/Data/TWAS/GTEx-V8_elastic_net_eqtls_models/en_',tissue,'.txt.gz',
                    ' --gwas_folder /public/home/yanganyi/Desktop/Data/gwas_sumstat/yyc_postgwas/ukb/BrainVolume',
                    ' --gwas_file_pattern "', pheno,'"',
                    ' --snp_column SNP',
                    ' --effect_allele_column A1',
                    ' --non_effect_allele_column A2',
                    ' --beta_column BETA',
                    ' --pvalue_column P',
                    ' --output_file /public/home/yanganyi/Desktop/genexpMR/results/TWAS/101IDP_13GTEx/',pheno,'.',tissue,'.MetaXcan.csv',
                    sep ='')
    
    system(command, intern=F)
    cat(pheno, '-', tissue, '\n')
  }
}


#22 WMMP & 13 GTEx_tissues-------------------------------------------------------
ssh bnode01 -X
cd /public/home/yanganyi/Desktop/Software/MetaXcan/software
conda activate imlabtools
export PATH=$PATH:/public/software/apps/R/4.2.0/bin


rm(list = ls())
setwd("/public/home/yanganyi/Desktop/genexpMR/")

trait_info = read.csv('data/Nsample.csv', header = T, stringsAsFactors = F)
tissue_all = trait_info[which(trait_info$Category == 'GTEx'), 'Trait']
pheno_all = list.files('/public/home/yanganyi/Desktop/Data/gwas_sumstat/zhuhongtu/')
pheno_all = pheno_all[grep('.fastGWA', pheno_all)]

for (i in 1:length(pheno_all)){
  for (j in 1:length(tissue_all)){
    pheno = pheno_all[i]
    tissue = tissue_all[j]
    
    command = paste('python /public/home/yanganyi/Desktop/Software/MetaXcan/software/SPrediXcan.py', 
                    ' --model_db_path /public/home/yanganyi/Desktop/Data/TWAS/GTEx-V8_elastic_net_eqtls_models/en_',tissue,'.db',
                    ' --covariance /public/home/yanganyi/Desktop/Data/TWAS/GTEx-V8_elastic_net_eqtls_models/en_',tissue,'.txt.gz',
                    ' --gwas_folder /public/home/yanganyi/Desktop/Data/gwas_sumstat/zhuhongtu',
                    ' --gwas_file_pattern "', pheno,'"',
                    ' --snp_column SNP',
                    ' --effect_allele_column A1',
                    ' --non_effect_allele_column A2',
                    ' --beta_column BETA',
                    ' --pvalue_column P',
                    ' --output_file /public/home/yanganyi/Desktop/genexpMR/results/TWAS/22WM_13GTEx/',pheno,'.',tissue,'.MetaXcan.csv',
                    sep ='')
    
    system(command, intern=F)
    cat(pheno, '-', tissue, '\n')
  }
}


#ubuntu: 26 phenos & DLPFC-------------------------------------------------------
conda activate imlabtools


rm(list = ls())
setwd("/home/u/桌面/temporary")

trait_info = read.csv('Nsample.csv', header = T, stringsAsFactors = F)
tissue_all = c('DLPFC')
pheno_all = trait_info[which(trait_info$Category == 'GWAS'), 'Trait']

for (i in 1:length(pheno_all)){
  for (j in 1:length(tissue_all)){
    pheno = pheno_all[i]
    tissue = tissue_all[j]
    
    command = paste('python /home/u/桌面/Software/MetaXcan/software/SPrediXcan.py', 
                    ' --model_db_path DLPFC_newMetax.db',
                    ' --covariance DLPFC.cov.txt.gz',
                    ' --gwas_folder brain_disorder',
                    ' --gwas_file_pattern "',pheno,'.format.txt"',
                    ' --snp_column SNP',
                    ' --effect_allele_column A1',
                    ' --non_effect_allele_column A2',
                    ' --beta_column b',
                    ' --pvalue_column p',
                    ' --output_file /home/u/桌面/temporary/results/',pheno,'.',tissue,'.MetaXcan.csv',
                    sep ='')
    system(command, intern=F)
    cat(pheno, '-', tissue, '\n')
  }
}
