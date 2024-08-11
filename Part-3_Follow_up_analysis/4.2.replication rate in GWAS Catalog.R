rm(list=ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)


signif_gene = fread('results/signif_genes.CellToDisorder.fdr.txt',data.table = F)
signif_gene$gene_pheno = paste(signif_gene$gene_name, signif_gene$pheno, sep = '-')
find_gene_pheno = unique(signif_gene$gene_pheno)


#gwas catalog association file
gwas_catalog0 = fread('F:/类脑/2_数据报告资料/数据库文件/GWAS catalog/gwas_catalog_v1.0.2-associations_e107_r2023-06-09.tsv',
                     data.table = F,quote="")
gwas_catalog0 = gwas_catalog0[,c('STUDY', 'DISEASE/TRAIT', 'SNPS', 'MAPPED_GENE')]
length(unique(gwas_catalog0$`DISEASE/TRAIT`))


paired_pheno = read.table('data/gwas_catalog与26 brain disorder配对表型.txt', header = T, sep='\t', stringsAsFactors=F)

gwas_catalog = data.frame(stringsAsFactors = F)
for(charac in paired_pheno$gwas_catalog_charac_include){
  gwas_catalog.tmp = gwas_catalog0[grep(charac, gwas_catalog0$`DISEASE/TRAIT`), ]
  gwas_catalog = rbind(gwas_catalog, gwas_catalog.tmp)
  cat(charac, '\n')
}
gwas_catalog$id = paste(gwas_catalog$STUDY, gwas_catalog$`DISEASE/TRAIT`, gwas_catalog$MAPPED_GENE, sep=',') 
gwas_catalog = gwas_catalog[!duplicated(gwas_catalog$id), ]
gwas_catalog = gwas_catalog[, -which(names(gwas_catalog) == 'id')]
gwas_catalog = subset(gwas_catalog, !is.na(gwas_catalog$MAPPED_GENE))


signif_gene$If_repl = 'N'
for (i in 1:nrow(signif_gene)){
  pheno = signif_gene$pheno[i]
  gene = signif_gene$gene_name[i]

  charac_include = paired_pheno[paired_pheno$Abbr == pheno, 'gwas_catalog_charac_include']
  
  gwas_catalog_pheno = data.frame(stringsAsFactors = F)
  for (charac in charac_include){
    gwas_catalog_pheno.tmp = gwas_catalog[grep(charac, gwas_catalog$`DISEASE/TRAIT`), ]
    gwas_catalog_pheno = rbind(gwas_catalog_pheno, gwas_catalog_pheno.tmp)
  }

  gwas_catalog_pheno_gene = gwas_catalog_pheno[grep(gene, gwas_catalog_pheno$MAPPED_GENE), ]
  
  if (nrow(gwas_catalog_pheno_gene) == 0) next
  if (nrow(gwas_catalog_pheno_gene) > 0) signif_gene$If_repl[i] = 'Y'
  
  if (i %% 10 == 0) cat(i, '\n')
}
length(unique(signif_gene[which(signif_gene$If_repl == 'Y'), 'gene_pheno']))
length(unique(signif_gene[which(signif_gene$If_repl == 'Y'), 'gene_pheno'])) / length(unique(signif_gene$gene_pheno))


write.table(signif_gene, 'results/signif_genes.repGWAS_catalog.CellToDisorder.fdr.txt', sep="\t", row.names=F, col.names=T, quote=F)


gwas_catalog_uni = gwas_catalog
for(i in 1:nrow(paired_pheno)){
  charac = paired_pheno[i,'gwas_catalog_charac_include']
  abbr = paired_pheno[i,'Abbr']
  gwas_catalog_uni$`DISEASE/TRAIT` = gsub(charac, abbr, gwas_catalog_uni$`DISEASE/TRAIT`)
  cat(charac, '\n')
}  
gwas_catalog_uni$gene_pheno = paste(gwas_catalog_uni$`DISEASE/TRAIT`, gwas_catalog_uni$MAPPED_GENE, sep='')


length(unique(gwas_catalog_uni$gene_pheno));length(unique(signif_gene[which(signif_gene$replicated == TRUE), 'gene_pheno']));length(unique(signif_gene$gene_pheno))
q = length(unique(signif_gene[which(signif_gene$replicated == TRUE), 'gene_pheno']))
m = length(unique(signif_gene$gene_pheno))
n = 26861*length(unique(signif_gene$pheno)) - m
k = length(unique(gwas_catalog_uni$gene_pheno))
phyper(q-1, m, n, k, lower.tail=F)

