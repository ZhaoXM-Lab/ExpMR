#Merge all DB TWAS result files#########################################################################################
rm(list =ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)

trait_info = read.csv('data/Nsample.csv', header = T, stringsAsFactors = F)
tissue_all = trait_info[which(trait_info$Category == 'GTEx'), 'Trait']
pheno_all = trait_info[which(trait_info$Category == 'GWAS'), 'Trait']

pheno = pheno_all[1]; tissue = tissue_all[1]
merge.data = data.frame(stringsAsFactors = F)
for (pheno in pheno_all){
  for (tissue in tissue_all){
    file_name = paste('results/TWAS/', pheno,'.',tissue, '.MetaXcan.csv', sep='')
    new.data = read.table(file_name, header = T, sep = ",", stringsAsFactors = F)
    new.data$pheno = pheno
    new.data$tissue = tissue
    new.data = new.data[,c(ncol(new.data)-1,ncol(new.data),1:(ncol(new.data)-2))]
    merge.data = rbind(merge.data, new.data)
  }
  cat(pheno, 'is OK.\n')
}
write.table(merge.data, paste('results/TWAS/','TWAS.all.txt',sep=''), sep="\t", row.names=F, col.names=T, quote=F)



#Merge all IDP TWAS result files#########################################################################################
rm(list =ls())
setwd('F:/类脑/1_实操文件/genexpMR/results/TWAS')
library(data.table)

files = list.files() 
files = files[grep('ukb',files)]

merge.data = data.frame(stringsAsFactors = F)
for (i in 1:length(files)){
  pheno = strsplit(files[i], "[.]")[[1]][1] 
  
  new.data = read.table(files[i], header = T, sep = ",", stringsAsFactors = F)
  new.data$pheno = pheno
  new.data$tissue = 'DLPFC'
  new.data = new.data[,c(ncol(new.data)-1,ncol(new.data),1:(ncol(new.data)-2))]
  
  merge.data = rbind(merge.data, new.data)
  cat(pheno, 'is OK.\n')
}
write.table(merge.data, 'TWAS.DLPFC.all.txt', sep="\t", row.names=F, col.names=T, quote=F)


#Bar chart
data = subset(merge.data, merge.data$pvalue < 0.05)
calcu_Ngene = function(data){
  Ngene = data.frame(stringsAsFactors = F)
  for (pheno_idx in unique(data$pheno)){
    for (tissue_idx in unique(data$tissue)){
      subset_file = subset(data, data$tissue == tissue_idx & data$pheno == pheno_idx)
      N = nrow(subset_file)
      Ngene.tmp = data.frame(pheno = pheno_idx, tissue = tissue_idx, N = N)
      Ngene = rbind(Ngene.tmp, Ngene)
    }
  }
  return(Ngene)
}
Ngene = calcu_Ngene(data)
Ngene = Ngene[order(Ngene$N, decreasing = F),]
Ngene$pheno = factor(Ngene$pheno, levels = unique(Ngene$pheno))

ggplot(Ngene, aes(x = pheno, y = N)) +
  geom_bar(stat = "identity", colour = "black", fill = "#FED9A6") +
  theme(axis.text.x = element_text(size=15, angle=90, hjust=1, vjust=1), 
        axis.text.y = element_text(size=12), 
        legend.title = element_text(), legend.position =  c(0.85, 0.1)) +
  ylab('Number of signature genes') +
  xlab('')  
  #coord_flip()




#Count nominal signif genes###########################################
rm(list =ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)

file = fread('results/TWAS/TWAS.all.txt', data.table = F)
tissue_all = unique(file$tissue)
pheno_all = unique(file$pheno)


Ngene = data.frame(stringsAsFactors = F)
for (pheno_idx in pheno_all){
  for (tissue_idx in tissue_all){
    subset_file = subset(file, file$tissue == tissue_idx & file$pheno == pheno_idx)
    N = nrow(subset_file)
    Ngene.tmp = data.frame(pheno_name = pheno_idx, tissue = tissue_idx, N = N)
    Ngene = rbind(Ngene.tmp, Ngene)
  }
}

#Heatmap
library(pheatmap)
data = Ngene[, c("pheno_name","tissue","N")]
matrix = matrix(data$N, length(unique(Ngene$tissue)), length(unique(Ngene$pheno_name)))
rownames(matrix) = unique(Ngene$tissue)
colnames(matrix) = unique(Ngene$pheno_name)
pheatmap(t(matrix), display_numbers = TRUE, number_format = "%.0f",fontsize_number = 10, number_color = 'black',
         cluster_row = FALSE, cluster_col = FALSE)



#Count nominal siginif genes###########################################
gene_list0 = read.table('data/gene_loc/gencode.v26.GRCh38.genes.loc', header = T, sep = '\t', stringsAsFactors = F)
#gene_list0 = subset(gene_list0, gene_list0$transcript_type ==' lincRNA')

file = subset(file, file$gene %in% gene_list0$gene_id)
signif_genes = data.frame(stringsAsFactors = F)
for (pheno_idx in pheno_all){
  for (tissue_idx in tissue_all){
    subset_file = subset(file, file$tissue == tissue_idx & file$pheno == pheno_idx)
    
    gene_num.1 =  length(which(subset_file$pvalue < 0.05))
    
    signif_genes.tmp = data.frame(pheno_name = pheno_idx, tissue = tissue_idx, gene_num.1 = gene_num.1)
    signif_genes = rbind(signif_genes.tmp, signif_genes)
  }
}

#Heatmap
library(pheatmap)
data = signif_genes[, c("pheno_name","tissue","gene_num.1")]
matrix = matrix(data$gene_num.1, length(unique(data$tissue)), length(unique(data$pheno_name)))
rownames(matrix) = unique(data$tissue)
colnames(matrix) = unique(data$pheno_name)
matrix = t(matrix)
pheatmap(matrix, display_numbers = TRUE, number_format = "%.0f",fontsize_number = 10, number_color = 'black',
         cluster_row = FALSE, cluster_col = FALSE)