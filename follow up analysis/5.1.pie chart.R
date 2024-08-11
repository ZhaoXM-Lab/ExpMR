rm(list=ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)
library(stats)

f1 = fread('results/signif_genes.CellToIDP.fdr.txt', data.table = F)
f1$TraitCategory3 = 'Brain regional volume'
f1[grep('FA', f1$pheno), 'TraitCategory3'] = 'White matter microstructure'


f2 = fread('results/signif_genes.CellToDisorder.fdr.txt', data.table = F)
trait_mapping = read.csv('gwas_meta_all.csv', header = T)
trait_mapping = trait_mapping[, c('TraitCategory1','TraitCategory2','TraitAbbr')]
f2$TraitCategory1 = trait_mapping[match(f2$pheno, trait_mapping$TraitAbbr), 'TraitCategory1']
f2$TraitCategory2 = trait_mapping[match(f2$pheno, trait_mapping$TraitAbbr), 'TraitCategory2']


IDP_all = unique(f1$TraitCategory3)
disorder_all = unique(f2$pheno)
cell_all = unique(c(f1$tissue, f2$tissue))

result0 = data.frame(stringsAsFactors = F)
for (IDP in IDP_all){
  for (disorder in disorder_all){
    for (cell in cell_all){
      sub_f1 = subset(f1, f1$TraitCategory3 == IDP & f1$tissue == cell)
      sub_f2 = subset(f2, f2$pheno == disorder & f2$tissue == cell)
      
      gene_list1 = sub_f1$gene_name
      gene_list2 = sub_f2$gene_name
      inter_gene = intersect(gene_list1, gene_list2)
      if (length(inter_gene) == 0) p_hyper = 1
      
      if (length(inter_gene) > 0){
        q = length(inter_gene)
        m = length(gene_list1)
        n = 26861 - m
        k = length(gene_list2)
        p_hyper = phyper(q-1, m, n, k, lower.tail=F)
      }
      
      result0.tmp = data.frame(IDP = IDP, Disorder = disorder, Cell = cell, Causal_genes1 = paste(gene_list1,collapse=','),
                               Causal_genes2 = paste(gene_list2, collapse=','),
                               Overlap_genes = paste(inter_gene,collapse=','), p_hyper = p_hyper, stringsAsFactors = F)
      result0 = rbind(result0, result0.tmp)
      cat(IDP, 'and', disorder, 'is OK.', '\n') 
    }
  }
}
result0$fdr = p.adjust(result0$p_hyper, method = 'fdr', n = nrow(result0))
result = subset(result0, result0$p_hyper < 1)


#画饼图--------------------
library(ggplot2)
library(jjPlot)
df = subset(result, result$fdr < 0.05)
df = df[, c("IDP","Disorder","Cell")]
df$per = 1
library(reshape2)
df2 = reshape2::dcast(data=df, Disorder + IDP ~ Cell, value.var = 'per')
df2[is.na(df2)] = 0

df3 = df2
df3[,3:ncol(df3)] <- df2[,3:ncol(df2)] %>% dplyr::mutate(across(everything(), ~./rowSums(across(everything()))))
df3$group = 1:nrow(df3)


df3.long = reshape2::melt(df3, id.vars = c('Disorder','IDP','group'),
               variable.name = 'Cell', value.name = 'per')
df3.long$Disorder = factor(df3$Disorder, levels = c("BIP","SCZ","ASD","PTSD","ADHD","ANX","MDD","TS",
                                                    "AD","MS","INS","SC","DPW","GRT","ASP","NEU","EA"))

df3.long$Disorder = factor(df3$Disorder, levels = c("BIP","SCZ","ASD","PTSD","ADHD","AUD","ANX","MDD","TS", #PD
                                                    "PD","MS","AD","IS","EPI","INS",   #ND
                                                    "SC","GRT","DPW","ASP","NEU","INT","EA")) #BCP

# mapping with type
my_colors = c("#F0027F","#BEAED4","#FDC086","#FFFF99","#386CB0","#7FC97F","#BF5B17")
ggplot(df3.long, aes(x = Disorder, y = IDP, group = group)) +
  geom_jjPointPie(aes(pievar = per,fill = Cell)) +
  #coord_fixed()+
  labs(x = "", y = "") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = my_colors) 



