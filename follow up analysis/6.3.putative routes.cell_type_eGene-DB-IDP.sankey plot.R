rm(list=ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)
library(stats)
library(tidyr)


f1 = fread('results/signif_genes.CellToIDP.fdr.txt', data.table = F)
f1 = f1[, c("pheno", "tissue","gene","gene_name","gene_biotype")]
f1$trait_type = 'IDP'

f2 = fread('results/signif_genes.CellToDisorder.fdr.txt', data.table = F)
f2 = f2[, c("pheno", "tissue","gene","gene_name","gene_biotype")]
trait_mapping = read.csv('gwas_meta_all.csv', header = T)
trait_mapping = trait_mapping[, c('TraitCategory1', 'TraitAbbr')]
f2$trait_type = trait_mapping[match(f2$pheno, trait_mapping$TraitAbbr), 'TraitCategory1']

dat = rbind(f1, f2)
trait_all = unique(dat$pheno)
cell_all = unique(dat$tissue)


causal_genes = data.frame(stringsAsFactors = F)
for (trait in trait_all){
  for (cell in cell_all){
    sub_dat = subset(dat, dat$pheno == trait & dat$tissue == cell)
    if (nrow(sub_dat) == 0) next
    causal_gene_list = paste(sub_dat$gene, collapse=',')
    causal_genes.tmp = data.frame(trait = trait, cell = cell, causal_gene = causal_gene_list, 
                                  trait_type = sub_dat$trait_type[1], stringsAsFactors = F)
    causal_genes = rbind(causal_genes, causal_genes.tmp)
    
    cat(trait, 'and', cell, 'is OK.', '\n')
  }
}



result = data.frame(stringsAsFactors = F)
for (i in 1:nrow(causal_genes)){
  trait1 = causal_genes[i,'trait']
  trait1_type = causal_genes[i,'trait_type']
  cell1 = causal_genes[i,'cell']
  gene_list1 = unlist(strsplit(causal_genes[i,'causal_gene'], split = ","))
  n1 = length(gene_list1)
  
  for (j in 1:nrow(causal_genes)){
    trait2 = causal_genes[j,'trait']
    trait2_type = causal_genes[j,'trait_type']
    cell2 = causal_genes[j,'cell']
    
    if (trait1 == trait2 | cell1 != cell2) next
    gene_list2 = unlist(strsplit(causal_genes[j,'causal_gene'], split = ","))
    n2 = length(gene_list2)
    
    overlap_gene = intersect(gene_list1, gene_list2)
    if (length(overlap_gene) == 0) next

    result.tmp = data.frame(trait1 = trait1, trait1_type = trait1_type, trait2 = trait2, trait2_type = trait2_type, cell = cell1, 
                            overlap_gene = overlap_gene, n_overlap = length(overlap_gene))
    result = rbind(result, result.tmp)
    cat(i, j, 'is OK.','\n')
  }
}
result = subset(result, result$n_overlap > 0)
result$pairs = paste(result$trait1, result$trait2, sep = '-')



signif_pairs1 = fread('results/signif_IDPToDisorder.fdr.new.txt', sep="\t", stringsAsFactors = F)
overlap1 = intersect(result$pairs, signif_pairs1$IDPToDisorder)

signif_pairs2 = fread('results/signif_DisorderToIDP.fdr.new.txt', sep="\t", stringsAsFactors = F)
overlap2 = intersect(result$pairs, signif_pairs2$DisorderToIDP)


combin1 = paste(signif_pairs1$pheno, signif_pairs1$tissue, sep='-')
combin2 = paste(signif_pairs2$tissue, signif_pairs2$pheno, sep='-')
length(intersect(combin1, combin2))


path1 = subset(result, result$pairs %in% overlap1)
path1 = path1[, c('cell', 'overlap_gene','trait1', 'trait2')]
names(path1) = c('cell', 'gene', 'IDP', 'Disorder')


path2 = subset(result, result$pairs %in% overlap2)
path2 = path2[, c('cell', 'overlap_gene','trait1', 'trait2')]
names(path2) = c('cell', 'gene', 'Disorder', 'IDP')


gene_info = fread('F:/类脑/1_实操文件/genexpMR/data/gene_loc/Homo_sapiens.GRCh37.87.gtf.gz', data.table = F, skip = 5)
gene_info = separate(gene_info, col = V9, 
                     into = c("gene_id", "gene_version","gene_name","gene_source","gene_biotype"),sep = "; ")
gene_info = subset(gene_info, gene_info$V3 == 'gene')
gene_info = gene_info[, c('gene_name', 'V1', 'V4', 'V5','V7', 'gene_id')]
gene_info$gene_id = gsub("gene_id ", "", gene_info$gene_id)
gene_info$gene_id = gsub('"', '', gene_info$gene_id)
gene_info$gene_name = gsub("gene_name ", "", gene_info$gene_name)
gene_info$gene_name = gsub('"', '', gene_info$gene_name)


path1 = merge(path1, gene_info[,c('gene_name', 'gene_id')], by.x = 'gene', by.y = 'gene_id', all.x = T)
path1 = path1[,c('cell', 'gene_name', 'IDP', 'Disorder')]
path2 = merge(path2, gene_info[,c('gene_name', 'gene_id')], by.x = 'gene', by.y = 'gene_id', all.x = T)
path2 = path2[,c('cell', 'gene_name', 'Disorder', 'IDP')]

#IDP
path1$IDP = gsub('[.]', ' ', path1$IDP)
path2$IDP = gsub('[.]', ' ', path2$IDP)

#cell
cells = c('Pericytes', 'Inhibitory neuros', 'OPCs','Microglia', 'Endothelial', 'Excitatory neurons','Astrocytes', 'Oligodendrocytes')
path1$cell = gsub('OPCs...COPs', 'OPCs', path1$cell)
path1$cell = gsub('Excitatory.neurons', 'Excitatory neurons', path1$cell)
path1$cell = gsub('Inhibitory.neurons', 'Inhibitory neurons', path1$cell)
path2$cell = gsub('OPCs...COPs', 'OPCs', path2$cell)
path2$cell = gsub('Excitatory.neurons', 'Excitatory neurons', path2$cell)
path2$cell = gsub('Inhibitory.neurons', 'Inhibitory neurons', path2$cell)


path1$id = paste(path1$cell, path1$gene_name, path1$IDP, path1$Disorder, sep=' -> ')
path2$id = paste(path2$cell, path2$gene_name, path2$IDP, path2$Disorder, sep=' -> ')
bidirect = intersect(path1$id, path2$id)



#sankey plot---------------
library(ggsankey)
library(ggplot2)
library(cols4all)
#IDP -> Disorder#################################################################
dt = path1
dt$freq = 1
df = make_long(dt, cell, gene_name, IDP, Disorder, value = freq)

df$node <- factor(df$node, levels = c(dt$Disorder %>% unique() %>% rev(),
                                      dt$IDP %>% unique() %>% rev(), 
                                      dt$gene_name %>% unique() %>% rev(),
                                      dt$cell %>% unique() %>% rev()))

Nnode = apply(dt[,-ncol(dt)], 2, unique) %>% unlist %>% length  #所需要的颜色数量=节点个数
mycol <- c4a('rainbow_wh_rd', Nnode)
set.seed(2024)
mycol1 <- sample(mycol,length(mycol), replace = F) #随机打乱配色顺序(后面细胞类型的颜色统一一下)

names = c(df$node[df$node %in% dt$cell] %>% unique() %>% as.vector(),
          df$node[df$node %in% dt$gene_name] %>% unique() %>% as.vector(),
          df$node[df$node %in% dt$IDP] %>% unique() %>% as.vector(),
          df$node[df$node %in% dt$Disorder] %>% unique() %>% as.vector())
names(mycol1) = names


mycol1[which(names(mycol1) == 'Oligodendrocytes')] = '#386CB0'
mycol1[which(names(mycol1) == 'OPCs')] = '#4DC97F'
mycol1[which(names(mycol1) == 'Astrocytes')] = '#F0027F'
mycol1[which(names(mycol1) == 'Excitatory neurons')] = '#BEAED4'   
mycol1[which(names(mycol1) == 'Inhibitory neurons')] = '#FDC086'
mycol1[which(names(mycol1) == 'Microglia')] = '#FFFF99'

BRV = dt$IDP[!grepl('FA', dt$IDP)]
WMM = dt$IDP[grep('FA', dt$IDP)]
mycol1[which(names(mycol1) %in% BRV)] = '#4AB74C'
mycol1[which(names(mycol1) %in% WMM)] = '#0069E9'
BCP = c('CI', 'EA','INT','NEU','ASP','DPW','GRT','SC')
NEUD = c('INS','EPI','ICH','IS','AD','ALS','MS','PD')
PSYD = c('OCD','TS','ANX','MDD','ADHD','AUD','ASD','PTSD','BIP','SCZ')
mycol1[which(names(mycol1) %in% BCP)] = '#C58CA9'
mycol1[which(names(mycol1) %in% NEUD)] = '#E44D00'
mycol1[which(names(mycol1) %in% PSYD)] = '#FEB500'
mycol1[which(names(mycol1) %in% dt$gene_name)] = 'grey'


ggplot(df, aes(x = x,
               next_x = next_x,
               node = node,
               next_node = next_node,
               fill = node,
               label = node)) +
  scale_fill_manual(values = mycol1) + #更改配色
  geom_sankey(flow.alpha = 0.5, #条带不透明度
              smooth = 8, #条带弯曲度
              width = 0.12) + #节点宽度
  geom_sankey_text(size = 5, color = "black") +
  theme_void() +
  theme(legend.position = 'none')



#Disorder -> IDP###############################################################
dt = path2
dt$freq = 1
df = make_long(dt, cell, gene_name, Disorder, IDP, value = freq)
df$node <- factor(df$node, levels = c(dt$Disorder %>% unique()%>% rev(),
                                      dt$IDP %>% unique() %>% rev(), 
                                      dt$gene_name %>% unique() %>% rev(),
                                      dt$cell %>% unique() %>% rev()))


Nnode = apply(dt[,-ncol(dt)], 2, unique) %>% unlist %>% length
mycol <- c4a('rainbow_wh_rd', Nnode)
set.seed(2024)
mycol2 <- sample(mycol, length(mycol), replace = F)

names = c(df$node[df$node %in% dt$cell] %>% unique() %>% as.vector(),
          df$node[df$node %in% dt$gene_name] %>% unique() %>% as.vector(),
          df$node[df$node %in% dt$Disorder] %>% unique() %>% as.vector(),
          df$node[df$node %in% dt$IDP] %>% unique() %>% as.vector())
names(mycol2) = names
diff_col = setdiff(mycol1, mycol2)
same_node = intersect(names(mycol1), names(mycol2))
mycol2[match(same_node, names(mycol2))] = mycol1[match(same_node, names(mycol1))]
mycol2[setdiff(1:length(mycol2), match(same_node, names(mycol2)))] = sample(diff_col, length(mycol2) - length(same_node))


mycol2[which(names(mycol2) == 'Oligodendrocytes')] = '#386CB0'
mycol2[which(names(mycol2) == 'OPCs')] = '#4DC97F'
mycol2[which(names(mycol2) == 'Astrocytes')] = '#F0027F'
mycol2[which(names(mycol2) == 'Excitatory neurons')] = '#BEAED4'   
mycol2[which(names(mycol2) == 'Inhibitory neurons')] = '#FDC086'
mycol2[which(names(mycol2) == 'Microglia')] = '#FFFF99'

BRV = dt$IDP[!grepl('FA', dt$IDP)]
WMM = dt$IDP[grep('FA', dt$IDP)]
mycol2[which(names(mycol2) %in% BRV)] = '#4AB74C'
mycol2[which(names(mycol2) %in% WMM)] = '#0069E9'
BCP = c('CI', 'EA','INT','NEU','ASP','DPW','GRT','SC')
NEUD = c('INS','EPI','ICH','IS','AD','ALS','MS','PD')
PSYD = c('OCD','TS','ANX','MDD','ADHD','AUD','ASD','PTSD','BIP','SCZ')
mycol2[which(names(mycol2) %in% BCP)] = '#C58CA9'
mycol2[which(names(mycol2) %in% NEUD)] = '#E44D00'
mycol2[which(names(mycol2) %in% PSYD)] = '#FEB500'
mycol2[which(names(mycol2) %in% dt$gene_name)] = 'grey'


ggplot(df, aes(x = x,
               next_x = next_x,
               node = node,
               next_node = next_node,
               fill = node,
               label = node)) +
  scale_fill_manual(values = mycol2) + 
  geom_sankey(flow.alpha = 0.5,
              smooth = 8,
              width = 0.12) +
  geom_sankey_text(size = 5,
                   color = "black") +
  theme_void() +
  theme(legend.position = 'none')

