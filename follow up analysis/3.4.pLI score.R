rm(list=ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)
library(ggplot2)


#pli score
pli_score = fread('data/gnomad.v2.1.1.lof_metrics.by_gene.txt', header = T, data.table = F)
pli_score = pli_score[,c('gene', 'pLI')]


f1 = fread('results/signif_genes.CellToIDP.fdr.txt', data.table = F)
f2 = fread('results/signif_genes.CellToDisorder.fdr.txt', data.table = F)
f1$tissue_gene = paste(f1$tissue, f1$gene, sep='.')
f2$tissue_gene = paste(f2$tissue, f2$gene, sep='.')
trait_mapping = read.csv('gwas_meta_all.csv', header = T)
trait_mapping = trait_mapping[, c('TraitCategory1', 'TraitAbbr')]
PD = trait_mapping[trait_mapping$TraitCategory1 == 'Psychiatric disorder', 'TraitAbbr']
ND = trait_mapping[trait_mapping$TraitCategory1 == 'Neurological disorder', 'TraitAbbr']
BP = trait_mapping[trait_mapping$TraitCategory1 == 'Behavioral-cognitive phenotype', 'TraitAbbr']

set1 = f1[!grepl('FA', f1$pheno), ]  #Brain regional volume
score1 = subset(pli_score, pli_score$gene %in% set1$gene_name)
score1$group = 'Brain regional volume'
set2 = f1[grepl('FA', f1$pheno), ]  #White matter microstructure
score2 = subset(pli_score, pli_score$gene %in% set2$gene_name)
score2$group = 'White matter microstructure'
set3 = subset(f2, f2$pheno %in% PD)  #Psychiatric disorder
score3 = subset(pli_score, pli_score$gene %in% set3$gene_name)
score3$group = 'Psychiatric disorder'
set4 = subset(f2, f2$pheno %in% ND)  #Neurological disorder
score4 = subset(pli_score, pli_score$gene %in% set4$gene_name)
score4$group = 'Neurological disorder'
set5 = subset(f2, f2$pheno %in% BP)  #Behavioral-cognitive disorder
score5 = subset(pli_score, pli_score$gene %in% set5$gene_name)
score5$group = 'Behavioral-cognitive disorder'


p = list()
p <- list(#`set1` = set1$tissue_gene,
          #`set2` = set2$tissue_gene,
          #`set3` = set3$tissue_gene,
          #`set4` = set4$tissue_gene,
          #`set5` = set5$tissue_gene,
          `set1 & set1` = names(table(set1$tissue_gene))[which(table(set1$tissue_gene) > 1)], #shared in at least two traits in BRV
          `set2 & set2` = names(table(set2$tissue_gene))[which(table(set2$tissue_gene) > 1)],
          `set3 & set3` = names(table(set3$tissue_gene))[which(table(set3$tissue_gene) > 1)],
          `set4 & set4` = names(table(set4$tissue_gene))[which(table(set4$tissue_gene) > 1)],
          `set5 & set5` = names(table(set5$tissue_gene))[which(table(set5$tissue_gene) > 1)],
          
          `set1 & set2` = Reduce(intersect, list(set1$tissue_gene, set2$tissue_gene)), 
          `set1 & set3` = Reduce(intersect, list(set1$tissue_gene, set3$tissue_gene)),
          `set1 & set4` = Reduce(intersect, list(set1$tissue_gene, set4$tissue_gene)),
          `set1 & set5` = Reduce(intersect, list(set1$tissue_gene, set5$tissue_gene)),
          `set2 & set3` = Reduce(intersect, list(set2$tissue_gene, set3$tissue_gene)),
          `set2 & set4` = Reduce(intersect, list(set2$tissue_gene, set4$tissue_gene)),
          `set2 & set5` = Reduce(intersect, list(set2$tissue_gene, set5$tissue_gene)),
          `set3 & set4` = Reduce(intersect, list(set3$tissue_gene, set4$tissue_gene)),
          `set3 & set5` = Reduce(intersect, list(set3$tissue_gene, set5$tissue_gene)),
          `set4 & set5` = Reduce(intersect, list(set4$tissue_gene, set5$tissue_gene)),
          `set1 & set2 & set3` = Reduce(intersect, list(set1$tissue_gene,set2$tissue_gene,set3$tissue_gene)),
          `set1 & set2 & set4` = Reduce(intersect, list(set1$tissue_gene,set2$tissue_gene,set4$tissue_gene)),
          `set1 & set2 & set5` = Reduce(intersect, list(set1$tissue_gene,set2$tissue_gene,set5$tissue_gene)),
          `set2 & set3 & set4` = Reduce(intersect, list(set2$tissue_gene,set3$tissue_gene,set4$tissue_gene)),
          `set2 & set3 & set5` = Reduce(intersect, list(set2$tissue_gene,set3$tissue_gene,set5$tissue_gene)),
          `set3 & set4 & set5` = Reduce(intersect, list(set3$tissue_gene,set4$tissue_gene,set5$tissue_gene)),
          `set1 & set2 & set3 & set4` = Reduce(intersect, list(set1$tissue_gene,set2$tissue_gene,set3$tissue_gene,set4$tissue_gene)),
          `set1 & set2 & set3 & set5` = Reduce(intersect, list(set1$tissue_gene,set2$tissue_gene,set3$tissue_gene,set5$tissue_gene)),
          `set2 & set3 & set4 & set5` = Reduce(intersect, list(set2$tissue_gene,set3$tissue_gene,set4$tissue_gene,set5$tissue_gene)),
          `set1 & set2 & set3 & set4 & set5` = Reduce(intersect, list(set1$tissue_gene,set2$tissue_gene,set3$tissue_gene,set4$tissue_gene,set5$tissue_gene))
          )


gene_set1 = unique(pli_score[which(pli_score$pLI >=0.95), 'gene'])


mapping = rbind(f1[,c("gene_name",'tissue_gene')], f2[, c("gene_name",'tissue_gene')])
mapping = mapping[!(duplicated(mapping$tissue_gene)), ]

result = data.frame(stringsAsFactors = F)
for (i in 1:length(p)){
  p[[i]] = subset(mapping, mapping$tissue_gene %in% p[[i]]) $ gene_name

  gene_set2 = p[[i]][which(p[[i]] %in% gene_set1)]
  q = length(which(gene_set1 %in% gene_set2))
  m = length(gene_set1)
  n = nrow(pli_score) - m
  k = length(gene_set2)
  
  p_hyper = phyper(q-1, m, n, k, lower.tail=F)
  result.tmp = data.frame(group = names(p)[i], p = p_hyper)
  result = rbind(result, result.tmp)
}
result$group = gsub('set1', 'BRV', result$group) #Brain regional volume
result$group = gsub('set2', 'WMM',result$group)  #White matter microstructure
result$group = gsub('set3', 'PD', result$group)  #Psychiatric disorder
result$group = gsub('set4', 'ND', result$group)  #Neurological disorder
result$group = gsub('set5', 'BCP', result$group) #Behavioral-cognitive disorder

group_high = result[which(result$p < 0.05), 'group']
group_low = result[which(result$p >= 0.05), 'group'] 


for (i in 1:length(p)){
  p[[i]] = subset(pli_score, pli_score$gene %in% p[[i]]) $ pLI
}

data = data.frame(stringsAsFactors = F)
for (i in 1:length(p)){
  data.tmp = data.frame(pLI = p[[i]])
  if(nrow(data.tmp) == 0) next
  data.tmp$group = names(p)[i]
  data = rbind(data.tmp, data)
}

data$group = gsub('set1', 'BRV', data$group)   #Brain regional volume
data$group = gsub('set2', 'WMM',data$group)  #White matter microstructure
data$group = gsub('set3', 'PD', data$group)  #Psychiatric disorder
data$group = gsub('set4', 'ND', data$group)  #Neurological disorder
data$group = gsub('set5', 'BCP', data$group) #Behavioral-cognitive Phenotypes




plot = aggregate(data$pLI, by=list(data$group), mean)
names(plot)[2] = 'mean'
plot$var = aggregate(data$pLI, by=list(data$group), var)$x
plot = plot[-which(is.na(plot$var)), ]
look = aggregate(data$pLI, by=list(data$group), min)
look$max = aggregate(data$pLI, by=list(data$group), max)$x

plot = plot[order(plot$mean, decreasing = T), ]
plot$Group.1 = factor(plot$Group.1, levels = unique(plot$Group.1))


color = rep('#F09F96', nrow(plot))
ggplot(plot, aes(x = Group.1, y = mean)) + 
  geom_segment(aes(x = Group.1, xend = Group.1, y = mean-var, yend = mean+var), linewidth = 1.2, color = color) +
  geom_segment(aes(x = 1:nrow(plot)-0.1, xend = 1:nrow(plot)+0.1, y = mean-var, yend = mean-var), linewidth = 1.2, color = color) +
  geom_segment(aes(x = 1:nrow(plot)-0.1, xend = 1:nrow(plot)+0.1, y = mean+var, yend = mean+var), linewidth = 1.2, color = color) +
  geom_point(size = 5, color = color) + 
  scale_shape_manual(values = c(19,17)) + #指定点的性状
  labs(y = "pLI score") +
  theme(plot.title = element_text(hjust = 0.5, size = 7.5), #标题居中, 字体设置
        panel.background = element_rect(fill = 'white', color = 'black', linewidth = 1.4), #图片周围加框
        panel.grid = element_blank(), #不要网格线
        legend.position = 'none', ##不要图例
        axis.title.x = element_blank(),  #不要x轴标题和标签
        axis.title.y = element_text(size=16),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = 12,color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = 'black'))

###Between class
p = list()
plot1 = subset(plot, plot$Group.1 %in% c('BRV & ND', 'WMM & ND','BRV & PD',
                                         'BRV & BCP','WMM & PD', 'WMM & BCP'))  #pairwise

color = '#675E88'
p[[1]] = ggplot(plot1, aes(x = Group.1, y = mean)) + 
  geom_segment(aes(x = Group.1, xend = Group.1, y = mean-var, yend = mean+var), linewidth = 1.2, color = color) +
  geom_segment(aes(x = 1:nrow(plot1)-0.1, xend = 1:nrow(plot1)+0.1, y = mean-var, yend = mean-var), linewidth = 1.2, color = color) +
  geom_segment(aes(x = 1:nrow(plot1)-0.1, xend = 1:nrow(plot1)+0.1, y = mean+var, yend = mean+var), linewidth = 1.2, color = color) +
  geom_point(size = 5, color = color) + 
  scale_shape_manual(values = c(19,17)) +
  labs(y = "pLI score") +
  scale_y_continuous(limits = c(0, 0.7), breaks = seq(0, 0.7, by = 0.1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 7.5), 
        panel.background = element_rect(fill = 'white', color = 'black', linewidth = 1.4),
        panel.grid = element_blank(),
        legend.position = 'none', 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = 12,color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = 'black'))


###Within class
plot2 = subset(plot, plot$Group.1 %in% c('PD & ND','ND & BCP', 'PD & BCP', 'BRV & WMM'))  #pairwise
color = '#F09F96'
p[[2]] = ggplot(plot2, aes(x = Group.1, y = mean)) + 
  geom_segment(aes(x = Group.1, xend = Group.1, y = mean-var, yend = mean+var), linewidth = 1.2, color = color) +
  geom_segment(aes(x = 1:nrow(plot2)-0.1, xend = 1:nrow(plot2)+0.1, y = mean-var, yend = mean-var), linewidth = 1.2, color = color) +
  geom_segment(aes(x = 1:nrow(plot2)-0.1, xend = 1:nrow(plot2)+0.1, y = mean+var, yend = mean+var), linewidth = 1.2, color = color) +
  geom_point(size = 5, color = color) + 
  scale_shape_manual(values = c(19,17)) +
  labs(y = "pLI score") +
  scale_y_continuous(limits = c(0, 0.7), breaks = seq(0, 0.7, by = 0.1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 7.5),
        panel.background = element_rect(fill = 'white', color = 'black', linewidth = 1.4),
        panel.grid = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=16),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = 12,color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = 'black'))


#within group
plot3 = subset(plot, plot$Group.1 %in% c( "BRV & BRV","WMM & WMM","PD & PD","BCP & BCP")) 
color = c('#C58CA9', '#4AB74C', '#0069E9', '#FEB500')
p[[3]] = ggplot(plot3, aes(x = Group.1, y = mean)) + 
  geom_segment(aes(x = Group.1, xend = Group.1, y = mean-var, yend = mean+var), linewidth = 1.2, color = color) +
  geom_segment(aes(x = 1:nrow(plot3)-0.1, xend = 1:nrow(plot3)+0.1, y = mean-var, yend = mean-var), linewidth = 1.2, color = color) +
  geom_segment(aes(x = 1:nrow(plot3)-0.1, xend = 1:nrow(plot3)+0.1, y = mean+var, yend = mean+var), linewidth = 1.2, color = color) +
  geom_point(size = 5, color = color) + 
  scale_shape_manual(values = c(19,17)) +
  labs(y = "pLI score") +
  scale_y_continuous(limits = c(0, 0.7), breaks = seq(0, 0.7, by = 0.1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 7.5),
        panel.background = element_rect(fill = 'white', color = 'black', linewidth = 1.4),
        panel.grid = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = 12,color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = 'black'))

library(gridExtra)
text = 'grid.arrange(p[[3]], p[[2]], p[[1]], ncol=3, nrow=1)'
eval(parse(text = text))
