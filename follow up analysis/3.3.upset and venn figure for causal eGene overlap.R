rm(list=ls())
setwd('F:/类脑/1_实操文件/genexpMR/')
library(data.table)
library(stats)


f1 = fread('results/signif_genes.CellToIDP.fdr.txt', data.table = F)
f2 = fread('results/signif_genes.CellToDisorder.fdr.txt', data.table = F)
trait_mapping = read.csv('gwas_meta_all.csv', header = T) 
trait_mapping = trait_mapping[, c('TraitCategory1', 'TraitAbbr')]
PD = trait_mapping[trait_mapping$TraitCategory1 == 'Psychiatric disorder', 'TraitAbbr']
ND = trait_mapping[trait_mapping$TraitCategory1 == 'Neurological disorder', 'TraitAbbr']
BP = trait_mapping[trait_mapping$TraitCategory1 == 'Behavioral-cognitive phenotype', 'TraitAbbr']

set1 = f1[!grepl('FA', f1$pheno), ]  #Brain regional volume
set1$tissue_gene = paste(set1$tissue, set1$gene, sep='.')
set1 = unique(set1$tissue_gene)
set2 = f1[grepl('FA', f1$pheno), ]  #White matter microstructure
set2$tissue_gene = paste(set2$tissue, set2$gene, sep='.')
set2 = unique(set2$tissue_gene)
set3 = subset(f2, f2$pheno %in% PD)  #Psychiatric disorder
set3$tissue_gene = paste(set3$tissue, set3$gene, sep='.')
set3 = unique(set3$tissue_gene)
set4 = subset(f2, f2$pheno %in% ND)  #Neurological disorder
set4$tissue_gene = paste(set4$tissue, set4$gene, sep='.')
set4 = unique(set4$tissue_gene)
set5 = subset(f2, f2$pheno %in% BP)  #Behavioral-cognitive disorder
set5$tissue_gene = paste(set5$tissue, set5$gene, sep='.')
set5 = unique(set5$tissue_gene)


#Upset figure--------------------------------------------------------------------
library(UpSetR)  
library(RColorBrewer)
listInput = list(set1 = set1, set2 = set2, set3 = set3, set4 = set4,  set5 = set5)

expressionInput <- c(
  `set1&set2` = length(Reduce(intersect, list(set1, set2))), 
  `set1&set3` = length(Reduce(intersect, list(set1, set3))),
  `set1&set4` = length(Reduce(intersect, list(set1, set4))),
  `set1&set5` = length(Reduce(intersect, list(set1, set5))),
  `set2&set3` = length(Reduce(intersect, list(set2, set3))),
  `set2&set4` = length(Reduce(intersect, list(set2, set4))),
  `set2&set5` = length(Reduce(intersect, list(set2, set4))),
  `set3&set4` = length(Reduce(intersect, list(set3, set4))),
  `set3&set5` = length(Reduce(intersect, list(set3, set4))),
  `set4&set5` = length(Reduce(intersect, list(set3, set4))),
  `set1&set2&set3` = length(Reduce(intersect, list(set1,set2,set3))),
  `set1&set2&set4` = length(Reduce(intersect, list(set1,set2,set4))),
  `set1&set2&set5` = length(Reduce(intersect, list(set1,set2,set5))),
  `set2&set3&set4` = length(Reduce(intersect, list(set2,set3,set4))),
  `set2&set3&set5` = length(Reduce(intersect, list(set2,set3,set5))),
  `set3&set4&set5` = length(Reduce(intersect, list(set3,set4,set5))),
  `set1&set2&set3&set4` = length(Reduce(intersect, list(set1,set2,set3,set4))),
  `set1&set2&set3&set5` = length(Reduce(intersect, list(set1,set2,set3,set5))),
  `set2&set3&set4&set5` = length(Reduce(intersect, list(set2,set3,set4,set5))),
  `set1&set2&set3&set4&set5` = length(Reduce(intersect, list(set1,set2,set3,set4,set5)))
  )

names(expressionInput) = gsub('set1', 'Brain regional volume', names(expressionInput))
names(expressionInput) = gsub('set2', 'White matter microstructure', names(expressionInput))
names(expressionInput) = gsub('set3', 'Psychiatric disorder', names(expressionInput))
names(expressionInput) = gsub('set4', 'Neurological disorder', names(expressionInput))
names(expressionInput) = gsub('set5', 'Behavioral-cognitive disorder', names(expressionInput))


upset(fromExpression(expressionInput), order.by = "freq",
      matrix.color = 'black', main.bar.color = '#675E88', 
      sets.bar.color = c('#FEB90E', '#C58CA9','#79B77A','#0069E9','#E46C2F'),#brewer.pal(5, "Set1"), #'#BEBADA',
      text.scale = 1.8, point.size = 3, 
      shade.col = 'grey', shade.alpha = 0.2, 
      sets.x.label = 'NsignifGenes',
      queries = list(
        list(query=intersects, params=list("Brain regional volume", "White matter microstructure"), color="#F09F96", active=T),
        list(query=intersects, params=list("Psychiatric disorder", "Neurological disorder"), color="#F09F96", active=T),
        list(query=intersects, params=list("Neurological disorder", "Behavioral-cognitive disorder"), color="#F09F96", active=T),
        list(query=intersects, params=list("Psychiatric disorder", "Behavioral-cognitive disorder"), color="#F09F96", active=T),
        list(query=intersects, params=list("Psychiatric disorder","Neurological disorder","Behavioral-cognitive disorder"), color="#F09F96", active=T),
        
        list(query=intersects, params=list("Brain regional volume", "Psychiatric disorder"), color="#675E88", active=T),
        list(query=intersects, params=list("Brain regional volume", "Neurological disorder"), color="#675E88", active=T),
        list(query=intersects, params=list("Brain regional volume", "Behavioral-cognitive disorder"), color="#675E88", active=T),
        list(query=intersects, params=list("White matter microstructure", "Psychiatric disorder"), color="#675E88", active=T),
        list(query=intersects, params=list("White matter microstructure", "Neurological disorder"), color="#675E88", active=T),
        list(query=intersects, params=list("White matter microstructure", "Behavioral-cognitive disorder"), color="#675E88", active=T),
        list(query=intersects, params=list("White matter microstructure", "Brain regional volume", "Psychiatric disorder"), color="#675E88", active=T),
        list(query=intersects, params=list("White matter microstructure", "Brain regional volume","Neurological disorder"), color="#675E88", active=T),
        list(query=intersects, params=list("White matter microstructure", "Brain regional volume","Behavioral-cognitive disorder"), color="#675E88", active=T),
        list(query=intersects, params=list("White matter microstructure","Psychiatric disorder","Behavioral-cognitive disorder"), color="#675E88", active=T),
        list(query=intersects, params=list("White matter microstructure", "Brain regional volume","Psychiatric disorder","Behavioral-cognitive disorder"), color="#675E88", active=T),
        list(query=intersects, params=list("White matter microstructure", "Psychiatric disorder","Neurological disorder"), color="#675E88", active=T),
        list(query=intersects, params=list("White matter microstructure", "Brain regional volume","Psychiatric disorder","Neurological disorder"), color="#675E88", active=T),
        list(query=intersects, params=list("White matter microstructure", "Psychiatric disorder","Neurological disorder","Behavioral-cognitive disorder"), color="#675E88", active=T),
        list(query=intersects, params=list("White matter microstructure", "Brain regional volume","Psychiatric disorder","Neurological disorder","Behavioral-cognitive disorder"), color="#675E88", active=T))
      )



#venn figure ---------------------------------------------------------
#'set1': Brain regional volume
#'set2': White matter microstructure
#'set3': Psychiatric disorder
#'set4': Neurological disorder
#'set5': Behavioral-cognitive phenotype
set_pair = data.frame(set_i = c('set1','set1','set1','set5','set3','set3','set1','set2','set2','set2'),
                      set_j = c('set2','set3','set5','set4','set4','set5','set4','set4','set5','set3'))
mapping = c(set1 = 'Brain regional volume',
            set2 = 'White matter microstructure',
            set3 = 'Psychiatric disorder',
            set4 = 'Neurological disorder',
            set5 = 'Behavioral-cognitive phenotype')
library(ggvenn)
p = list()
count = 0
for (n_row in 1:nrow(set_pair)){
  count = count + 1
  set_i = set_pair[n_row, 'set_i']
  set_j = set_pair[n_row, 'set_j']
  name_i = as.character(mapping[set_i])
  name_j = as.character(mapping[set_j])

  q = length(which(eval(parse(text = set_i)) %in% eval(parse(text = set_j))) )
  m = length(eval(parse(text = set_i)))
  n = 26861*8 - m
  k = length(eval(parse(text = set_j)))
  p_phyer = phyper(q-1, m, n, k, lower.tail=F)
  
  
  a = list(set_i = eval(parse(text = set_i)), set_j =  eval(parse(text = set_j)))
  names(a) = c(name_i, name_j)
  
  p[[count]] = ggvenn(a, names(a), fill_color = c("pink", "grey"), 
                      show_percentage = F, text_size = 4,
                      stroke_size = 1, set_name_size = 2.5) +
    ggtitle(paste('p_phyer', p_phyer, sep=' = ')) +
    theme(plot.title = element_text(size = 8, hjust = 0.5))
}

library(gridExtra)
text = 'grid.arrange(p[[1'
for (i in 2:10){
  text = paste(text, ']], p[[', i, sep='')
}	
text = paste(text, ']], ncol=5, nrow=2)', sep='')
eval(parse(text = text))

