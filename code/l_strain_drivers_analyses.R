library(mgcv)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(reshape2)
library(matrixStats)
library(gtools)
library(phytools)
library(ape)
library(ggtree)
library(ggtreeExtra)
library(forcats)
library(viridis)
library(ggnewscale)
library(doMC)
library(foreach)
library(ggridges)

#install.packages('rankdist',repos = "http://cran.us.r-project.org")
library(rankdist)

clean_names = c('bioclim_PC6'='PC6', 'bioclim_PC5'='PC5', 'bioclim_PC4'='PC4', 'bioclim_PC3'='PC3', 'bioclim_PC2'='PC2', 'bioclim_PC1'='PC1',
                'clim_tas'='Monthly temperature', 'clim_scd'='Annual snow cover', 'clim_pr'='Monthly precipitation', 'clim_tasmin'='Min. daily temperature', 'clim_tasmax'='Max. daily temperature',
                'gl_distance'='Distance to the glacier', 'gl_coverage'='Glacier coverage', 'gl_area'='Glacier area', 'slope'='Slope', 
                'abs_latitude'='Absolute latitude', 'latitude'='Latitude', 'longitude'='Longitude', 'elevation'='Elevation',
                'min_quartz'='Quartz', 'min_feldspar'='Feldspar', 'min_clays'='Clays', 'min_calcite'='Calcite',
                'nut_srp_predicted'='Soluble reactive phosphate', 'nut_din_predicted'='Dissolved inorganic nitrogen', 'chla_predicted'='Chlorophyll-a',
                'pc_water_temp_predicted'='Water temperature', 'pc_turbidity_predicted'='Turbidity', 'pc_ph_predicted'='pH', 'pc_conductivity_predicted'='Conductivity')

LoadData <- function(){
  model_out = read.csv('stats/model_res_all_scenarios_final.csv')
  min_non_zero = min(model_out$median_future[model_out$median_future > 0])
  model_out$wilcox_p[is.na(model_out$wilcox_p)] = 1 # NAs for models that selected only variables that do not change (0 median change, etc.)
  model_out$median_future[model_out$median_future < 0] = min_non_zero/2
  model_out$median_present[model_out$median_present < 0] = min_non_zero/2
  model_out$log2fc = log2(model_out$median_future) - log2(model_out$median_present)

  model_out$Category = 'Not significant'
  model_out$Category[(model_out$wilcox_p < 0.05) & (model_out$median_change > 0)] = 'Increase'
  model_out$Category[(model_out$wilcox_p < 0.05) & (model_out$median_change < 0)] = 'Decrease'

  importance_tab = model_out %>% filter(r2 >= 0.05) %>% group_by(MAG) %>% mutate(ImpRank=rank(Freq)/length(Freq)) %>% select(MAG, vars_selection_tab, ImpRank, mean_rel_ab)

  importance_tab$Category = 'Others'
  importance_tab$Category[startsWith(importance_tab$vars_selection_tab, 'bioclim_')] = 'Bioclimatic'
  importance_tab$Category[startsWith(importance_tab$vars_selection_tab, 'min_')] = 'Minerals'
  importance_tab$Category[startsWith(importance_tab$vars_selection_tab, 'pc_')] = 'Phy.-Chem.'
  importance_tab$Category[startsWith(importance_tab$vars_selection_tab, 'clim_')] = 'Climatic'
  importance_tab$Category[startsWith(importance_tab$vars_selection_tab, 'gl_')] = 'Glacio.'
  print(unique(importance_tab$vars_selection_tab))

  return(list(imp=importance_tab, mod=model_out))}

LoadTree <- function(){
  tree = read.tree('data/raw/bacteria/treeBacteria.tree')
  tree = midpoint.root(tree)
  tree$tip.label = vapply(tree$tip.label, function(x) strsplit(x, '.fa')[[1]][1], FUN.VALUE = character(1))
  return(tree)}

PlotDrivers <- function(importance_tab, model_out){
sp1 = importance_tab %>% group_by(vars_selection_tab, Category) %>% summarise(Importance = weightedMedian(ImpRank, w = mean_rel_ab)) %>% 
  ggplot(aes(x=Importance*100, y=reorder(vars_selection_tab, Importance), fill=Importance)) + geom_bar(stat='identity') + 
  facet_grid(Category~., scales = 'free', space = 'free') + theme_linedraw() + ylab('') + xlab('Median rel. rank [%]') + 
  scale_y_discrete(labels=as_labeller(clean_names)) + theme(legend.position = 'none', panel.grid = element_blank(), axis.text=element_text(size=12)) + 
  scale_fill_gradientn(colours = c('#E04D55', '#15356B'))

# Difference between decreasing / increasing MAGs
model_out$Category = 'Not significant and increase'
model_out$Category[(model_out$wilcox_p < 0.05) & (model_out$median_change < 0)] = 'Decrease'

vars_comp = model_out %>% group_by(MAG) %>% mutate(ImpRank=rank(Freq)/length(Freq)) %>% 
  select(MAG, vars_selection_tab, ImpRank, r2, Category) %>% group_by(vars_selection_tab) %>% 
  do(w = wilcox.test(ImpRank~Category, data=., paired=FALSE)$p.value) 

vars_comp$median_dec = map_dbl(vars_comp$vars_selection_tab, function(x) model_out %>% filter(Category == 'Decrease') %>% group_by(MAG) %>% 
                                 mutate(ImpRank=rank(Freq)/length(Freq)) %>% ungroup() %>% filter(vars_selection_tab == x) %>% pull(ImpRank) %>% median())
vars_comp$median_others = map_dbl(vars_comp$vars_selection_tab, function(x) model_out %>% filter(Category != 'Decrease') %>% group_by(MAG) %>% 
                                 mutate(ImpRank=rank(Freq)/length(Freq)) %>% ungroup() %>% filter(vars_selection_tab == x) %>% pull(ImpRank) %>% median())
vars_comp$median_diff = vars_comp$median_dec - vars_comp$median_others
vars_comp$Significant = 'No'
vars_comp$Significant[vars_comp$w < 0.05] = 'Yes'

sp2 = vars_comp %>% ggplot(aes(x=median_diff*100, y=reorder(vars_selection_tab, median_diff), fill=Significant)) + 
  geom_bar(stat='identity') + theme_linedraw() + ylab('') + xlab('Median diff. in rank [%]') + scale_y_discrete(labels=as_labeller(clean_names)) + 
  theme(legend.position = 'none', panel.grid = element_blank(), axis.text=element_text(size=12)) + scale_fill_manual(values = c('#E04D55', '#15356B'))

p = ggarrange(sp1, sp2, nrow = 1, ncol = 2, labels = c('A', 'B'), widths = c(0.5,0.5))
ggsave(p, filename = 'plots/Fig_5_Drivers.pdf', width = 9, height = 7)}

PlotPhyloDrivers <- function(model_out, tree){
  taxonomy = read.csv('data/raw/bacteria/gtdbtk.bac120.summary.tsv', sep='\t')
  taxonomy$Genus = vapply(taxonomy$classification, function(x) strsplit( strsplit(x, ';s__')[[1]][1], ';g__')[[1]][2], FUN.VALUE = character(1))

  model_out = model_out %>% filter(scenario == 370) %>% group_by(MAG) %>% mutate(ImpRank=rank(Freq)/length(Freq)) %>% ungroup()
  features_tab = data.frame(MAG=unique(model_out$MAG))
  for (feat in unique(model_out$vars_selection_tab)){
    features_tab[feat] = 0
    for (i in 1:nrow(features_tab)){value = model_out$ImpRank[(model_out$vars_selection_tab == feat) & (model_out$MAG == features_tab$MAG[i])]
      if (length(value) > 0){features_tab[i,feat] = value}
      else {features_tab[i,feat] = 0}}}
  rownames(features_tab) = features_tab$MAG
  features_tab$MAG = NULL

  registerDoMC(9)
  dists = melt(as.matrix(cophenetic.phylo(tree)), value.name = 'PD') %>% filter(Var1 != Var2) %>% mutate(weights = (1 / PD)^2)
  dists = left_join(dists, taxonomy %>% select(user_genome, Genus), by = c('Var1' = 'user_genome'))
  dists = left_join(dists, taxonomy %>% select(user_genome, Genus), by = c('Var2' = 'user_genome'))
  median_same_genus = median(dists$PD[(dists$Genus.x == dists$Genus.y) & (!is.na(dists$Genus.x))], na.rm = T)

  dists_sub = dists %>% slice_sample(n = 10000, weight_by = weights)
  dists_sub$Distance = foreach(i = 1:nrow(dists_sub), .combine = c) %dopar% {rank1 = rank(features_tab[rownames(features_tab) == dists_sub$Var1[i],])
                                                                            rank2 = rank(features_tab[rownames(features_tab) == dists_sub$Var2[i],])
                                                                            return(DistancePair(rank1, rank2))}

  p = ggplot(dists_sub, aes(x=PD, y=Distance)) + geom_point(alpha=0.2) + geom_vline(xintercept = median_same_genus, colour='red') + 
    geom_smooth(method='loess', se=F) + theme_bw() + scale_x_log10() + xlab('Phylogenetic distance') + ylab('Kendall distance of variable ranks') + stat_cor(method='spearman')
  ggsave(p, filename = 'plots/Fig_S8_Distance_variable_ranks.pdf')
  print(cor.test(dists_sub$Distance, dists_sub$PD, method='spearman'))}

MainL <- function(){
  data = LoadData()
  importance_tab = data$imp
  model_out = data$mod
  PlotDrivers(importance_tab, model_out)

  tree = LoadTree()
  changmodel_outes_tab = model_out[model_out$MAG %in% tree$tip.label,]
  tree = keep.tip(tree, tree$tip.label[tree$tip.label %in% model_out$MAG])

  PlotPhyloDrivers(model_out, tree)}




























