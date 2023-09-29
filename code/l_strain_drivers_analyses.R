Sys.setenv("LANGUAGE"="EN")
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
  model_out = read.csv('stats/model_res_all_scenarios_final.csv') %>% filter(scenario == 370)
  min_non_zero = min(model_out$median_future[model_out$median_future > 0])
  model_out$wilcox_p[is.na(model_out$wilcox_p)] = 1 # NAs for models that selected only variables that do not change (0 median change, etc.)
  model_out$median_future[model_out$median_future < 0] = min_non_zero/2
  model_out$median_present[model_out$median_present < 0] = min_non_zero/2
  model_out$log2fc = log2(model_out$median_future) - log2(model_out$median_present)

  model_out$Change = 'Others'
  model_out$Change[(model_out$wilcox_p < 0.05) & (model_out$median_change < 0)] = 'Decrease'

  importance_tab = model_out %>% filter(r2 >= 0.05) %>% group_by(MAG, Change) %>% mutate(ImpRank=rank(Freq)/length(Freq), log2fc=median(log2fc)) %>% 
                                 select(MAG, vars_selection_tab, ImpRank, mean_rel_ab, r2, Change, median_change)
  change_tab = model_out %>% filter(r2 >= 0.05) %>% group_by(MAG, Change) %>% select(MAG, mean_rel_ab, r2, Change, median_change) %>% distinct()
  
  tax_changes_vars_stats = data.frame()
  
  pred_krusk = kruskal.test(change_tab, r2 ~ Class)
  chng_krusk = kruskal.test(change_tab, log2fc ~ Class)
  
  tax_changes_vars_stats = rbind(tax_changes_vars_stats, data.frame(Variable='Predictability', stat=pred_krusk$statistic, p=pred_krusk$p.value))
  tax_changes_vars_stats = rbind(tax_changes_vars_stats, data.frame(Variable='Predicted_change', stat=chng_krusk$statistic, p=chng_krusk$p.value))
  
  for (var in unique(importance_tab$vars_selection_tab)){
    krusk = kruskal.test(importance_tab %>% filter(vars_selection_tab == var), ImpRank ~ Class)
    kstat = krusk$statistic
    kp = krusk$p.value
    tax_changes_vars_stats = rbind(tax_changes_vars_stats, data.frame(Variable=var, stat=kstat, p=kp))}
  write.csv(tax_changes_vars_stats, 'stats/tax_kruskal_stats.csv', quote = F, row.names = F)
  

  importance_tab$Category = 'Others'
  importance_tab$Category[startsWith(importance_tab$vars_selection_tab, 'bioclim_')] = 'Bioclimatic'
  importance_tab$Category[startsWith(importance_tab$vars_selection_tab, 'min_')] = 'Minerals'
  importance_tab$Category[startsWith(importance_tab$vars_selection_tab, 'pc_')] = 'Phy.-Chem.'
  importance_tab$Category[startsWith(importance_tab$vars_selection_tab, 'clim_')] = 'Climatic'
  importance_tab$Category[startsWith(importance_tab$vars_selection_tab, 'gl_')] = 'Glacio.'
  print(unique(importance_tab$vars_selection_tab))
  
  taxonomy = read.csv('data/raw/bacteria/gtdbtk.bac120.summary.tsv', sep='\t')
  importance_tab = importance_tab[importance_tab$MAG %in% taxonomy$user_genome,]
  importance_tab$Taxonomy = map_chr(importance_tab$MAG, function(x) taxonomy$classification[taxonomy$user_genome == x])
  importance_tab$Phylum = map_chr(importance_tab$Taxonomy, function(x) strsplit(x, split = ';')[[1]][2])
  importance_tab$Class = map_chr(importance_tab$Taxonomy, function(x) gsub('c__', '', strsplit(x, split = ';')[[1]][3]))
  return(list(imp=importance_tab, mod=model_out))}

LoadTree <- function(){
  tree = read.tree('data/raw/bacteria/treeBacteria.tree')
  tree = midpoint.root(tree)
  tree$tip.label = vapply(tree$tip.label, function(x) strsplit(x, '.fa')[[1]][1], FUN.VALUE = character(1))
  return(tree)}

PlotDrivers <- function(importance_tab, model_out){
print(colnames(importance_tab))
importance_tab_1 = importance_tab
importance_tab_1$Change = 'All'
sp1 = ggplot(importance_tab_1, aes(x=ImpRank*100, y=reorder(vars_selection_tab, ImpRank), fill=after_stat(x))) + 
  geom_density_ridges_gradient(scale = 1.5, rel_min_height = 0.01, quantile_lines=TRUE, quantiles=2) + 
  facet_grid(Category~Change, scales = 'free', space = 'free') + theme_linedraw() + ylab('') + xlab('Relative rank [%]') + 
  scale_y_discrete(labels=as_labeller(clean_names)) + theme(legend.position = 'none', panel.grid = element_blank(), axis.text=element_text(size=12)) + 
  scale_fill_gradientn(colours = c('#E04D55', '#15356B')) +  xlim(0,100) +
  theme(strip.text.y = element_blank()) + scale_x_continuous(breaks=c(0, 50, 100)) +
  theme(strip.background.y = element_blank())

importance_tab_2 = importance_tab
importance_tab_2$variable_all_median = map_dbl(importance_tab_2$vars_selection_tab, function(x) median(importance_tab_1$ImpRank[importance_tab_1$vars_selection_tab == x]*100))
importance_tab_2$variable_all_median_diff = importance_tab_2$variable_all_median - (importance_tab_2$ImpRank * 100)
sp2 = ggplot(importance_tab_2, aes(x=variable_all_median_diff, y=reorder(vars_selection_tab, ImpRank), fill=after_stat(x))) + geom_vline(xintercept = 0, colour='dimgrey') + 
  geom_density_ridges_gradient(scale = 1.5, rel_min_height = 0.01, quantile_lines=TRUE, quantiles=2) + 
  facet_grid(Category~Change, scales = 'free', space = 'free') + theme_linedraw() + ylab('') + xlab('Difference to All median [%]') + 
  scale_y_discrete(labels=as_labeller(clean_names)) + theme(legend.position = 'none', panel.grid = element_blank(), axis.text=element_text(size=12),axis.text.y = element_blank()) + # 
  xlim(-100,100) + scale_x_continuous(breaks=c(-50, 0, 50)) + 
  scale_fill_gradientn(colours = c('#E04D55', 'lightgrey', '#15356B')) 

p = ggarrange(sp1, sp2, nrow = 1, ncol = 2, labels = c('A', 'B'), widths = c(0.5,0.5))
ggsave(p, filename = 'plots/Fig_5_Drivers.pdf', width = 9, height = 7)}

PlotDriversClass <- function(importance_tab, model_out){
  class_table = table(importance_tab$Class)
  class_to_keep = names(class_table)[class_table > 333]
  importance_tab$Class[!(importance_tab$Class %in% class_to_keep)] = 'Others'

  importance_tab$Class = factor(importance_tab$Class, levels=c("Gammaproteobacteria","Bacteroidia","Alphaproteobacteria",
                                                               "Planctomycetia","Verrucomicrobiae","Acidimicrobiia",
                                                               "Actinomycetia","Gemmatimonadetes","Acidobacteriae",
                                                               "Polyangia","Bdellovibrionia","Others"))
  p = ggplot(importance_tab, aes(x=ImpRank*100, y=reorder(vars_selection_tab, ImpRank), fill=after_stat(x))) + 
    geom_density_ridges_gradient(scale = 1.5, rel_min_height = 0.01, quantile_lines=TRUE, quantiles=2) + 
    facet_grid(Category~Class, scales = 'free', space = 'free') + theme_linedraw() + ylab('') + xlab('Relative rank [%]') + 
    scale_y_discrete(labels=as_labeller(clean_names)) + theme(legend.position = 'none', panel.grid = element_blank(), axis.text=element_text(size=12)) + 
    scale_fill_gradientn(colours = c('#E04D55', '#15356B')) +  xlim(0,100) +
    theme(strip.text.y = element_blank()) + scale_x_continuous(breaks=c(0, 50, 100)) +
    theme(strip.background.y = element_blank())
  
  ggsave(p, filename = 'plots/Fig_S10_taxa_drivers.pdf', width = 20, height = 7)}


PlotPredClass <- function(importance_tab, model_out){
  class_table = table(importance_tab$Class)
  class_to_keep = names(class_table)[class_table > 333]
  importance_tab$Class[!(importance_tab$Class %in% class_to_keep)] = 'Others'
  
  importance_tab$Class = factor(importance_tab$Class, levels=c("Gammaproteobacteria","Bacteroidia","Alphaproteobacteria",
                                                               "Planctomycetia","Verrucomicrobiae","Acidimicrobiia",
                                                               "Actinomycetia","Gemmatimonadetes","Acidobacteriae",
                                                               "Polyangia","Bdellovibrionia","Others"))
  p = ggplot(importance_tab, aes(x=r2, y=Class, fill=after_stat(x))) + 
    geom_density_ridges_gradient(scale = 1.5, rel_min_height = 0.01, quantile_lines=TRUE, quantiles=2) + 
    theme_linedraw() + ylab('') + xlab(expression('r'^2*''[prediction])) + 
    theme(legend.position = 'none', panel.grid = element_blank(), axis.text=element_text(size=12)) + 
    scale_fill_gradientn(colours = c('#E04D55', '#15356B')) +  xlim(0,1) +
    theme(strip.text.y = element_blank()) + scale_x_continuous(breaks=c(0, 0.1, 0.25, 0.5, 1), trans='log10') +
    theme(strip.background.y = element_blank())
  
  ggsave(p, filename = 'plots/Fig_S8_taxa_predictability.pdf', width = 7, height = 7)}

PlotPhyloDrivers <- function(model_out, tree){
  taxonomy = read.csv('data/raw/bacteria/gtdbtk.bac120.summary.tsv', sep='\t')
  taxonomy$Genus = vapply(taxonomy$classification, function(x) strsplit( strsplit(x, ';s__')[[1]][1], ';g__')[[1]][2], FUN.VALUE = character(1))

  model_out = model_out %>% group_by(MAG) %>% arrange(desc(Freq)) %>% slice(1:5)
  #features_tab = data.frame(MAG=unique(model_out$MAG))
  #for (feat in unique(model_out$vars_selection_tab)){
  #  features_tab[feat] = 0
  #  for (i in 1:nrow(features_tab)){value = model_out$ImpRank[(model_out$vars_selection_tab == feat) & (model_out$MAG == features_tab$MAG[i])]
  #    if (length(value) > 0){features_tab[i,feat] = value}
  #    else {features_tab[i,feat] = 0}}}
  #rownames(features_tab) = features_tab$MAG
  #features_tab$MAG = NULL

  registerDoMC(9)
  dists = melt(as.matrix(cophenetic.phylo(tree)), value.name = 'PD') %>% filter(Var1 != Var2) %>% mutate(weights = (1 / PD)^2)
  dists = left_join(dists, taxonomy %>% select(user_genome, Genus), by = c('Var1' = 'user_genome'))
  dists = left_join(dists, taxonomy %>% select(user_genome, Genus), by = c('Var2' = 'user_genome'))
  median_same_genus = median(dists$PD[(dists$Genus.x == dists$Genus.y) & (!is.na(dists$Genus.x))], na.rm = T)

  dists_sub = dists %>% slice_sample(n = 100000, weight_by = weights)
  dists_sub$N_shared = foreach(i = 1:nrow(dists_sub), .combine = c) %dopar% {vars1 = model_out %>% filter(MAG == dists_sub$Var1[i]) %>% pull(vars_selection_tab)
                                                                             vars2 = model_out %>% filter(MAG == dists_sub$Var2[i]) %>% pull(vars_selection_tab)
                                                                             return(length(intersect(vars1,vars2)))}

  p = ggplot(dists_sub, aes(x=PD, y=N_shared)) + geom_point(alpha=0.2) + geom_vline(xintercept = median_same_genus, colour='red') + 
    geom_smooth(method='loess', se=F) + theme_bw() + scale_x_log10() + xlab('Phylogenetic distance') + ylab('Number of shared top 5 variables') + stat_cor(method='spearman')
  ggsave(p, filename = 'plots/Fig_S9_Distance_variable_ranks.pdf')
  print(cor.test(dists_sub$N_shared, dists_sub$PD, method='spearman'))}

MainL <- function(){
  data = LoadData()
  importance_tab = data$imp
  model_out = data$mod
  
  PlotDrivers(importance_tab, model_out)
  PlotDriversClass(importance_tab, model_out)
  PlotPredClass(importance_tab, model_out)
  
  tree = LoadTree()
  changmodel_outes_tab = model_out[model_out$MAG %in% tree$tip.label,]
  tree = keep.tip(tree, tree$tip.label[tree$tip.label %in% model_out$MAG])

  PlotPhyloDrivers(model_out, tree)
}




























