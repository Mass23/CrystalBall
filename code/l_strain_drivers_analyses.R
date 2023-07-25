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
library(lsa)
library(ggridges)
library(rankdist)

setwd('~/Documents/PhD/SpeciesDistributions')

# 2842 MAGs models
model_out = read.csv('data/processed/model_res_all_scenarios_final.csv')
min_non_zero = min(model_out$median_future[model_out$median_future > 0])
model_out$median_future[model_out$median_future < 0] = min_non_zero/2
model_out$median_present[model_out$median_present < 0] = min_non_zero/2
model_out$log2fc = log2(model_out$median_future) - log2(model_out$median_present)

taxonomy = read.csv('data/raw/bacteria/gtdbtk.bac120.summary.tsv', sep='\t')
taxonomy$Phylum = vapply(taxonomy$classification, function(x) strsplit( strsplit(x, ';c__')[[1]][1], ';p__')[[1]][2], FUN.VALUE = character(1))
taxonomy$Class = vapply(taxonomy$classification, function(x) strsplit( strsplit(x, ';o__')[[1]][1], ';c__')[[1]][2], FUN.VALUE = character(1))
taxonomy$Order = vapply(taxonomy$classification, function(x) strsplit( strsplit(x, ';f__')[[1]][1], ';o__')[[1]][2], FUN.VALUE = character(1))
taxonomy$Family = vapply(taxonomy$classification, function(x) strsplit( strsplit(x, ';g__')[[1]][1], ';f__')[[1]][2], FUN.VALUE = character(1))
taxonomy$Genus = vapply(taxonomy$classification, function(x) strsplit( strsplit(x, ';s__')[[1]][1], ';g__')[[1]][2], FUN.VALUE = character(1))

########################################################################################################################################
changes_tab = model_out %>% select(-vars_selection_tab, -Freq, scenario) %>% distinct()
changes_tab$Taxonomy = map_chr(changes_tab$MAG, function(x) taxonomy$classification[taxonomy$user_genome == x])

weighted.mean(changes_tab$median_change[changes_tab$scenario == 126] / changes_tab$median_present[changes_tab$scenario == 126], changes_tab$r2[changes_tab$scenario == 126])
# 0.1674692, sba is 21%

weighted.mean(changes_tab$median_change[changes_tab$scenario == 370] / changes_tab$median_present[changes_tab$scenario == 370], changes_tab$r2[changes_tab$scenario == 370])
# 0.455727, sba is 47%

weighted.mean(changes_tab$median_change[changes_tab$scenario == 585] / changes_tab$median_present[changes_tab$scenario == 585], changes_tab$r2[changes_tab$scenario == 585])
# 0.6599204, sba is 58%

ggplot(changes_tab, aes(x=median_present, y=median_future, colour=as.factor(scenario))) + geom_point(alpha=0.1) + geom_abline(slope = 1) + 
  geom_smooth(method='lm', se=F) + scale_x_log10() + scale_y_log10() + theme_linedraw()

summary(lm(data = changes_tab %>% filter(scenario == 126), formula = median_future ~ median_present))
# Coefficients:
#                 Estimate  Std. Error  t value Pr(>|t|)    
# (Intercept)     4.752510    2.424587    1.96   0.0501 .  
# median_present  1.060369    0.004511  235.06   <2e-16 ***
# Residual standard error: 109.4 on 2328 degrees of freedom
# Multiple R-squared:  0.9596,	Adjusted R-squared:  0.9596 
# F-statistic: 5.525e+04 on 1 and 2328 DF,  p-value: < 2.2e-16

summary(lm(data = changes_tab %>% filter(scenario == 370), formula = median_future ~ median_present))
# Coefficients:
#                 Estimate  Std. Error  t value Pr(>|t|)    
# (Intercept)     17.58818    5.47682   3.211  0.00134 **  
# median_present  1.09731    0.01019 107.684  < 2e-16 ***
# Residual standard error: 247.1 on 2328 degrees of freedom
# Multiple R-squared:  0.8328,	Adjusted R-squared:  0.8327 
# F-statistic: 1.16e+04 on 1 and 2328 DF,  p-value: < 2.2e-16

summary(lm(data = changes_tab %>% filter(scenario == 585), formula = median_future ~ median_present))
# Coefficients:
#                 Estimate  Std. Error  t value Pr(>|t|)    
# (Intercept)     25.85155     7.04743   3.668  0.00025 *** 
# median_present   1.10099     0.01311  83.966  < 2e-16 ***
# Residual standard error: 318 on 2328 degrees of freedom
# Multiple R-squared:  0.7518,	Adjusted R-squared:  0.7517 
# F-statistic: 7050 on 1 and 2328 DF,  p-value: < 2.2e-16

changes_tab$Category = 'Not significant'
changes_tab$Category[(changes_tab$wilcox_p < 0.05) & (changes_tab$median_change >= 0)] = 'Increase'
changes_tab$Category[(changes_tab$wilcox_p < 0.05) & (changes_tab$median_change < 0)] = 'Decrease'

# RCP 2.6, proportion of categories
mean(changes_tab$Category[changes_tab$scenario == 126] == 'Increase') # 0.5463519
mean(changes_tab$Category[changes_tab$scenario == 126] == 'Decrease') # 0.2957082
mean(changes_tab$Category[changes_tab$scenario == 126] == 'Not significant') # 0.1579399

# RCP 4.5, proportion of categories
mean(changes_tab$Category[changes_tab$scenario == 370] == 'Increase') # 0.5193133
mean(changes_tab$Category[changes_tab$scenario == 370] == 'Decrease') # 0.311588
mean(changes_tab$Category[changes_tab$scenario == 370] == 'Not significant') # 0.1690987

# RCP 8.5, proportion of categories
mean(changes_tab$Category[changes_tab$scenario == 585] == 'Increase') # 0.5094421
mean(changes_tab$Category[changes_tab$scenario == 585] == 'Decrease') # 0.3223176
mean(changes_tab$Category[changes_tab$scenario == 585] == 'Not significant') # 0.1682403

# Quantiles of log2FC
quantile(changes_tab$log2fc[changes_tab$scenario == 126], probs = c(0.25,0.5,0.75)) # -0.1234130  0.1302959  0.4506936
quantile(changes_tab$log2fc[changes_tab$scenario == 370], probs = c(0.25,0.5,0.75)) # -0.3159652  0.2027560  0.8329553 
quantile(changes_tab$log2fc[changes_tab$scenario == 585], probs = c(0.25,0.5,0.75)) # -0.4217135  0.1995929  0.9512945 

########################################################################################################################################
# 1. MODEL ASSESSEMENT & ABUNDANCE/SCENARIOS

# Correlation between abundance r2 values
cor.test(changes_tab$mean_rel_ab, changes_tab$r2, method = 'spearman')
#S = 4.7479e+10, p-value < 2.2e-16, rho = 0.1658951 

# Correlation between scenarios
cor.test(changes_tab$median_future[changes_tab$scenario == 370], changes_tab$median_future[changes_tab$scenario == 126], method = 'spearman')
# S = 46303085, p-value < 2.2e-16, rho = 0.9780369
cor.test(changes_tab$median_future[changes_tab$scenario == 370], changes_tab$median_future[changes_tab$scenario == 585], method = 'spearman')
# S = 9916759, p-value < 2.2e-16, rho = 0.9952962

##############################
# Plotting
# Plots models r2 (A) and models output proportion of increase/decrease/not significant
changes_tab$R2_cat = '< 0.1'
changes_tab$R2_cat[changes_tab$r2 >= 0.1] = '< 0.2'
changes_tab$R2_cat[changes_tab$r2 >= 0.2] = '< 0.3'
changes_tab$R2_cat[changes_tab$r2 >= 0.3] = '>= 0.3'

changes_tab$ab_cat = '< 1e-5 %'
changes_tab$ab_cat[changes_tab$mean_rel_ab >= 1e-5] = '< 1e-4 %'
changes_tab$ab_cat[changes_tab$mean_rel_ab >= 1e-4] = '< 1e-3 %'
changes_tab$ab_cat[changes_tab$mean_rel_ab >= 1e-3] = '>= 1e-3 %'

table(changes_tab$ab_cat, changes_tab$R2_cat)
r2_label = as.expression(bquote(""*Cross-validation~r^2*""))

p1 = ggplot(changes_tab, aes(x=..count.., y=as.character(scenario), fill = R2_cat)) + geom_bar() + facet_wrap(~scenario, scales = 'free') +
  theme_linedraw() + theme(axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title.y = element_blank()) +
  xlab('Strain count') + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_manual(name = r2_label, values = c('#CC2F50','#E3C78D','#4BC992','#1A5A61'))

changes_tab <- changes_tab %>% group_by(scenario) %>% mutate(med = quantile(log2fc, probs = 0.5),
                                                             q25 = quantile(log2fc, probs = 0.25),
                                                             q75 = quantile(log2fc, probs = 0.75))
p2 = ggplot(changes_tab, aes(x=mean_rel_ab, y=log2fc)) + geom_point(aes(fill = Category, colour = Category), alpha=0.5) + scale_x_continuous(trans = 'log10') +
  scale_fill_manual(values = c('#E04D55', '#15356B', '#CAD1E8')) + scale_colour_manual(values = c('#E04D55', '#15356B', '#CAD1E8')) +
  geom_hline(aes(yintercept = 0), colour='dimgrey', linetype='dashed') + 
  geom_hline(aes(yintercept = med, group = scenario)) + facet_wrap(~scenario) +
  geom_hline(aes(yintercept = q25, group = scenario)) + facet_wrap(~scenario) +
  geom_hline(aes(yintercept = q75, group = scenario)) + facet_wrap(~scenario) +
  theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab('Mean relative abundance') + scale_y_continuous(name = bquote(""*log[2]~fold-change*""))

p = ggarrange(p1, p2, nrow = 2, align = 'v', labels = c('A', 'B'), heights = c(0.333,0.666))
ggsave(p, filename = 'Plots/SFig_5_models_results.pdf', width = 7, height = 5)


########################################################################################################################################
# 2. PHYLOGENETIC ANALYSES

tree = read.tree('Data/treeBacteria.tree')
tree = midpoint.root(tree)
tree$tip.label = vapply(tree$tip.label, function(x) strsplit(x, '.fa')[[1]][1], FUN.VALUE = character(1))
changes_tab = changes_tab[changes_tab$MAG %in% tree$tip.label,]
tree = keep.tip(tree, tree$tip.label[tree$tip.label %in% changes_tab$MAG])
changes_tab = changes_tab %>% arrange(match(MAG, tree$tip.label))

##############################
# Changes : log2 FC
phylosig(tree, changes_tab$log2fc[changes_tab$scenario == 126], method = 'lambda', test=T)
# Phylogenetic signal lambda : 0.91958 
# logL(lambda) : -807.818
# LR(lambda=0) : 1522.21 
# P-value (based on LR test) : 0

phylosig(tree, changes_tab$log2fc[changes_tab$scenario == 370], method = 'lambda', test=T)
# Phylogenetic signal lambda : 0.915828 
# logL(lambda) : -2548.7
# LR(lambda=0) : 1567.76 
# P-value (based on LR test) : 0

phylosig(tree, changes_tab$log2fc[changes_tab$scenario == 585], method = 'lambda', test=T)
# Phylogenetic signal lambda : 0.916226 
# logL(lambda) : -3061.33
# LR(lambda=0) : 1517.79 
# P-value (based on LR test) : 0

##############################
# Predictability : r2
phylosig(tree, changes_tab$r2[changes_tab$scenario == 126], method = 'lambda', test=T)
# Phylogenetic signal lambda : 0.781509 
# logL(lambda) : 1517.03 
# LR(lambda=0) : 447.374 
# P-value (based on LR test) : 2.68953e-99

phylosig(tree, changes_tab$r2[changes_tab$scenario == 370], method = 'lambda', test=T)
# Phylogenetic signal lambda : 0.781509 
# logL(lambda) : 1517.03 
# LR(lambda=0) : 447.374 
# P-value (based on LR test) : 2.68953e-99 

phylosig(tree, changes_tab$r2[changes_tab$scenario == 585], method = 'lambda', test=T)
# Phylogenetic signal lambda : 0.781509 
# logL(lambda) : 1517.03 
# LR(lambda=0) : 447.374 
# P-value (based on LR test) : 2.68953e-99


sub_trees = subtrees(tree)
subt_res = data.frame()
for (i in 2:length(sub_trees)){
  tips = as.vector(sub_trees[[i]]$tip.label)
  categories = table(changes_tab$MAG[changes_tab$scenario == 370] %in% tips, changes_tab$Category[changes_tab$scenario == 370])
  prop_decreasing = categories[2,1] / sum(categories[2,])
  mrca = getMRCA(tree, tips)
  mrca_depth = node.depth.edgelength(tree)[mrca]
  subt_res = rbind(subt_res, data.frame(Depth = mrca_depth, Ntips = length(tips), prop_decrease = prop_decreasing, tips = paste0(tips, collapse = ',')))}
subt_res = subt_res %>% filter(Ntips >= 3, prop_decrease == 1)

# REMOVE SUBTRESS WITHIN OTHER SUBTREES
to_remove = c()
for (subt1 in subt_res$tips){
  for (subt2 in subt_res$tips){
    if (subt1 != subt2){
      subt1_tips = strsplit(subt1, ',')[[1]]
      subt2_tips = strsplit(subt2, ',')[[1]]
      if (all(subt1_tips %in% subt2_tips)){to_remove = c(to_remove, subt1)}
      if (all(subt2_tips %in% subt1_tips)){to_remove = c(to_remove, subt2)}}}}

subt_res = subt_res[!(subt_res$tips %in% to_remove),]

# CONSENSUS TAXONOMY
# source: Roland @ https://stackoverflow.com/questions/26285010/r-find-largest-common-substring-starting-at-the-beginning
fun_sbst <- function(words) {
  #extract substrings from length 1 to length of shortest word
  subs <- sapply(seq_len(min(nchar(words))), 
                 function(x, words) substring(words, 1, x), 
                 words=words)
  #max length for which substrings are equal
  neqal <- max(cumsum(apply(subs, 2, function(x) length(unique(x)) == 1L)))
  #return substring
  substring(words[1], 1, neqal)
}

subt_res$Taxonomy = NA
for (i in 1:nrow(subt_res)){
  subt_tips = strsplit(subt_res$tips[i], ',')[[1]]
  taxs = c()
  for (mag in subt_tips){taxs = c(taxs, taxonomy$classification[taxonomy$user_genome == mag])}
  subt_res$Taxonomy[i] = fun_sbst(taxs)}
subt_res$RelativeDepth = (1.622841 - subt_res$Depth) / 1.622841

quantile(subt_res$RelativeDepth)
#        0%       25%       50%       75%      100% 
# 0.1010911 0.2245602 0.2890822 0.3579994 0.5377389 

quantile(subt_res$Ntips)
#   0%   25%   50%   75%  100% 
# 3.00  3.00  4.00  5.25 14.00 

# Notably:??
#      Depth Ntips                                                                                                              Taxonomy   RelativeDepth
# 1 1.042761    12                       d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Acetobacterales;f__Acetobacteraceae;g__     0.3574470
# 2 1.375389    14              d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Chitinophagales;f__Chitinophagaceae;g__Ferruginibacter;s__     0.1524810
# 3 1.358744    13  d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Actinomycetales;f__Microbacteriaceae;g__Lacisediminihabitans;s__     0.1627373

sum(subt_res$Ntips) # 271
write.csv(subt_res, file = 'Statistics/Decreasing_clades.csv', quote = F)


##############################
changes_list = list()
changes_list$Decrease = changes_tab$MAG[(changes_tab$Category == 'Decrease') & (changes_tab$scenario == 370)]
changes_list$`Not significant and increase` = changes_tab$MAG[(changes_tab$Category != 'Decrease') & (changes_tab$scenario == 370)]

tree	<- groupOTU(tree, changes_list)

changes_tab = changes_tab[order(match(changes_tab$MAG, tree$tip.label)),]
changes_tab$Class = vapply(changes_tab$MAG, function(x) ifelse(x %in% taxonomy$user_genome, 
                                                                taxonomy$Class[taxonomy$user_genome == x], 
                                                                'Others'), FUN.VALUE = character(1))

top_class = taxonomy %>% group_by(Class) %>% summarise(n=n()) %>% top_n(12) %>% pull(Class)
changes_tab$Class[!(changes_tab$Class %in% top_class)] = 'Others'

sp1 = ggtree(tree, layout="fan", size=0.5, open.angle=5,
             aes(color = group)) + scale_colour_manual(values=c('#E04D55', '#212224')) + new_scale_colour() +
      geom_fruit(data=changes_tab[changes_tab$scenario == 370,] , geom=geom_point, 
                 mapping=aes(y = MAG, x = log2fc, colour = log2fc, size = mean_rel_ab*100),
                 offset = 0.15,
                 axis.params=list(axis = "x",text.size  = 2,hjust = 1,vjust = 0.5, nbreak=3), grid.params = list()) + 
      scale_colour_gradientn(colours = c('#E04D55', '#FA8E8E', '#CAD1E8', '#15356B'), name = bquote(""*log[2]~fold-change*"")) +
      scale_size_continuous(name='Relative abundance (%)', limits = c(0, 1.25), breaks = c(0,0.4,0.8,1.2)) +
      theme(legend.position="bottom", legend.box="vertical", legend.margin=margin())

# Part 2: Taxonomy
p2 = changes_tab %>% mutate(Class = fct_reorder(Class, log2fc)) %>% ggplot(aes(x=log2fc, fill=after_stat(x), y=Class)) + 
  geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01, quantile_lines=TRUE, quantiles=2) + theme_bw() + xlim(-5,5) +
  scale_fill_gradientn(colours = c('#E04D55', '#CAD1E8', '#15356B'), values = scales::rescale(c(-10, -5, 0, 5, 10))) +
  ylab('') + xlab(bquote(""*log[2]~fold-change*"")) + geom_vline(xintercept = 0, colour='#212224', linetype='dashed')  + theme_linedraw() + 
  theme(legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

p3 = changes_tab %>% group_by(Class) %>% summarise(abs_change = sum(median_change)/1000) %>% mutate(Class = fct_reorder(Class, abs_change)) %>%
  ggplot(aes(x=abs_change, y=Class, fill=abs_change>=0)) + ylab('') + xlab('Absolute change') +
  geom_bar(stat = 'identity') + theme_linedraw() + theme(legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c('#E04D55', '#15356B'))

sp2 = ggarrange(p2, p3, heights = c(2, 2), ncol = 1, labels = c('B','C'))

p = ggarrange(sp1, sp2, ncol = 2, nrow = 1, widths = c(0.6, 0.4), labels = c('A',''))
ggsave(p, filename = 'Plots/Fig_4_Abundance_changes.pdf', width = 9, height = 7)



########################################################################################################################################
# 3. DRIVERS (Variables selected in the models)
clean_names = c('bioclim_PC6'='PC6', 'bioclim_PC5'='PC5', 'bioclim_PC4'='PC4', 'bioclim_PC3'='PC3', 'bioclim_PC2'='PC2', 'bioclim_PC1'='PC1',
                'clim_tas'='Monthly temperature', 'clim_scd'='Annual snow cover', 'clim_pr'='Monthly precipitation',
                'gl_dist'='Distance to the glacier', 'gl_coverage'='Glacier coverage', 'gl_area'='Glacier area',
                'min_quartz'='Quartz', 'min_feldspar'='Feldspar', 'min_clays'='Clays', 'min_calcite'='Calcite',
                'nut_srp_predicted'='Soluble reactive phosphate', 'nut_din_predicted'='Dissolved inorganic nitrogen', 'chla_predicted'='Chlorophyll-a',
                'pc_water_temp_predicted'='Water temperature', 'pc_turbidity_predicted'='Turbidity', 'pc_ph_predicted'='pH', 'pc_conductivity_predicted'='Conductivity')

importance_tab = model_out %>% filter(r2 >= 0.05) %>% group_by(MAG) %>% mutate(ImpRank=rank(Freq)/length(Freq)) %>% select(MAG, vars_selection_tab, ImpRank, mean_rel_ab)

importance_tab$Category = 'Others'
importance_tab$Category[startsWith(importance_tab$vars_selection_tab, 'bioclim_')] = 'Bioclimatic'
importance_tab$Category[startsWith(importance_tab$vars_selection_tab, 'min_')] = 'Minerals'
importance_tab$Category[startsWith(importance_tab$vars_selection_tab, 'pc_')] = 'Physico-chemical'
importance_tab$Category[startsWith(importance_tab$vars_selection_tab, 'clim_')] = 'Climatic'
importance_tab$Category[startsWith(importance_tab$vars_selection_tab, 'gl_')] = 'Glaciological'

sp1 = importance_tab %>% group_by(vars_selection_tab, Category) %>% summarise(Importance = weightedMean(ImpRank, w = mean_rel_ab)) %>% 
  ggplot(aes(x=Importance*100, y=reorder(vars_selection_tab, Importance), fill=Importance)) + geom_bar(stat='identity') + 
  facet_grid(Category~., scales = 'free', space = 'free') + theme_linedraw() + ylab('') + xlab('Mean rel. rank [%]') + 
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
ggsave(p, filename = 'Plots/Fig_5_Drivers.pdf', width = 9, height = 7)


# Phylogeny vs Variables used
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

p = ggplot(dists_sub, aes(x=PD, y=Distance)) + geom_point() + geom_vline(xintercept = median_same_genus, colour='red') + 
  geom_smooth(method='loess', se=F) + theme_bw() + scale_x_log10() + xlab('Phylogenetic distance') + ylab('Kendall distance of variable ranks')
ggsave(p, filename = 'Plots/sFig_6_Distance_variable_ranks.pdf')
cor.test(dists_sub$Distance, dists_sub$PD, method='spearman') # S = 1.3771e+11, p-value < 2.2e-16, 0.1737688


































asv_tax = read.csv('RawData/16s/taxonomy.tsv', sep = '\t')
asv_tax$Genus = vapply(asv_tax$Taxon, function(x) strsplit(x, split = '; ')[[1]][6], FUN.VALUE = character(1))
asv_tax$Family = vapply(asv_tax$Taxon, function(x) strsplit(x, split = '; ')[[1]][5], FUN.VALUE = character(1))
asv_tax$Order = vapply(asv_tax$Taxon, function(x) strsplit(x, split = '; ')[[1]][4], FUN.VALUE = character(1))
asv_tax$Class = vapply(asv_tax$Taxon, function(x) strsplit(x, split = '; ')[[1]][3], FUN.VALUE = character(1))
asv_tax$Phylum = vapply(asv_tax$Taxon, function(x) strsplit(x, split = '; ')[[1]][2], FUN.VALUE = character(1))

data_1 = read.csv('Data/SDMs_results_ASVs_370.csv', sep=',')
data_2 = read.csv('Data/SDMs_results_ASVs_126_585.csv', sep=',')
data_all = rbind(data_1, data_2)

data_all = data_all %>%
  group_by(ASV, scenario) %>% 
  summarise_all(mean)

data_all = na.omit(data_all)
data_all$Model_r2[data_all$Model_r2 <= 0] = 0.0000001

data_all$r2_cat[data_all$Model_r2 >= 0.3] = '>= 0.3'
data_all$r2_cat[data_all$Model_r2 < 0.3] = '< 0.3'
data_all$r2_cat[data_all$Model_r2 < 0.2] = '< 0.2'
data_all$r2_cat[data_all$Model_r2 < 0.1] = '< 0.1'
data_all$r2_cat[data_all$Model_r2 < 0.05] = '< 0.05'

data_all$preval_cat[data_all$Prevalence >= 0.75] = '>= 0.75'
data_all$preval_cat[data_all$Prevalence < 0.75] = '< 0.75'
data_all$preval_cat[data_all$Prevalence < 0.5] = '< 0.5'
data_all$preval_cat[data_all$Prevalence < 0.3] = '< 0.3'
data_all$preval_cat[data_all$Prevalence < 0.2] = '< 0.2'

ggplot(data_all) + geom_bar(aes(fill=r2_cat, y=preval_cat), stat='count') + xlab('ASV #') + ylab('Prevalence') + 
  theme_linedraw() + facet_grid(~scenario) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave('Plotting/Models_performance.pdf', width = 7, height = 5)


###################################################
# Projection of ASVs abundance
###################################################
# 1. Community level
ra_df = data.frame(ASV_abundance = c(data_all$Abundance_predicted_median, data_all$Abundance_projected_median),
                   ASV_rank = c(rank(-data_all$Abundance_predicted_median) / nrow(data_all), rank(-data_all$Abundance_projected_median) / nrow(data_all)),
                   Date = c(rep('Present', nrow(data_all)), rep('Future', nrow(data_all))))
p2 = ggplot(data=ra_df, aes(y=ASV_abundance, x=ASV_rank, colour=Date)) + geom_line(size=3) + theme_bw() + scale_colour_manual(values = c('#4292BD','#585A55')) +
  ylab('Predicted median log(abundance+1)') + xlab('Relative rank')

p1 = ggplot(data_all) + 
  geom_point(aes(x=log1p(Abundance_predicted_median), y=log1p(Abundance_projected_median), size=Model_r2), colour='#4292BD', alpha=0.2) + 
  geom_abline(slope = 1, linetype='dashed', colour='dimgrey') +
  geom_smooth(aes(x=log1p(Abundance_predicted_median), y=log1p(Abundance_projected_median), weight=Model_r2), method='gam', colour='#585A55') + 
  theme_bw() + xlab('Predicted median log(abundance+1)') + ylab('Projected median log(abundance+1)') + xlim(0,2.5) + ylim(0,2.5) + scale_size_continuous(name='R2')
  
ggarrange(p1, p2, nrow = 1, ncol = 2, labels = c('A', 'B'))
ggsave('Plotting/Community_level.pdf', width = 10, height = 4)

# median diff in log10(abundance) = 0.07777381, Wilcoxon test: V = 124433372, p-value < 2.2e-16
median(data_all$MedianAbundanceLog10fc)
wilcox.test(data_all$MedianAbundanceLog10fc)
mean(data_all$MedianAbundanceLog10fc < 0) # 0.1897521
mean(data_all$MedianAbundanceLog10fc >= 0) # 0.8102479

# 2. Genus level
plot_data = data.frame(ASV=data_all$ASV)
plot_data$logfc = vapply(plot_data$ASV, function(x) data_all$MedianAbundanceLog10fc[data_all$ASV == x] / log1p(data_all$Abundance_predicted[data_all$ASV == x]), FUN.VALUE = numeric(1))
plot_data$Model_r2 = vapply(plot_data$ASV, function(x) data_all$Model_r2[data_all$ASV == x], FUN.VALUE = numeric(1))
plot_data$Prevalence = vapply(plot_data$ASV, function(x) data_all$Prevalence[data_all$ASV == x], FUN.VALUE = numeric(1))

plot_data$Taxonomy = vapply(plot_data$ASV, function(x) asv_tax$Taxon[asv_tax$Feature.ID == x], FUN.VALUE = character(1))
plot_data$Phylum = vapply(plot_data$Taxonomy, function(x) strsplit(x, split = '; ')[[1]][2], FUN.VALUE = character(1))

phyla_to_keep = names(table(plot_data$Phylum))[table(plot_data$Phylum) > 315]
plot_data$Phylum[!(plot_data$Phylum %in% phyla_to_keep)] = 'Others'
plot_data$Phylum = as.factor(plot_data$Phylum)

plot_data$Genus = vapply(plot_data$Taxonomy, function(x) strsplit(x, split = '; ')[[1]][6], FUN.VALUE = character(1))
plot_data$Genus[is.na(plot_data$Genus)] = 'g__na'
plot_data = plot_data[plot_data$Genus != 'g__na',]
plot_data$Genus = factor(plot_data$Genus)
plot_data$Genus = fct_reorder(plot_data$Genus, -plot_data$logfc , median)
plot_data$n_asvs = vapply(plot_data$Genus, function(x) sum(plot_data$Genus == x), FUN.VALUE = numeric(1))

plot_data$Direction = 'Increase'
plot_data$Direction[plot_data$logfc < 0] = 'Decrease'

plot_data$abundance_present = vapply(plot_data$ASV, function(x) data_all$Abundance_predicted[data_all$ASV == x], FUN.VALUE = numeric(1))
plot_data$abundance_genus_present = vapply(1:nrow(plot_data), function(i) sum(plot_data$abundance_present[plot_data$Genus == plot_data$Genus[i]]), FUN.VALUE = numeric(1))

kruskal.test(plot_data$logfc ~ as.factor(plot_data$Genus)) # Kruskal-Wallis chi-squared = 839.78, df = 474, p-value < 2.2e-16
genus_medians = plot_data %>% group_by(Genus) %>% summarise(median = weightedMedian(logfc, Model_r2), N = n()) %>% filter(N > 9)

# Genera statistical testing
compared_to_0 = data.frame()
compared_to_others = data.frame()
for (genus in unique(plot_data$Genus)[unique(plot_data$Genus) != 'g__uncultured']){
  if (sum(plot_data$Genus == genus) > 6){
    
    test_0 = wilcox.test(plot_data$logfc[plot_data$Genus == genus], weights = plot_data$Model_r2)
    compared_to_0 = rbind(compared_to_0, data.frame(Genus=genus, p=test_0$p.value, median=weightedMedian(plot_data$logfc[plot_data$Genus == genus])))
    
    test_other = wilcox.test(plot_data$logfc[plot_data$Genus == genus], plot_data$logfc[plot_data$Genus != genus],  weights = plot_data$Model_r2)
    compared_to_others = rbind(compared_to_others, data.frame(Genus=genus, p=test_other$p.value, median_diff=weightedMedian(plot_data$logfc[plot_data$Genus == genus]) - weightedMedian(plot_data$logfc[plot_data$Genus != genus])))}}
compared_to_0$padj = p.adjust(compared_to_0$p, method = 'fdr')
compared_to_others$padj = p.adjust(compared_to_others$p, method = 'fdr')

ggplot(compared_to_0, aes(x=median, y=-log(p), color=padj<0.05, label=Genus)) + geom_point() + ggrepel::geom_label_repel()
ggplot(compared_to_others, aes(x=median_diff, y=-log(p), color=padj<0.05, label=Genus)) + geom_point() + ggrepel::geom_label_repel()

mean(compared_to_0$padj < 0.05)

plot_data$Genus = as.character(plot_data$Genus)
plot_data$Phylum = as.character(plot_data$Phylum)
plot_data = plot_data[!((plot_data$n_asvs < 25 ) | (plot_data$Genus %in% c('g__uncultured', 'g__na'))),]
plot_data$Genus = vapply(plot_data$Genus, function(x) ifelse(startsWith(x, prefix = 'g__'), strsplit(x, split = 'g__')[[1]][2], x), FUN.VALUE = character(1))
plot_data$Phylum = vapply(plot_data$Phylum, function(x) ifelse(startsWith(x, prefix = 'p__'), strsplit(x, split = 'p__')[[1]][2], x), FUN.VALUE = character(1))

p1 = plot_data  %>%
  ggplot(aes(x=logfc, y=Genus, fill=Phylum, height = ..count.., weight=Model_r2**2)) + geom_density_ridges(stat = "density", scale=0.8, trim = TRUE) + 
  ylab('') + xlab('Relative change in log(abundance)') + geom_vline(xintercept = median(plot_data$logfc), linetype='dashed') + 
  theme_bw() + scale_x_continuous(trans = pseudolog10_trans) + facet_grid(Phylum~., scales = 'free', space = 'free') + 
  scale_fill_manual(values=c("#FAD72A","#E3985B","dimgrey","#FA7085", "#A45BE3","#4288FF","#3FFFC0")) + theme(strip.background = element_blank(),
                                                                                                              strip.text = element_blank(),
                                                                                                              legend.position = 'top',
                                                                                                              axis.text.y = element_text(face = "italic"))
l1 = get_legend(p1)

p2 = plot_data  %>%
  dplyr::count(Genus, Direction, Phylum) %>% group_by(Genus) %>%  mutate(Proportion = n / sum(n), Rank = mean(rank)) %>% 
  ggplot(aes(y=Genus, x=Proportion, fill=Direction)) + geom_bar(position = position_fill(reverse = TRUE), stat = "identity") +
  facet_grid(Phylum~., scales = 'free', space = 'free') + xlab('Prop.') + ylab('') + scale_fill_manual(values=c('#9E2D37','#326E8F')) + 
  scale_x_continuous(breaks=c(0, 1)) + theme_bw() + theme(strip.background = element_blank(),
                                                          strip.text = element_blank(),
                                                          axis.text.y = element_blank(),
                                                          legend.position = 'top')
l2 = get_legend(p2)
leg = ggarrange(l1, l2, nrow = 2, ncol = 1)

r1 = ggarrange(p1 + theme(legend.position="none") , p2 + theme(legend.position="none"), ncol = 2, nrow = 1, labels = c('A','B'), widths = c(0.87, 0.13))
p = ggarrange(r1, leg, nrow = 2, ncol = 1, heights = c(0.8, 0.2))
ggsave('Plotting/Genera_logfc_decrease.pdf', width = 6.7, height = 7.5)







###################################################
# Drivers of bacterial abundance
###################################################
SameTax <- function(tax1, tax2){
  if ((tax1 == tax2) &
     (!(is.na(tax1) | is.na(tax2))) &
     (!(grepl('uncultured', tax1)) & !(grepl('uncultured', tax2)))){return('Yes')}
  else{return('No')}}
asv_comp = expand.grid(ASV1 = data_all$ASV[data_all$Model_r2 >= 0.2], 
                       ASV2 = data_all$ASV[data_all$Model_r2 >= 0.2])
asv_comp = asv_comp[asv_comp$ASV1 != asv_comp$ASV2,]
asv_comp = asv_comp[sample(1:nrow(asv_comp), 20000),]
asv_comp$weight = vapply(1:nrow(asv_comp), function(i) data_all$Model_r2[data_all$ASV == asv_comp$ASV1[i]] * data_all$Model_r2[data_all$ASV == asv_comp$ASV2[i]], FUN.VALUE = numeric(1))
asv_comp$ImpDist_present = vapply(1:nrow(asv_comp), function(i) cor(cbind(as.numeric(data_all[data_all$ASV == asv_comp$ASV1[i], startsWith(colnames(data_all), 'present_VarImp')]),
                                                                            as.numeric(data_all[data_all$ASV == asv_comp$ASV2[i], startsWith(colnames(data_all), 'present_VarImp')])), method = 'kendall')[1,2],
                             FUN.VALUE = numeric(1))
asv_comp$ImpDist_future = vapply(1:nrow(asv_comp), function(i) cor(cbind(as.numeric(data_all[data_all$ASV == asv_comp$ASV1[i], startsWith(colnames(data_all), 'future_VarImp')]),
                                                                            as.numeric(data_all[data_all$ASV == asv_comp$ASV2[i], startsWith(colnames(data_all), 'future_VarImp')])), method = 'kendall')[1,2],
                               FUN.VALUE = numeric(1))

asv_comp$same_genus = vapply(1:nrow(asv_comp), function(i) SameTax(asv_tax$Genus[asv_tax$Feature.ID == asv_comp$ASV1[i]], 
                                                                   asv_tax$Genus[asv_tax$Feature.ID == asv_comp$ASV2[i]]), FUN.VALUE = character(1))
asv_comp$same_family = vapply(1:nrow(asv_comp), function(i) SameTax(asv_tax$Family[asv_tax$Feature.ID == asv_comp$ASV1[i]], 
                                                                    asv_tax$Family[asv_tax$Feature.ID == asv_comp$ASV2[i]]), FUN.VALUE = character(1))
asv_comp$same_order = vapply(1:nrow(asv_comp), function(i) SameTax(asv_tax$Order[asv_tax$Feature.ID == asv_comp$ASV1[i]], 
                                                                   asv_tax$Order[asv_tax$Feature.ID == asv_comp$ASV2[i]]), FUN.VALUE = character(1))
asv_comp$same_class = vapply(1:nrow(asv_comp), function(i) SameTax(asv_tax$Class[asv_tax$Feature.ID == asv_comp$ASV1[i]], 
                                                                   asv_tax$Class[asv_tax$Feature.ID == asv_comp$ASV2[i]]), FUN.VALUE = character(1))
asv_comp$same_phylum = vapply(1:nrow(asv_comp), function(i) SameTax(asv_tax$Phylum[asv_tax$Feature.ID == asv_comp$ASV1[i]], 
                                                                    asv_tax$Phylum[asv_tax$Feature.ID == asv_comp$ASV2[i]]), FUN.VALUE = character(1))

asv_comp$Shared_taxonomic_level = 'Domain'
asv_comp$Shared_taxonomic_level[asv_comp$same_phylum == 'Yes'] = 'Phylum'
asv_comp$Shared_taxonomic_level[asv_comp$same_class == 'Yes'] = 'Class'
asv_comp$Shared_taxonomic_level[asv_comp$same_order == 'Yes'] = 'Order'
asv_comp$Shared_taxonomic_level[asv_comp$same_family == 'Yes'] = 'Family'
asv_comp$Shared_taxonomic_level[asv_comp$same_genus == 'Yes'] = 'Genus'
table(asv_comp$Shared_taxonomic_level)

ggplot(na.omit(asv_comp), aes(x=factor(Shared_taxonomic_level, levels= c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus')), 
                              y=ImpDist_present, fill=Shared_taxonomic_level, weight=weight)) + 
  geom_boxplot() + stat_compare_means(comparisons = list(c('Domain', 'Phylum'), 
                                                         c('Phylum', 'Class'), 
                                                         c('Class', 'Order'), 
                                                         c('Order', 'Family'),
                                                         c('Family', 'Genus'))) + ylab('Feature importance correlation') + xlab('') + theme_bw() 

ggplot(na.omit(asv_comp), aes(x=factor(Shared_taxonomic_level, levels= c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus')), 
                              y=ImpDist_future, fill=Shared_taxonomic_level, weight=weight)) + 
  geom_boxplot() + stat_compare_means(comparisons = list(c('Domain', 'Phylum'), 
                                                         c('Phylum', 'Class'), 
                                                         c('Class', 'Order'), 
                                                         c('Order', 'Family'),
                                                         c('Family', 'Genus'))) + ylab('Feature importance correlation') + xlab('') + theme_bw() 

# to add facets for scenarios                                                                                  
# Importance of features
data_rcp45 = data_all[data_all$scenario == 370,]
feat_imp = melt(cbind(data_rcp45[, startsWith(colnames(data_rcp45), 'present_VarImp')], ASV=data_rcp45$ASV), 'ASV')
feat_imp = feat_imp[feat_imp$value > 0,]
feat_imp$model_r2 = vapply(feat_imp$ASV, function(x) data_rcp45$Model_r2[data_rcp45$ASV == x], FUN.VALUE = numeric(1))

feat_imp_rel =  feat_imp %>%
    group_by(variable) %>%    
    summarise(sum = sum(value))
feat_imp_rel$sum = feat_imp_rel$sum / sum(feat_imp_rel$sum)
p = ggplot(feat_imp_rel, aes(y = reorder(variable, sum), x = sum)) + 
           geom_bar(stat = 'identity') + theme_bw() + xlab('Average relative importance') + ylab('')
ggsave(p, file ='Plotting/Feature_importance_RCP45.pdf')

# Future importance of features
feat_imp_pres = melt(cbind(data_rcp45[, startsWith(colnames(data_rcp45), 'present_VarImp')], ASV=data_rcp45$ASV), 'ASV')
feat_imp_pres = feat_imp_pres[feat_imp_pres$value > 0,]
feat_imp_pres$model_r2 = vapply(feat_imp_pres$ASV, function(x) data_rcp45$Model_r2[data_rcp45$ASV == x], FUN.VALUE = numeric(1))
feat_imp_pres$id = paste0(feat_imp_pres$ASV, feat_imp_pres$variable)

feat_imp_futu = melt(cbind(data_rcp45[, startsWith(colnames(data_rcp45), 'future_VarImp')], ASV=data_rcp45$ASV), 'ASV')
feat_imp_futu = feat_imp_futu[feat_imp_futu$value > 0,]
feat_imp_futu$model_r2 = vapply(feat_imp_futu$ASV, function(x) data_rcp45$Model_r2[data_rcp45$ASV == x], FUN.VALUE = numeric(1))
feat_imp_futu$id = paste0(feat_imp_futu$ASV, feat_imp_futu$variable)
feat_imp_futu$id = gsub('future', 'present', feat_imp_futu$id)

feat_imp_futu = feat_imp_futu[feat_imp_futu$id %in% feat_imp_pres$id,]
feat_imp_pres = feat_imp_pres[feat_imp_pres$id %in% feat_imp_futu$id,]

feat_imp_time_comp = data.frame(ASV = feat_imp_pres$ASV, r2 = feat_imp_pres$model_r2, Var = feat_imp_pres$variable, Diff = (feat_imp_futu$value - feat_imp_pres$value))

erfc = function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
wwilcox = function(x, wx){
  # https://stackoverflow.com/questions/24648052/weighted-wilcoxon-test-in-r
  U = 0
  ## Loop over the selection branches
  y = rep(0, length(x))
  for( iy in y ){
    ## Neutral branches smaller or equal
    smaller = which(  x > iy )
    equal = which( x == iy )
    
    ## Count
    sumSmaller = sum(wx[smaller])
    sumEqual = sum(wx[equal]/2)
    sumTot = sumSmaller + sumEqual
    
    ## Total rank
    U = U + sumTot}
  
  ## U statistics
  nY = length(y)
  nX = sum(wx)
  
  ## Large sample: U follows a Gaussian
  mU = nY * nX / 2
  sigU = sqrt( ( nY * nX * ( 1 + nY + nX ) ) / 12 )
  zU = ( U - mU ) / sigU
  
  ## p-value, one-sided
  pU = erfc( zU / sqrt(2) ) /2
  
  return(pU)}
vars_to_plot = feat_imp_time_comp %>% group_by(Var) %>% summarize(p = wilcox.test(Diff)$p.value, median = median(Diff)) %>% filter(p < 0.05) %>% arrange(median)
                                                                       
# to add facets for scenarios                                                                                  
ggplot(feat_imp_time_comp, aes(x=Diff, y=reorder(Var, Diff), weight=r2)) + geom_boxplot() + theme_bw() + 
  xlab('Difference in feature importance (future - present)') + ylab('')
ggsave('Plotting/Difference_importance.pdf')
















features = as.character(vapply(colnames(data_all)[startsWith(colnames(data_all), 'present_VarImp')], function(x) strsplit(x, split='VarImp')[[1]][2], FUN.VALUE = character(1)))

imp_df = expand.grid(Feature=features, ASV=data_all$ASV)
imp_df$Importance = vapply(1:nrow(imp_df), function(i) data_all[data_all$ASV == imp_df$ASV[i], paste0('present_VarImp', imp_df$Feature[i], collapse = '')] %>% pull(), FUN.VALUE = numeric(1))
imp_df$Taxonomy = vapply(imp_df$ASV, function(x) asv_tax$Taxon[asv_tax$Feature.ID == x], FUN.VALUE = character(1))
imp_df$r2 = vapply(imp_df$ASV, function(x) data_all$Model_r2[data_all$ASV == x], FUN.VALUE = numeric(1))
imp_df$Phylum = vapply(imp_df$Taxonomy, function(x) strsplit(x, split = '; ')[[1]][2], FUN.VALUE = character(1))
phyla_to_keep = names(table(imp_df$Phylum))[table(imp_df$Phylum) > 1000]
imp_df$Phylum[!(imp_df$Phylum %in% phyla_to_keep)] = 'Others'
imp_df$Phylum = as.factor(imp_df$Phylum)


acido_best_vars =  imp_df[(imp_df$Phylum == 'p__Acidobacteriota'),] %>% group_by(Feature) %>% summarise(mean(Importance)) %>% filter(`mean(Importance)` > quantile(`mean(Importance)`, 0.85)) %>% pull(Feature)
p_acido = imp_df[(imp_df$Phylum == 'p__Acidobacteriota') & (imp_df$Feature %in% acido_best_vars),] %>% mutate(Feature = reorder(Feature, Importance, mean)) %>%
  ggplot(aes(y=Feature,x=Importance)) + geom_bar(aes(weight=r2), fill = '#FF7521', position = "dodge", stat = "summary", fun.y = "mean") + 
  theme_bw() + ggtitle('Acidobacteriota') + ylab('') + xlab('Mean relative imp.')

actin_best_vars =  imp_df[(imp_df$Phylum == 'p__Actinobacteriota'),] %>% group_by(Feature) %>% summarise(mean(Importance)) %>% filter(`mean(Importance)` > quantile(`mean(Importance)`, 0.85)) %>% pull(Feature)
p_actin = imp_df[(imp_df$Phylum == 'p__Actinobacteriota') & (imp_df$Feature %in% actin_best_vars),] %>% mutate(Feature = reorder(Feature, Importance, mean)) %>%
  ggplot(aes(y=Feature,x=Importance)) + geom_bar(aes(weight=r2), fill = '#B30AFB', position = "dodge", stat = "summary", fun.y = "mean") + 
  theme_bw() + ggtitle('Actinobacteriota') + ylab('') + xlab('Mean relative imp.')

bacte_best_vars =  imp_df[(imp_df$Phylum == 'p__Bacteroidota'),] %>% group_by(Feature) %>% summarise(mean(Importance)) %>% filter(`mean(Importance)` > quantile(`mean(Importance)`, 0.85)) %>% pull(Feature)
p_bacte = imp_df[(imp_df$Phylum == 'p__Bacteroidota') & (imp_df$Feature %in% bacte_best_vars),] %>% mutate(Feature = reorder(Feature, Importance, mean)) %>% 
  ggplot(aes(y=Feature,x=Importance)) + geom_bar(aes(weight=r2), fill = '#FB0A22', position = "dodge", stat = "summary", fun.y = "mean") + 
  theme_bw() + ggtitle('Bacteroidota') + ylab('') + xlab('Mean relative imp.')

chlor_best_vars =  imp_df[(imp_df$Phylum == 'p__Chloroflexi'),] %>% group_by(Feature) %>% summarise(mean(Importance)) %>% filter(`mean(Importance)` > quantile(`mean(Importance)`, 0.85)) %>% pull(Feature)
p_chlor = imp_df[(imp_df$Phylum == 'p__Chloroflexi') & (imp_df$Feature %in% chlor_best_vars),] %>% mutate(Feature = reorder(Feature, Importance, mean)) %>% 
  ggplot(aes(y=Feature,x=Importance)) + geom_bar(aes(weight=r2), fill = '#23DED5', position = "dodge", stat = "summary", fun.y = "mean") + 
  theme_bw() + ggtitle('Chloroflexi') + ylab('') + xlab('Mean relative imp.')

pates_best_vars =  imp_df[(imp_df$Phylum == 'p__Patescibacteria'),] %>% group_by(Feature) %>% summarise(mean(Importance)) %>% filter(`mean(Importance)` > quantile(`mean(Importance)`, 0.85)) %>% pull(Feature)
p_pates = imp_df[(imp_df$Phylum == 'p__Patescibacteria') & (imp_df$Feature %in% pates_best_vars),] %>% mutate(Feature = reorder(Feature, Importance, mean)) %>% 
  ggplot(aes(y=Feature,x=Importance)) + geom_bar(aes(weight=r2), fill = '#E6C619', position = "dodge", stat = "summary", fun.y = "mean") + 
  theme_bw() + ggtitle('Patescibacteria') + ylab('') + xlab('Mean relative imp.')

planc_best_vars =  imp_df[(imp_df$Phylum == 'p__Planctomycetota'),] %>% group_by(Feature) %>% summarise(mean(Importance)) %>% filter(`mean(Importance)` > quantile(`mean(Importance)`, 0.85)) %>% pull(Feature)
p_planc = imp_df[(imp_df$Phylum == 'p__Planctomycetota') & (imp_df$Feature %in% planc_best_vars),] %>% mutate(Feature = reorder(Feature, Importance, mean)) %>% 
  ggplot(aes(y=Feature,x=Importance)) + geom_bar(aes(weight=r2), fill = '#21CF6F', position = "dodge", stat = "summary", fun.y = "mean") + 
  theme_bw() + ggtitle('Planctomycetota') + ylab('') + xlab('Mean relative imp.')

prote_best_vars =  imp_df[(imp_df$Phylum == 'p__Proteobacteria'),] %>% group_by(Feature) %>% summarise(mean(Importance)) %>% filter(`mean(Importance)` > quantile(`mean(Importance)`, 0.85)) %>% pull(Feature)
p_prote = imp_df[(imp_df$Phylum == 'p__Proteobacteria') & (imp_df$Feature %in% prote_best_vars),] %>% mutate(Feature = reorder(Feature, Importance, mean)) %>% 
  ggplot(aes(y=Feature,x=Importance)) + geom_bar(aes(weight=r2), fill = '#4A64E6', position = "dodge", stat = "summary", fun.y = "mean") + 
  theme_bw() + ggtitle('Proteobacteria') + ylab('') + xlab('Mean relative imp.')

verru_best_vars =  imp_df[(imp_df$Phylum == 'p__Verrucomicrobiota'),] %>% group_by(Feature) %>% summarise(mean(Importance)) %>% filter(`mean(Importance)` > quantile(`mean(Importance)`, 0.85)) %>% pull(Feature)
p_verru = imp_df[(imp_df$Phylum == 'p__Verrucomicrobiota') & (imp_df$Feature %in% verru_best_vars),] %>% mutate(Feature = reorder(Feature, Importance, mean)) %>% 
  ggplot(aes(y=Feature,x=Importance)) + geom_bar(aes(weight=r2), fill = '#FF4F70', position = "dodge", stat = "summary", fun.y = "mean") + 
  theme_bw() + ggtitle('Verrucomicrobiota') + ylab('') + xlab('Mean relative imp.')

other_best_vars =  imp_df[(imp_df$Phylum == 'Others'),] %>% group_by(Feature) %>% summarise(mean(Importance)) %>% filter(`mean(Importance)` > quantile(`mean(Importance)`, 0.85)) %>% pull(Feature)
p_other = imp_df[(imp_df$Phylum == 'Others') & (imp_df$Feature %in% other_best_vars),] %>% mutate(Feature = reorder(Feature, Importance, mean)) %>% 
  ggplot(aes(y=Feature,x=Importance)) + geom_bar(aes(weight=r2), fill = 'dimgrey', position = "dodge", stat = "summary", fun.y = "mean") + 
  theme_bw() + ggtitle('Others') + ylab('') + xlab('Mean relative imp.')

p = ggarrange(p_acido, p_actin, p_bacte, p_chlor, p_pates, p_planc, p_prote, p_verru, p_other, align = 'v', ncol=2, nrow = 5)
ggsave(p, filename = 'Plotting/FeatureImportance.pdf', width = 7, height = 10)



ggplot(imp_df, aes(y=Feature,x=Importance)) + geom_bar(stat = "summary", fun.y = "median")


data_all$Taxonomy = vapply(data_all$ASV, function(x) asv_tax$Taxon[asv_tax$Feature.ID == x], FUN.VALUE = character(1))
data_all$Phylum = vapply(data_all$Taxonomy, function(x) strsplit(x, split = '; ')[[1]][2], FUN.VALUE = character(1))
phyla_to_keep = names(table(data_all$Phylum))[table(data_all$Phylum) > 10]
data_all$Phylum[!(data_all$Phylum %in% phyla_to_keep)] = 'Others'
data_all$Phylum = as.factor(data_all$Phylum)

ggplot(data_all, aes(x=future_VarImppredicted_chla-present_VarImppredicted_chla, y=Phylum)) + geom_boxplot()


ggplot(imp_df, aes(y=Feature,x=as.integer(Importance>0))) + geom_bar(stat = "summary", fun.y = "median")



pdists = cophenetic(tree)
data_all$log1p_val = log1p(data_all$MedianAbundanceDiff / data_all$TotalAbundancePresent_median)

pairs_data = expand.grid(data_all$ASV, data_all$ASV)
pairs_data = pairs_data[pairs_data$Var1 != pairs_data$Var2,]
pairs_data = pairs_data[sample(1:nrow(pairs_data)),]
pairs_data$DeltaDiff = vapply(1:nrow(pairs_data), function(i) abs(data_all$log1p_val[data_all$ASV == pairs_data$Var1[i]] -
                                                                  data_all$log1p_val[data_all$ASV == pairs_data$Var2[i]]), FUN.VALUE = numeric(1))
pairs_data$PD = vapply(1:nrow(pairs_data), function(i) pdists[pairs_data$Var1[i], pairs_data$Var2[i]], FUN.VALUE = numeric(1))

ggplot(pairs_data, aes(x=PD, y=DeltaDiff)) + geom_point() + geom_smooth()





px = ggplot(pairs_data, aes(x=Var1,y=Var2,fill=FeatCor,alpha=R2_geom_mean)) + geom_tile() + facet_grid(P2~P1, scales = 'free', space = 'free') +
  scale_fill_gradient(low = 'salmon', high = 'turquoise') + scale_alpha_continuous(range = c(0.3,1)) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())
ggsave(px, filename = 'Plotting/test_heatmap.png', width = 11, height = 10)

data_all$PropChange = data_all$MedianAbundanceDiff / data_all$TotalAbundancePresent_median
data_all = data_all[data_all$Phylum %in% names(table(data_all$Phylum))[table(data_all$Phylum) >= 10],]

ggplot(data_all, aes(x=TotalAbundancePresent_median,y=TotalAbundanceFuture_median,colour=Phylum)) + geom_point() + geom_smooth(method='lm', se=F) + geom_abline(slope = 1)
ggplot(data_all, aes(x=PropChange, y=Phylum, fill=Phylum, weight=Model_r2)) + geom_boxplot() + scale_x_log10()






CustomDist <- function(v1, v2){
  shared = sum(v1[(v1 == v2) & (v1 == 1)])
  tot = min(c(sum(v1), sum(v2)))
  return(shared/tot)}

pairs_data = expand.grid(data_all$ASV, data_all$ASV)
pairs_data = pairs_data[pairs_data$Var1 != pairs_data$Var2,]
pairs_data = pairs_data[sample(1:nrow(pairs_data), 10000),]
phylo_dists = cophenetic(tree)
pairs_data$PD = vapply(1:nrow(pairs_data), function(i) phylo_dists[pairs_data$Var1[i], pairs_data$Var2[i]], FUN.VALUE = numeric(1))
feats = colnames(data_all)[startsWith(colnames(data_all), 'present_VarImp')]
pairs_data$FeatD = vapply(1:nrow(pairs_data), function(i) CustomDist(data_all[data_all$ASV == pairs_data$Var1[i], feats] > 0,
                                                                     data_all[data_all$ASV == pairs_data$Var2[i], feats] > 0), FUN.VALUE = numeric(1))

ggplot(pairs_data, aes(x=PD, y=FeatD)) + geom_point() + geom_smooth()

