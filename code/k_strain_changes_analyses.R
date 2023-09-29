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
library(ggridges)

LoadModelOutput <- function(){
  model_out = read.csv('stats/model_res_all_scenarios_final.csv')
  min_non_zero = min(model_out$median_future[model_out$median_future > 0])
  model_out$wilcox_p[is.na(model_out$wilcox_p)] = 1 # NAs for models that selected only variables that do not change (0 median change, etc.)
  model_out$wilcox_p = p.adjust(model_out$wilcox_p, method='holm')
  model_out$median_future[model_out$median_future < 0] = min_non_zero/2
  model_out$median_present[model_out$median_present < 0] = min_non_zero/2
  model_out$log2fc = log2(model_out$median_future) - log2(model_out$median_present)

  model_out$Category = 'Not significant'
  model_out$Category[(model_out$wilcox_p < 0.05) & (model_out$median_change > 0)] = 'Increase'
  model_out$Category[(model_out$wilcox_p < 0.05) & (model_out$median_change < 0)] = 'Decrease'

  taxonomy = read.csv('data/raw/bacteria/gtdbtk.bac120.summary.tsv', sep='\t')
  model_out = model_out[model_out$MAG %in% taxonomy$user_genome,]
  taxonomy$Phylum = vapply(taxonomy$classification, function(x) strsplit( strsplit(x, ';c__')[[1]][1], ';p__')[[1]][2], FUN.VALUE = character(1))
  taxonomy$Class = vapply(taxonomy$classification, function(x) strsplit( strsplit(x, ';o__')[[1]][1], ';c__')[[1]][2], FUN.VALUE = character(1))
  taxonomy$Order = vapply(taxonomy$classification, function(x) strsplit( strsplit(x, ';f__')[[1]][1], ';o__')[[1]][2], FUN.VALUE = character(1))
  taxonomy$Family = vapply(taxonomy$classification, function(x) strsplit( strsplit(x, ';g__')[[1]][1], ';f__')[[1]][2], FUN.VALUE = character(1))
  taxonomy$Genus = vapply(taxonomy$classification, function(x) strsplit( strsplit(x, ';s__')[[1]][1], ';g__')[[1]][2], FUN.VALUE = character(1))

  changes_tab = model_out %>% select(-vars_selection_tab, -Freq, scenario) %>% distinct()
  changes_tab$Taxonomy = map_chr(changes_tab$MAG, function(x) ifelse(x %in% taxonomy$user_genome, taxonomy$classification[taxonomy$user_genome == x], ';'))
  write.csv(changes_tab, 'stats/changes_tab_all_scenarios_final.csv', quote=F, row.names=F)
  return(list(changes=changes_tab, taxonomy=taxonomy))}

LoadTree <- function(){
  tree = read.tree('data/raw/bacteria/treeBacteria.tree')
  tree = midpoint.root(tree)
  tree$tip.label = vapply(tree$tip.label, function(x) strsplit(x, '.fa')[[1]][1], FUN.VALUE = character(1))
  return(tree)}

CompareScenariosChanges <- function(changes_tab){
  for (scenario in c(126, 370, 585)){
    msum = summary(lm(data = changes_tab %>% filter(scenario == scenario), formula = median_future ~ median_present))
    capture.output(msum, file = paste0('stats/changes_scenario_comparison_',scenario,'.txt', collapse = ''))}

  changes_tab$scenario[changes_tab$scenario == 126] = 'RCP2.6'
  changes_tab$scenario[changes_tab$scenario == 370] = 'RCP4.5'
  changes_tab$scenario[changes_tab$scenario == 585] = 'RCP8.5'
  p = ggplot(changes_tab, aes(x=median_present, y=median_future, colour=scenario)) + geom_point(alpha=0.1) + geom_abline(slope = 1) + 
    geom_smooth(method='lm', se=F) + scale_x_log10() + scale_y_log10() + theme_bw() +  stat_regline_equation() + 
    xlab('Median present predicted abundance')  + ylab('Median future projected abundance') + labs(colour='Scenario')
  ggsave(p, filename = 'plots/Fig_S7_compare_scenarios_evenness.pdf', width=6, height=4)}

ProportionIncrease <- function(changes_tab){
  prop_increase_tab = data.frame()
  for (scenario in c(126, 370, 585)){
    mean_incr = mean(changes_tab$Category[changes_tab$scenario == scenario] == 'Increase')
    mean_decr = mean(changes_tab$Category[changes_tab$scenario == scenario] == 'Decrease')
    mean_nots = mean(changes_tab$Category[changes_tab$scenario == scenario] == 'Not significant')
    quant = quantile(changes_tab$log2fc[changes_tab$scenario == scenario], probs = c(0.25,0.5,0.75))
    prop_increase_tab = rbind(prop_increase_tab, data.frame(Scenario=scenario, MeanIncrease=mean_incr, MeanDecrease=mean_decr, MeanNotSignificant=mean_nots,
                                                            Q25=quant[1], Median=quant[2], Q75=quant[3]))}
  write.csv(prop_increase_tab, 'stats/proportion_increase.csv', quote=F, row.names=F)}

PlotModelPerformances <- function(changes_tab){
  changes_tab$R2_cat = '< 0.1'
  changes_tab$R2_cat[changes_tab$r2 >= 0.1] = '< 0.2'
  changes_tab$R2_cat[changes_tab$r2 >= 0.2] = '< 0.3'
  changes_tab$R2_cat[changes_tab$r2 >= 0.3] = '>= 0.3'

  changes_tab$ab_cat = '< 1e-5 %'
  changes_tab$ab_cat[changes_tab$mean_rel_ab >= 1e-5] = '< 1e-4 %'
  changes_tab$ab_cat[changes_tab$mean_rel_ab >= 1e-4] = '< 1e-3 %'
  changes_tab$ab_cat[changes_tab$mean_rel_ab >= 1e-3] = '>= 1e-3 %'

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
  ggsave(p, filename = 'plots/Fig_S6_models_performances.pdf', width = 7, height = 5)
  write.csv(table(changes_tab$ab_cat, changes_tab$R2_cat), 'stats/table_r2_abundance.csv', quote=F, row.names=F)}

PhyloSignal <- function(changes_tab, tree){
  sign_tab = data.frame()
  for (scenario in c(126, 370, 585)){
    sub_data = changes_tab %>% filter(scenario == 126)
    for (var in c('log2fc', 'r2')){
      phylo_test = phylosig(tree, sub_data[,var], method = 'lambda', test=T)
      pval = phylo_test$P
      lambda = phylo_test$lambda
      logl = phylo_test$logL
      logl0 = phylo_test$logL0
      sign_tab = rbind(sign_tab, data.frame(Scenario=scenario, Variable=var, p=pval, logl=logl, logl0=logl0, lambda=lambda))}}
  write.csv(sign_tab, 'stats/phylogenetic_signal.csv', quote=F, row.names=F)}

fun_sbst <- function(words) { # source: Roland @ https://stackoverflow.com/questions/26285010/r-find-largest-common-substring-starting-at-the-beginning
    #extract substrings from length 1 to length of shortest word
    subs <- sapply(seq_len(min(nchar(words))), 
                  function(x, words) substring(words, 1, x), 
                  words=words)
    #max length for which substrings are equal
    neqal <- max(cumsum(apply(subs, 2, function(x) length(unique(x)) == 1L)))
    #return substring
    substring(words[1], 1, neqal)}

MonophyleticClades <- function(tree, changes_tab, taxonomy){
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

  subt_res$Taxonomy = NA
  for (i in 1:nrow(subt_res)){
    subt_tips = strsplit(subt_res$tips[i], ',')[[1]]
    taxs = c()
    for (mag in subt_tips){taxs = c(taxs, taxonomy$classification[taxonomy$user_genome == mag])}
    subt_res$Taxonomy[i] = fun_sbst(taxs)}
  subt_res$RelativeDepth = (1.622841 - subt_res$Depth) / 1.622841

  sink(file = 'stats/mono_decreasing_clades_description.txt')
  print(paste0('Number of monophyletic decreasing clades with at least 3 members: ', nrow(subt_res)))
  print(paste0('Number of strains included: ', sum(subt_res$Ntips)))
  print('Relative depth of monophyletic decreasing clades: ')
  print(quantile(subt_res$RelativeDepth))
  print('Number of tips of monophyletic decreasing clades: ')
  print(quantile(subt_res$Ntips))
  sink(file = NULL)

  write.csv(subt_res, file = 'stats/decreasing_clades.csv', quote = F)}

PlotTreeChanges <- function(changes_tab, tree, taxonomy){
  changes_list = list()
  changes_list$Decrease = changes_tab$MAG[(changes_tab$Category == 'Decrease') & (changes_tab$scenario == 370)]
  changes_list$`Not significant and increase` = changes_tab$MAG[(changes_tab$Category != 'Decrease') & (changes_tab$scenario == 370)]

  tree	<- groupOTU(tree, changes_list)

  changes_tab = changes_tab[order(match(changes_tab$MAG, tree$tip.label)),]
  changes_tab$Class = vapply(changes_tab$MAG, function(x) ifelse(x %in% taxonomy$user_genome, 
                                                                  taxonomy$Class[taxonomy$user_genome == x], 
                                                                  'Others'), FUN.VALUE = character(1))

  top_class = taxonomy %>% group_by(Class) %>% summarise(n=n()) %>% top_n(11) %>% pull(Class)
  changes_tab$Class[!(changes_tab$Class %in% top_class)] = 'Others'
  print(quantile(changes_tab$log2fc))

  class_cols = c('#4A2B94','#3A6EF0','#93CDFA','#A957B3', # blue
                 '#256316','#61CFC7','#AAC4A5','#E3E8F0', # green
                 '#96022E','#ED1C09','#E0C16A','#BF5915') # red
  sp1 = ggtree(tree, layout="fan", size=0.2, open.angle=5,
              aes(color = group)) + labs(color='Group') + scale_colour_manual(values=c('#E04D55', '#212224')) + new_scale_colour() +
        geom_fruit(data=changes_tab[changes_tab$scenario == 370,], pwidth	= 0.04, geom=geom_bar, offset = 0.02,
                   mapping=aes(y=MAG,	x = 5, fill = Class), orientation="y", stat="identity") +
        labs(fill='') + guides(fill = guide_legend(ncol = 3)) + scale_fill_manual(values=class_cols) + new_scale_colour() + new_scale_fill() +
        geom_fruit(data=changes_tab[changes_tab$scenario == 370,] , geom=geom_point, 
                  mapping=aes(y = MAG, x = log2fc, colour = log2fc, size = mean_rel_ab*100),
                  offset = 0.3, pwidth=0.4,
                  axis.params=list(axis = "x",text.size  = 2,hjust = 1,vjust = 0.5, nbreak=3), grid.params = list()) + 
        scale_colour_gradientn(colours = c('#E04D55', '#CAD1E8', '#15356B'), name = bquote(""*log[2]~fold-change*"")) +
        scale_size_continuous(name='Present relative abundance (%)', limits = c(0, 1.25), breaks = c(0.4,0.8,1.2)) +
        theme(legend.key.size = unit(0.2, 'cm'), legend.position="bottom", legend.box="vertical", legend.margin=margin(),
              legend.key.width = unit(0.5, "cm"))

  p2 = changes_tab %>% mutate(Class = fct_reorder(Class, log2fc)) %>% ggplot(aes(x=log2fc, fill=after_stat(x), y=Class)) + 
    geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01, quantile_lines=TRUE, quantiles=2) + theme_bw() + xlim(-4.35,5.3) +
    scale_fill_gradientn(colours = c('#E04D55', '#CAD1E8', '#15356B'), values = scales::rescale(c(-10, -5, 0, 5, 10))) +
    ylab('') + xlab(bquote(""*log[2]~fold-change*"")) + geom_vline(xintercept = 0, colour='#212224', linetype='dashed')  + theme_linedraw() + 
    theme(legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

  p3 = changes_tab %>% group_by(Class) %>% summarise(abs_change = sum(median_change)/1000) %>% mutate(Class = fct_reorder(Class, abs_change)) %>%
    ggplot(aes(x=abs_change, y=Class, fill=abs_change>=0)) + ylab('') + xlab('Absolute change') +
    geom_bar(stat = 'identity') + theme_linedraw() + theme(legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_fill_manual(values = c('#E04D55', '#15356B'))

  sp2 = ggarrange(p2, p3, heights = c(2, 2), ncol = 1, labels = c('B','C'))

  p = ggarrange(sp1, sp2, ncol = 2, nrow = 1, widths = c(0.6, 0.4), labels = c('A',''))
  ggsave(p, filename = 'plots/Fig_4_Abundance_changes.pdf', width = 9, height = 7)}


MainK <- function(){
  loaded_data = LoadModelOutput()
  changes_tab = loaded_data$changes
  taxonomy = loaded_data$taxonomy
  tree = LoadTree()

  changes_tab = changes_tab[changes_tab$MAG %in% tree$tip.label,]
  tree = keep.tip(tree, tree$tip.label[tree$tip.label %in% changes_tab$MAG])

  changes_tab = changes_tab %>% arrange(match(MAG, tree$tip.label))

  # Correlation between scenarios
  sink(file = 'stats/strain_scenarios_correlations.txt')
  print('Correlation between RCP2.6 and 4.5: ')
  print(cor.test(changes_tab$median_future[changes_tab$scenario == 370], changes_tab$median_future[changes_tab$scenario == 126], method = 'pearson'))
  print('Correlation between RCP8.5 and 4.5: ')
  print(cor.test(changes_tab$median_future[changes_tab$scenario == 370], changes_tab$median_future[changes_tab$scenario == 585], method = 'pearson'))
  sink(file = NULL)
  
  PlotModelPerformances(changes_tab)
  CompareScenariosChanges(changes_tab)
  ProportionIncrease(changes_tab)
  PhyloSignal(changes_tab, tree)
  MonophyleticClades(tree, changes_tab, taxonomy)
  PlotTreeChanges(changes_tab, tree, taxonomy)
  }
