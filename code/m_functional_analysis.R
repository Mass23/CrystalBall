library(ggplot2)
library(ggpubr)
library(tidyverse)
library(reshape2)
library(matrixStats)
library(phytools)
library(ape)
library(ggtree)
library(ggtreeExtra)
library(forcats)
library(viridis)
library(ggnewscale)
library(ranger)
library(pROC)
library(doMC)
library(poolr)
library(phytools)

LoadData <- function(){
  model_out = read.csv('stats/model_res_all_scenarios_final.csv')
  model_out = model_out %>% filter(scenario == 370, r2 >= 0.1)
  model_out$Category = 'Not significant'
  model_out$Category[(model_out$wilcox_p < 0.05) & (model_out$median_change < 0)] = 'Decrease'
  model_out$Category[(model_out$wilcox_p < 0.05) & (model_out$median_change > 0)] = 'Increase'

  # Functions to table
  func_mags = read.csv('data/raw/bacteria/NOMIS_MAGs_Contigs_Genes_KO.txt', sep='\t')
  func_mags = func_mags %>% filter(MAGs %in% model_out$MAG)
  func_mags$value = 1
  func_tab = dcast(func_mags, formula = MAGs ~ KEGG_ko, value.var = 'value', fill = 0)
  func_tab$Category = map_chr(func_tab$MAGs, function(x) unique(model_out$Category[model_out$MAG == x]))

  # Phylogeny
  tree = read.tree('data/raw/bacteria/treeBacteria.tree')
  tree = midpoint.root(tree)
  tree$tip.label = vapply(tree$tip.label, function(x) strsplit(x, '.fa')[[1]][1], FUN.VALUE = character(1))

  tree = keep.tip(tree, unique(func_tab$MAGs[func_tab$MAGs %in% tree$tip.label]))
  func_tab = func_tab[func_tab$MAGs %in% tree$tip.label,]

  phylo_clusters = cutree(hclust(as.dist(cophenetic(tree)), method = 'ward.D2'), 20)
  func_tab$PhyloCluster = map_chr(func_tab$MAGs, function(x) as.character(phylo_clusters[names(phylo_clusters) == x]))
  func_tab$ClusterWeight = map_dbl(func_tab$PhyloCluster, function(x) 1/sum(func_tab$PhyloCluster == x))

  # Completeness
  checkm2_out = read.delim('data/raw/bacteria/checkm2_final.tsv')
  func_tab$Completeness = map_dbl(func_tab$MAGs, function(x) checkm2_out$Completeness[checkm2_out$Name == x])

  # Final weight is: completeness * cluster weight (smaller clusters weight as much as larger ones) * class_weight
  func_tab$ClassWeight = map_dbl(func_tab$Category, function(x) 1/(sum(func_tab$Category == x)))
  func_tab$Weight = func_tab$ClassWeight * func_tab$ClusterWeight * func_tab$Completeness

  func_tab$Category = ifelse(func_tab$Category == 'Decrease', 0, 1)
  func_tab = func_tab %>% filter(PhyloCluster %in% names(which(rowSums(table(func_tab$PhyloCluster, func_tab$Category) >= 3) == 2))) # 34 clusters kept
  return(func_tab)}

FunctionalRandomForests <- function(func_tab, n_cores=8){
  # Leave one cluster out cross validation, computing feature importance on hold out set.
  registerDoMC(n_cores)
  rf_out = foreach(clu = unique(func_tab$PhyloCluster), .combine = rbind) %dopar% {
    func_tab$Category = as.factor(func_tab$Category)
    train = func_tab %>% filter(PhyloCluster != clu) %>% select(-ClassWeight, -Completeness, -ClusterWeight, -PhyloCluster, -MAGs)
    test = func_tab %>% filter(PhyloCluster == clu) %>% select(-ClassWeight, -Completeness, -ClusterWeight, -PhyloCluster, -MAGs)
    
    train_stats = data.frame()
    for (i in 1:50){
      reg_factor = sample(c(0.2,0.4,0.6,0.8,1),1)[1]
      reg_depth = sample(c(TRUE, FALSE), 1)[1]
      n_trees = sample(c(250,500,750,1000,1500), 1)[1]
      mtr = sample(c(20,30,40,50), 1)[1]
      max_d = sample(c(5,15,30,50), 1)[1]
      s_frac = sample(c(0.5,0.6,0.7,0.8))[1]
      splitr = sample(c('hellinger','extratrees','gini'), 1)[1]
      rf_mod = ranger(y = train$Category, x = train %>% select(-Weight, -Category), importance = "none", case.weights = train$Weight, 
                      num.trees = n_trees, max.depth = max_d, mtry = mtr, sample.fraction = s_frac, splitrule = splitr,
                      regularization.factor = reg_factor, regularization.usedepth = reg_depth)
      train_stats = rbind(train_stats, data.frame(reg_factor=reg_factor, reg_depth=reg_depth, mtr=mtr, splitr=splitr,
                                                  n_trees=n_trees, max_d=max_d, s_frac=s_frac, splitr=splitr, oob_err=rf_mod$prediction.error))}
    best_model = train_stats %>% arrange(oob_err) %>% head(1)
    
    rf_final = ranger(y = train$Category, x = train %>% select(-Weight, -Category), importance = "permutation", case.weights = train$Weight, 
                      num.trees = best_model$n_trees[1], max.depth = best_model$max_d[1], mtry = best_model$mtr[1], sample.fraction = best_model$s_frac[1], 
                      regularization.factor = best_model$reg_factor[1], regularization.usedepth = best_model$reg_depth[1], splitrule = best_model$splitr[1])
    imp_pvals = as.data.frame(importance_pvalues(rf_final, data = test %>% select(-Weight), formula = Category ~ ., method = "altmann", num.permutations = 999))
    
    imp_pvals$KEGG_ko = train %>% select(-Weight, -Category) %>% colnames()
    imp_pvals$cluster = clu
    return(imp_pvals)}
  write.csv(rf_out, file = 'data/processed/functional_random_forests.csv', quote = F, row.names = F)
  return(rf_out)}

GetTopKOs <- function(rf_out){
  kegg_tab = read.csv('data/raw/bacteria/keggPathwayGood.txt', sep='\t')

  sign_kos = data.frame()
  for (clu in unique(rf_out$cluster)){
    rf_clu = rf_out %>% filter(cluster == clu)
    top_kos = rf_clu %>% filter(pvalue < 0.01, importance >= quantile(rf_clu$importance, probs = 0.95)) %>% pull(KEGG_ko)
    sign_kos = rbind(sign_kos, data.frame(Cluster=rep(clu, length(top_kos)), KEGG_ko=as.vector(top_kos)))}

  top_kos = sign_kos %>% group_by(KEGG_ko) %>% summarise(n=n()) %>% filter(n > 15) %>% arrange(-n)
  write.csv(top_kos, 'stats/functional_top_kos.csv', quote=F, row.names=F)
  return(sign_kos)}
#   KEGG_ko       n
#     <chr>    <int>
#  1 K14215      19 - trans,polycis-decaprenyl diphosphate synthase - biosynthesis of cell wall in Mycobacterium (https://journals.asm.org/doi/10.1128/JB.186.22.7564-7570.2004)
#  2 K18455      19 - mycothiol S-conjugate amidase - Detoxification, stress related (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0115075)
#  3 K00697      18 - trehalose 6-phosphate synthase - sugar metabolism and related to cold stress (https://www.tandfonline.com/doi/full/10.1080/21505594.2020.1809326)
#  4 K21273      18 - trans,polycis-polyprenyl diphosphate synthase - biosynthesis of cell wall in Mycobacterium (https://journals.asm.org/doi/10.1128/JB.186.22.7564-7570.2004)
#  5 K03588      17 - cell division protein FtsW - transporter of lipid-linked cell wall precursor across the membrane (https://www.embopress.org/doi/full/10.1038/emboj.2011.61)
#  6 K11263      17 - acetyl-CoA/propionyl-CoA carboxylase, biotin carboxylase, biotin carboxyl carrier protein - carbon storage (https://journals.asm.org/doi/10.1128/AEM.03167-14)
#  7 K15327      17 - polyketide biosynthesis malonyl-CoA-[acyl-carrier-protein] transacylase - fatty acids, phospholipids (https://pubmed.ncbi.nlm.nih.gov/18384517/)
#  8 K01716      16 - 3-hydroxyacyl-[acyl-carrier protein] dehydratase / trans-2-decenoyl-[acyl-carrier protein] isomerase - fatty acids (https://www.sciencedirect.com/science/article/pii/S0021925818688301?via%3Dihub)
#  9 K03071      16 - preprotein translocase subunit SecB - chaperone linked to temperature sensitivity of growth (https://www.pnas.org/doi/full/10.1073/pnas.0402398101)
#  10 K03666     16 - host factor-I protein - biofilm formation (https://journals.asm.org/doi/10.1128/JB.183.6.1997-2005.2001)
#  11 K07320     16 - ribosomal protein L3 glutamine methyltransferase - transcription 
# polyketides and fatty acids, related to cold adaptation of fatty acids! : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4095770/

FunctionalEnrichment <- function(func_tab, sign_kos){
  kegg_tab = read.csv('data/raw/bacteria/keggPathwayGood.txt', sep='\t')
  keggs_to_keep = unique(func_tab$KEGG_ko)
  kegg_tab_clean = kegg_tab %>% filter(ko %in% keggs_to_keep)
  kegg_tab_clean$Sign = ifelse(kegg_tab_clean$ko %in% sign_kos$KEGG_ko, 1, 0)

  enrich_analysis = data.frame()
  for (cat in unique(kegg_tab_clean$Path_2)){
    if ((sum(kegg_tab_clean$Sign) >= 0) & (sum(kegg_tab_clean$Path_2 == cat) >= 5)){
      cont_tab = table(kegg_tab_clean$Path_2 == cat, kegg_tab_clean$Sign == 1)
      ftest = fisher.test(cont_tab)
      enrich_analysis = rbind(enrich_analysis, data.frame(Path_2 = cat, p=ftest$p.value, OR=ftest$estimate))}}
  enrich_analysis$padj_holm = p.adjust(enrich_analysis$p, method = 'holm')
  enrich_analysis %>% filter(padj_holm < 0.05, OR > 1)
  #                        Path_2             p         OR     padj_holm
  # 09101 Carbohydrate metabolism  3.934662e-10  2.1794938  7.082392e-09
  #       09102 Energy metabolism  9.376012e-05  1.8888975  1.500162e-03
  write.csv(enrich_analysis %>% filter(padj_holm < 0.05, OR > 1), 'stats/functional_enrichment.csv', quote=F, row.names=F)

  sign_categories = enrich_analysis %>% filter(padj_holm< 0.05, OR > 1) %>% pull(Path_2)
  sign_categories_df = expand.grid(MAG=func_tab$MAGs, Category=sign_categories)

  mags_summary = read.csv('Data/genomeInformation.csv')

  sign_categories_df$Group = map_chr(sign_categories_df$MAG, function(x) ifelse(func_tab$Category[func_tab$MAGs == x] == 1, 'Decrease', 'Not significant or increase'))
  sign_categories_df$Count = map_int(1:nrow(sign_categories_df), function(i) sum(func_tab[func_tab$MAGs == sign_categories_df$MAG[i], colnames(func_tab) %in% kegg_tab$ko[kegg_tab$Path_2 == sign_categories_df$Category[i]]]))
  sign_categories_df$Count_total = map_int(1:nrow(sign_categories_df), function(i) sum(func_tab[func_tab$MAGs == sign_categories_df$MAG[i], startsWith(colnames(func_tab), 'K')]))
  sign_categories_df$Completeness = map_dbl(sign_categories_df$MAG, function(x) checkm2_out$Completeness[checkm2_out$Name == x]/100)
  sign_categories_df$Contamination = map_dbl(sign_categories_df$MAG, function(x) checkm2_out$Contamination[checkm2_out$Name == x]/100)
  sign_categories_df$Size = map_dbl(sign_categories_df$MAG, function(x) mags_summary$length[mags_summary$genome == paste0(x, '.fa', collapse = '')])

  sign_categories_df$Category = as.character(sign_categories_df$Category)
  sign_categories_df$Category[sign_categories_df$Category == '09101 Carbohydrate metabolism'] = 'Carbohydrate metabolism'
  sign_categories_df$Category[sign_categories_df$Category == '09102 Energy metabolism'] = 'Energy metabolism'

  p1 = sign_categories_df %>% select(Completeness, Group, MAG, Count, Category) %>% distinct() %>%
    ggplot(aes(x=Completeness*100, y=Count, colour=Group)) + geom_point(alpha=0.1) + geom_smooth(method = 'lm', formula = y ~ x -1) + facet_grid(~Category) + 
    theme_linedraw() + ylab('KO [#]') + xlab('Completeness [%]') + scale_colour_manual(values=c('#E04D55', '#15356B')) + theme(panel.grid = element_blank())

  p2 = sign_categories_df %>% select(Completeness, Group, MAG, Count_total) %>% distinct() %>%
    ggplot(aes(x=Group, y=Count_total*(1/Completeness), colour=Group)) + geom_boxplot() + theme_linedraw() + stat_compare_means(label.y = 640) +
    ylab('Total KO [#]') + xlab('') + scale_colour_manual(values=c('#E04D55', '#15356B')) + ylim(100, 660) +
    theme(panel.grid = element_blank(), legend.position = 'none', axis.text.x = element_blank(), axis.ticks.x = element_blank())

  p3 = sign_categories_df %>% select(Completeness, Group, MAG, Size) %>% distinct() %>%
    ggplot(aes(x=Group, y=Size/1000000*(1/Completeness), colour=Group)) + geom_boxplot() + theme_linedraw() + stat_compare_means(label.y = 11.6) +
    ylab('Genome size [mbp]') + xlab('') + scale_colour_manual(values=c('#E04D55', '#15356B')) + ylim(1, 12) +
    theme(panel.grid = element_blank(), legend.position = 'none', axis.text.x = element_blank(), axis.ticks.x = element_blank())

  pb = ggarrange(p2, p3, ncol = 2, nrow = 1, labels = c('B', 'C'))
  p = ggarrange(p1, pb, ncol = 1, nrow = 2, align = 'v', common.legend = T, legend = 'bottom', labels = c('A','',''))
  ggsave('Plots/Fig_6_functional.pdf', width = 5, height = 5.5)
  return(sign_categories_df)}

FunctionalCategories <- function(sign_categories_df){
  model_carb_met = lm(data = sign_categories_df[sign_categories_df$Category == 'Carbohydrate metabolism',], formula = Count ~ Completeness:Group -1)
  print(summary(model_carb_met))
  #                                               Estimate Std. Error t value Pr(>|t|)    
  # Completeness:GroupDecrease                    0.511031   0.005193   98.41   <2e-16 ***
  # Completeness:GroupNot significant or increase 0.569076   0.008317   68.43   <2e-16 ***
  # Residual standard error: 15.28 on 1754 degrees of freedom
  # Multiple R-squared:  0.8912,	Adjusted R-squared:  0.8911 
  # F-statistic:  7184 on 2 and 1754 DF,  p-value: < 2.2e-16
  print(predict(model_carb_met, newdata = data.frame(Group=c('Decrease','Not significant or increase'), Completeness=c(1,1)), se.fit = T))
  # Fit : Decrease = 51.10307 +/- 0.5192785, Increase and NS = 56.90761 +/- 0.8316619

  model_en_met = lm(data = sign_categories_df[sign_categories_df$Category == 'Energy metabolism',], formula = Count ~ Completeness:Group -1)
  print(summary(model_en_met))
  #                                               Estimate Std. Error t value Pr(>|t|)    
  # Completeness:GroupDecrease                    0.363233   0.004334   83.82   <2e-16 ***
  # Completeness:GroupNot significant or increase 0.393663   0.006940   56.72   <2e-16 ***
  # Residual standard error: 12.75 on 1754 degrees of freedom
  # Multiple R-squared:  0.8538,	Adjusted R-squared:  0.8536 
  # F-statistic:  5121 on 2 and 1754 DF,  p-value: < 2.2e-16
  print(predict(model_en_met, newdata = data.frame(Group=c('Decrease','Not significant or increase'), Completeness=c(1,1)), se.fit = T))
  # Fit : Decrease = 36.32331 +/- 0.4333552, Increase and NS = 39.36633 +/- 0.6940496
  }

FunctionalGenomeStats <- function(sign_categories_df){
  print('KO number:')
  print(sign_categories_df %>% select(Completeness, Group, MAG, Count_total) %>% distinct() %>% group_by(Group) %>% summarise(median = median(Count_total*(1/Completeness))))
  #   Group                       median
  # 1 Decrease                      320.
  # 2 Not significant or increase   335.
  # difference in median: 15
  print(sign_categories_df %>% select(Completeness, Group, MAG, Count_total) %>% distinct() %>% summarise(wp = wilcox.test(Count_total*(1/Completeness) ~ Group)$p.value))
  # 4.292379e-07

  print('Genome length:')
  print(sign_categories_df %>% select(Completeness, Group, MAG, Size) %>% distinct() %>% group_by(Group) %>% summarise(median = median(Size*(1/Completeness))))
  # 1 Decrease                    3809740.
  # 2 Not significant or increase 4187851.
  # differece in median: 378111
  print(sign_categories_df %>% select(Completeness, Group, MAG, Size) %>% distinct() %>% summarise(wp = wilcox.test(Size*(1/Completeness) ~ Group)$p.value))
  # 2.244707e-07

  # KO per mbp decrease
  print('KO / mbp (decrease):')
  print(median(sign_categories_df %>% filter(Group == 'Decrease') %>% mutate(density = Count_total/Size*1000000) %>% pull(density)))
  # 84.33736

  # KO per mbp increase/not sign
  print('KO / mbp (increase):')
  print(median(sign_categories_df %>% filter(Group == 'Not significant or increase') %>% mutate(density = Count_total/Size*1000000) %>% pull(density)))
  # 81.64844
  }

MainM <- function(){
  func_tab = LoadData()
  rf_res = FunctionalRandomForests(func_tab)
  sign_kos = GetTopKOs(rf_res)
  sign_categories_df = FunctionalEnrichment(func_tab, sign_kos)
  FunctionalCategories(sign_categories_df)
  FunctionalGenomeStats(sign_categories_df)
  }


