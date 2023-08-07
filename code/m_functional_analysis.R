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
#install.packages('poolr',repos = "http://cran.us.r-project.org")
library(poolr)
library(phytools)
library(mgcv)
library(performance)

LoadData <- function(){
  model_out = read.csv('stats/model_res_all_scenarios_final.csv')
  model_out = model_out %>% select(MAG, wilcox_p, median_change, r2, scenario, mean_rel_ab) %>% distinct()
  decreasing_mags = model_out %>% filter(wilcox_p < 0.05, median_change < 0) %>% group_by(MAG) %>% summarise(n=n()) %>% filter(n > 2) %>% pull(MAG)
  model_out = model_out %>% filter(scenario == 370, r2 > 0.1)
  print(paste0('Number of MAGs kept: ',nrow(model_out)))
  model_out$Category = 'Others'
  model_out$Category[model_out$MAG %in% decreasing_mags] = 'Decrease'

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

  phylo_clusters = cutree(hclust(as.dist(cophenetic(tree)), method = 'ward.D2'), 10)
  func_tab$PhyloCluster = map_chr(func_tab$MAGs, function(x) as.character(phylo_clusters[names(phylo_clusters) == x]))
  func_tab$ClusterWeight = map_dbl(func_tab$PhyloCluster, function(x) 1/sum(func_tab$PhyloCluster == x))

  print(table(func_tab$Category, func_tab$PhyloCluster))
  print(sum(table(func_tab$Category, func_tab$PhyloCluster)))

  # Completeness
  checkm2_out = read.delim('data/raw/bacteria/checkm2_final.tsv')
  func_tab$Completeness = map_dbl(func_tab$MAGs, function(x) checkm2_out$Completeness[checkm2_out$Name == x])

  # Final weight is: completeness * cluster weight (smaller clusters weight as much as larger ones) * class_weight
  func_tab$ClassWeight = map_dbl(func_tab$Category, function(x) 1/(sum(func_tab$Category == x)))
  func_tab$Weight = func_tab$ClassWeight * func_tab$ClusterWeight * func_tab$Completeness

  func_tab$Category = ifelse(func_tab$Category == 'Decrease', 1, 0)
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
      set.seed(i)
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
    top_kos = rf_clu %>% filter(pvalue < 0.05, importance >= quantile(rf_clu$importance, probs = 0.95)) %>% pull(KEGG_ko)
    sign_kos = rbind(sign_kos, data.frame(Cluster=rep(clu, length(top_kos)), KEGG_ko=as.vector(top_kos)))}

  top_kos = sign_kos %>% group_by(KEGG_ko) %>% summarise(n=n()) %>% filter(n > 5) %>% arrange(-n)
  write.csv(top_kos, 'stats/functional_top_kos.csv', quote=F, row.names=F)
  return(sign_kos)}
#    KEGG_ko     n
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
  keggs_to_keep = unique(colnames(func_tab))[startsWith(unique(colnames(func_tab)), 'K')]
  kegg_tab_clean = kegg_tab %>% filter(ko %in% keggs_to_keep) # Perform the enrichment only for the KOs in the dataset
  kegg_tab_clean$Sign = ifelse(kegg_tab_clean$ko %in% sign_kos$KEGG_ko, 1, 0)

  enrich_analysis = data.frame()
  for (cat in unique(kegg_tab_clean$Path_2)){
    if ((sum(kegg_tab_clean$Sign == 0) >= 1) & (sum(kegg_tab_clean$Sign == 1) >= 1) & (sum(kegg_tab_clean$Path_2 == cat) >= 5)){
      cont_tab = table(kegg_tab_clean$Path_2 == cat, kegg_tab_clean$Sign == 1)
      ftest = fisher.test(cont_tab)
      enrich_analysis = rbind(enrich_analysis, data.frame(Path_2 = cat, p=ftest$p.value, OR=ftest$estimate))}}
  enrich_analysis$padj = p.adjust(enrich_analysis$p, method = 'bonferroni')
  enrich_analysis %>% filter(padj < 0.05, OR > 1)
  write.csv(enrich_analysis, 'stats/functional_enrichment.csv', quote=F, row.names=F)}

FunctionalGenomeStats <- function(func_tab){
  enrich_analysis = read.csv('stats/functional_enrichment.csv')
  sign_categories = enrich_analysis %>% filter(padj < 0.05, OR > 1) %>% pull(Path_2)
  sign_categories_df = expand.grid(MAG=func_tab$MAGs, Category=sign_categories)

  mags_summary = read.csv('data/raw/bacteria/genomeInformation.csv')
  mags_summary$genome = map_chr(mags_summary$genome, function(x) gsub('.fa','',x))
  checkm2_out = read.delim('data/raw/bacteria/checkm2_final.tsv')
  kegg_tab = read.csv('data/raw/bacteria/keggPathwayGood.txt', sep='\t')
  spec_tab = read.csv('data/raw/bacteria/spec_gen.tsv', sep='\t')

  model_out = read.csv('stats/model_res_all_scenarios_final.csv')
  model_out = model_out %>% select(MAG, wilcox_p, median_change, r2, scenario, mean_rel_ab) %>% distinct()
  decreasing_mags = model_out %>% filter(wilcox_p < 0.05, median_change < 0) %>% group_by(MAG) %>% summarise(n=n()) %>% filter(n > 2) %>% pull(MAG)

  model_out$Category = 'Others'
  model_out$Category[model_out$MAG %in% decreasing_mags] = 'Decrease'

  sign_categories_df$Group = map_chr(sign_categories_df$MAG, function(x) unique(model_out$Category[model_out$MAG == x]))
  sign_categories_df$PhyloCluster =  map_chr(sign_categories_df$MAG, function(x) paste0('PC', unique(func_tab$PhyloCluster[func_tab$MAGs == x])))
  sign_categories_df$Count = map_int(1:nrow(sign_categories_df), function(i) sum(func_tab[func_tab$MAGs == sign_categories_df$MAG[i], colnames(func_tab) %in% kegg_tab$ko[kegg_tab$Path_2 == sign_categories_df$Category[i]]]))
  sign_categories_df$Count_total = map_int(1:nrow(sign_categories_df), function(i) sum(func_tab[func_tab$MAGs == sign_categories_df$MAG[i], startsWith(colnames(func_tab), 'K')]))
  sign_categories_df$Count_unique = map_int(1:nrow(sign_categories_df), function(i) sum(func_tab[func_tab$MAGs == sign_categories_df$MAG[i], startsWith(colnames(func_tab), 'K')] > 0))
  sign_categories_df$Count_prop = sign_categories_df$Count / sign_categories_df$Count_total
  sign_categories_df$Genome_length = map_int(sign_categories_df$MAG, function(x) mags_summary$length[mags_summary$genome == x])
  sign_categories_df$N50 = map_int(sign_categories_df$MAG, function(x) mags_summary$N50[mags_summary$genome == x])
  sign_categories_df$MeanRelAb = map_dbl(sign_categories_df$MAG, function(x) unique(model_out$mean_rel_ab[model_out$MAG == x]))

  colnames(checkm2_out)[colnames(checkm2_out) == 'Name'] = 'MAG'
  checkm2_out = checkm2_out %>% filter(MAG %in% sign_categories_df$MAG)
  checkm2_out$Completeness = checkm2_out$Completeness/100
  checkm2_out$Contamination = checkm2_out$Contamination/100
  sign_categories_df = inner_join(sign_categories_df, checkm2_out, by = 'MAG')

  sign_categories_df$Category = as.character(sign_categories_df$Category)
  sign_categories_df$Category[sign_categories_df$Category == '09101 Carbohydrate metabolism'] = 'Carbohydrate metabolism'
  sign_categories_df$Category[sign_categories_df$Category == '09102 Energy metabolism'] = 'Energy metabolism'

  write.csv(sign_categories_df, 'stats/functional_stats.csv', quote=F, row.names=F)}

FitModel <- function(mod_data, resp_var){
  
  form = paste0(resp_var, " ~ Group + s(Completeness, k=-1) + s(Contamination, k=-1)", collapse='')
  gaussian = bam(data = mod_data, formula = eval(parse(text=form)), family=gaussian(), weights=pc_weight*Completeness*MeanRelAb)
  return(gaussian)}

FunctionalModels <- function(){
  sign_categories_df = read.csv('stats/functional_stats.csv') %>% na.omit()
  sign_categories_genomes = sign_categories_df %>% select(-Category, -Count, -Count_prop) %>% distinct()
  sign_categories_df$pc_weight = map_dbl(sign_categories_df$PhyloCluster, function(x) 1/sum(sign_categories_genomes$PhyloCluster == x))
  sign_categories_genomes$pc_weight = map_dbl(sign_categories_genomes$PhyloCluster, function(x) 1/sum(sign_categories_genomes$PhyloCluster == x))
  sign_categories_genomes$Genome_length = sign_categories_genomes$Genome_length / 1000000

  sign_categories_df$Group <- factor(sign_categories_df$Group, levels = c('Others', 'Decrease'))
  sign_categories_genomes$Group <- factor(sign_categories_genomes$Group, levels = c('Others', 'Decrease'))

  print(table(sign_categories_genomes$Group))
  sink('stats/functional_bulk_features.txt')

  print('-----------------------------------------------------------------------------------------------')
  print('Carbohydrate metabolism:')
  cm_count = bam(data = sign_categories_df %>% filter(Category == 'Carbohydrate metabolism'), 
                                            formula = Count ~ Group + s(Completeness, k=-1, bs='ts') + s(Contamination, k=3, bs='ts'), 
                                            family=gaussian(), weights=pc_weight*Completeness*MeanRelAb)
  cm_prop = bam(data = sign_categories_df %>% filter(Category == 'Carbohydrate metabolism'), 
                                            formula = Count_prop ~ Group + s(Completeness, k=-1, bs='ts') + s(Contamination, k=3, bs='ts'), 
                                            family=quasibinomial(), weights=pc_weight*Completeness*MeanRelAb)
  print(summary(cm_count))
  print(summary(cm_prop))

  print('-----------------------------------------------------------------------------------------------')
  print('Energy metabolism:')
  em_count = bam(data = sign_categories_df %>% filter(Category == 'Energy metabolism'), 
                                            formula = Count ~ Group + te(Completeness, Contamination, N50, k=3, bs='cs'),
                                            family=gaussian(), weights=pc_weight*Completeness*MeanRelAb)
  em_prop = bam(data = sign_categories_df %>% filter(Category == 'Energy metabolism'), 
                                            formula = Count_prop ~ Group + te(Completeness, Contamination, N50, k=3, bs='cs'),
                                            family=quasibinomial(), weights=pc_weight*Completeness*MeanRelAb)
  print(summary(em_count))
  print(summary(em_prop))

  print('-----------------------------------------------------------------------------------------------')
  print('KO number:')
  ko_count = bam(data = sign_categories_genomes, 
                 formula = Count_total ~ Group + te(Completeness, Contamination, N50, k=3, bs='cs'), 
                 family=gaussian(), weights=pc_weight*Completeness*MeanRelAb)
  print(summary(ko_count))

  print('-----------------------------------------------------------------------------------------------')
  print('Genome length:')
  len_count = bam(data = sign_categories_genomes, 
                 formula = Genome_length ~ Group + te(Completeness, Contamination, N50, k=3, bs='cs'),
                 family=gaussian(), weights=pc_weight*Completeness*MeanRelAb)
  print(summary(len_count))

  print('-----------------------------------------------------------------------------------------------')
  print('Count_unique / Count_total:')
  len_count = bam(data = sign_categories_genomes, 
                 formula = Count_total/Count_unique ~ Group + te(Completeness, Contamination, N50, k=3, bs='cs'),
                 family=gaussian(), weights=pc_weight*Completeness*MeanRelAb)
  print(summary(len_count))

  print('-----------------------------------------------------------------------------------------------')
  print('Comparison:')
  comp_len_count = lm(data = sign_categories_genomes, formula = Count_total ~ Group:log(Genome_length), 
                      weights=pc_weight*Completeness*MeanRelAb)
  print(summary(comp_len_count))
  print(anova(comp_len_count))

  print('Comparison with phylo. structure:')
  comp_len_count_phylo = lm(data = sign_categories_genomes, formula = Count_total ~ PhyloCluster + Group:log(Genome_length), 
                            weights=pc_weight*Completeness*MeanRelAb)
  print(summary(comp_len_count_phylo))
  print(anova(comp_len_count_phylo))
  sink()}


FunctionalPlots <- function(){
  sign_categories_df = read.csv('stats/functional_stats.csv') %>% na.omit()
  sign_categories_genomes = sign_categories_df %>% select(-Category, -Count, -Count_prop) %>% distinct()
  sign_categories_df$pc_weight = map_dbl(sign_categories_df$PhyloCluster, function(x) 1/sum(sign_categories_genomes$PhyloCluster == x))
  sign_categories_genomes$pc_weight = map_dbl(sign_categories_genomes$PhyloCluster, function(x) 1/sum(sign_categories_genomes$PhyloCluster == x))
  sign_categories_genomes$Genome_length = sign_categories_genomes$Genome_length / 1000000

  sign_categories_df$Group <- factor(sign_categories_df$Group, levels = c('Others', 'Decrease'))
  sign_categories_genomes$Group <- factor(sign_categories_genomes$Group, levels = c('Others', 'Decrease'))

  p1 = sign_categories_df %>% select(Completeness, Contamination, pc_weight, MeanRelAb, Group, MAG, Count, Count_total, Category) %>% distinct() %>% filter(Category %in% c('Energy metabolism', 'Carbohydrate metabolism')) %>%
    ggplot(aes(x=Group, y=(Count/Count_total)*100, colour=Group)) + geom_boxplot(aes(weight=Completeness*pc_weight*MeanRelAb)) + facet_grid(~Category) + 
    theme_linedraw() + xlab('') + ylab('KO proportion [%]') + scale_colour_manual(values=c('#E04D55', '#15356B')) + theme(panel.grid = element_blank())

  p2 = sign_categories_genomes %>% select(Completeness, Contamination, pc_weight, MeanRelAb, Group, MAG, Count_total, Genome_length) %>% 
    ggplot(aes(x=log(Genome_length), y=Count_total, colour=Group, group=Group)) + geom_point(alpha=0.1) + geom_smooth(method='lm', formula = y ~ x, aes(weight=Completeness*pc_weight*MeanRelAb)) + 
    theme_linedraw() + ylab('Total KO [#]') + xlab('Genome size [ln mbp]') + scale_colour_manual(values=c('#E04D55', '#15356B')) + 
    theme(panel.grid = element_blank(), legend.position = 'none')

  p = ggarrange(p1, p2, ncol = 1, nrow = 2, align = 'v', common.legend = T, legend = 'bottom', labels = c('A','B'))
  ggsave('plots/Fig_6_functional.pdf', width = 3.5, height = 5.5)
  }

MainM <- function(){
  #func_tab = LoadData()
  #rf_res = FunctionalRandomForests(func_tab)
  #rf_res = read.csv('data/processed/functional_random_forests.csv')
  #sign_kos = GetTopKOs(rf_res)
  #sign_categories_df = FunctionalEnrichment(func_tab, sign_kos)
  #FunctionalGenomeStats(func_tab)
  FunctionalModels()
  FunctionalPlots()
  }


