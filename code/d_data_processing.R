library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(reshape2)
library(matrixStats)
library(phytools)
library(ape)
library(forcats)
library(ggnewscale)
library(vegan)
library(purrr)
library(picante)
library(maps)
library(mapproj)

BioclimPCA <- function(all_data) {
  # Bioclimatic variables PCA
  bioclim_data = all_data %>% select(c(
    'Sample',
    'Mountain_range',
    'Scenario',
    'Date',
    "clim_bio1",
    "clim_bio10",
    "clim_bio11",
    "clim_bio12",
    "clim_bio13",
    "clim_bio14",
    "clim_bio15",
    "clim_bio16",
    "clim_bio17",
    "clim_bio18",
    "clim_bio19",
    "clim_bio2",
    "clim_bio3",
    "clim_bio4",
    "clim_bio5",
    "clim_bio6",
    "clim_bio7",
    "clim_bio8",
    "clim_bio9"
  ))
  
  bioclim_pca = prcomp(
    bioclim_data %>% select(-Sample,-Mountain_range,-Scenario,-Date),
    center = T,
    scale. = T
  )
  summary(bioclim_pca) #6 variables = 96.447% variance explained
  
  bioclim_pca_data = as.data.frame(bioclim_pca$x[, 1:6])
  bioclim_pca_data$Sample = bioclim_data$Sample
  bioclim_pca_data$Mountain_range = bioclim_data$Mountain_range
  bioclim_pca_data$Scenario = bioclim_data$Scenario
  bioclim_pca_data$Date = bioclim_data$Date
  
  bioclim_pca_variables = as.data.frame(bioclim_pca$rotation[, 1:6])
  bioclim_pca_variables$Variable = rownames(bioclim_pca_variables)
  
  for (i in 1:nrow(all_data)) {
    all_data$bioclim_PC1[i] = bioclim_pca_data$PC1[(bioclim_pca_data$Sample == all_data$Sample[i]) &
                                                     (bioclim_pca_data$Scenario == all_data$Scenario[i]) &
                                                     (bioclim_pca_data$Date == all_data$Date[i])]
    all_data$bioclim_PC2[i] = bioclim_pca_data$PC2[(bioclim_pca_data$Sample == all_data$Sample[i]) &
                                                     (bioclim_pca_data$Scenario == all_data$Scenario[i]) &
                                                     (bioclim_pca_data$Date == all_data$Date[i])]
    all_data$bioclim_PC3[i] = bioclim_pca_data$PC3[(bioclim_pca_data$Sample == all_data$Sample[i]) &
                                                     (bioclim_pca_data$Scenario == all_data$Scenario[i]) &
                                                     (bioclim_pca_data$Date == all_data$Date[i])]
    all_data$bioclim_PC4[i] = bioclim_pca_data$PC4[(bioclim_pca_data$Sample == all_data$Sample[i]) &
                                                     (bioclim_pca_data$Scenario == all_data$Scenario[i]) &
                                                     (bioclim_pca_data$Date == all_data$Date[i])]
    all_data$bioclim_PC5[i] = bioclim_pca_data$PC5[(bioclim_pca_data$Sample == all_data$Sample[i]) &
                                                     (bioclim_pca_data$Scenario == all_data$Scenario[i]) &
                                                     (bioclim_pca_data$Date == all_data$Date[i])]
    all_data$bioclim_PC6[i] = bioclim_pca_data$PC6[(bioclim_pca_data$Sample == all_data$Sample[i]) &
                                                     (bioclim_pca_data$Scenario == all_data$Scenario[i]) &
                                                     (bioclim_pca_data$Date == all_data$Date[i])]
  }
  return(list(all_data=all_data, pca_vars=bioclim_pca_variables))
}

PlotBioclimPCA <- function(all_data_with_bioclim_pca, pca_variables) {
  p1 = ggplot() + geom_point(
    all_data_with_bioclim_pca,
    mapping = aes(x = bioclim_PC1, y = bioclim_PC2, colour = Mountain_range)
  ) + xlab('PC1') + ylab('PC2') +
    geom_segment(
      pca_variables,
      mapping = aes(
        x = 0,
        xend = PC1 * 7,
        y = 0,
        yend = PC2 * 7
      )
    ) +
    geom_text(
      pca_variables,
      mapping = aes(
        x = PC1 * 7,
        y = PC2 * 1.2 * 7,
        label = Variable
      )
    ) + theme_minimal() +
    scale_colour_manual(
      values = c(
        '#CC2F50',
        '#1A5A61',
        '#E3C78D',
        '#4BC992',
        '#806A2A',
        '#177FCF',
        '#982DA1',
        '#A7BFE8',
        'dimgrey',
        'tomato',
        'grey'
      ),
      name = 'Mountain range'
    )
  leg = get_legend(p1)
  
  p2 = ggplot() + geom_point(
    all_data_with_bioclim_pca,
    mapping = aes(x = bioclim_PC3, y = bioclim_PC4, colour = Mountain_range)
  ) + xlab('PC3') + ylab('PC4') +
    geom_segment(
      pca_variables,
      mapping = aes(
        x = 0,
        xend = PC3 * 5,
        y = 0,
        yend = PC4 * 5
      )
    ) +
    geom_text(
      pca_variables,
      mapping = aes(
        x = PC3 * 5,
        y = PC4 * 1.2 * 5,
        label = Variable
      )
    ) + theme_minimal() + theme(legend.position = 'none') +
    scale_colour_manual(
      values = c(
        '#CC2F50',
        '#1A5A61',
        '#E3C78D',
        '#4BC992',
        '#806A2A',
        '#177FCF',
        '#982DA1',
        '#A7BFE8',
        'dimgrey',
        'tomato',
        'grey'
      ),
      name = 'Mountain range'
    )
  
  p3 = ggplot() + geom_point(
    all_data_with_bioclim_pca,
    mapping = aes(x = bioclim_PC5, y = bioclim_PC6, colour = Mountain_range)
  ) + xlab('PC5') + ylab('PC6') +
    geom_segment(
      pca_variables,
      mapping = aes(
        x = 0,
        xend = PC5 * 3,
        y = 0,
        yend = PC6 * 3
      )
    ) +
    geom_text(
      pca_variables,
      mapping = aes(
        x = PC5 * 3,
        y = PC6 * 1.2 * 3,
        label = Variable
      )
    ) + theme_minimal() + theme(legend.position = 'none') +
    scale_colour_manual(
      values = c(
        '#CC2F50',
        '#1A5A61',
        '#E3C78D',
        '#4BC992',
        '#806A2A',
        '#177FCF',
        '#982DA1',
        '#A7BFE8',
        'dimgrey',
        'tomato',
        'grey'
      ),
      name = 'Mountain range'
    )
  
  p = ggarrange(p1 + theme(legend.position = 'none'),
                p2,
                p3,
                leg,
                nrow = 2,
                ncol = 2)
  ggsave(
    plot = p,
    'plots/Fig_S4_biolimatic_pca.pdf',
    width = 12,
    height = 12
  )
}

log_transform <- function(x) {
  min_non_zero = min(x[x > 0], na.rm = T)
  return(log(x + (min_non_zero / 2)))
}

TransformData <- function(data_raw) {
  data = data_raw %>% distinct()
  data = data[!(is.na(data$min_quartz) & is.na(data$min_calcite) & is.na(data$min_feldspar) & is.na(data$min_clays)), ]
  
  colmins = data %>% summarise_if(is.numeric, function(x) {
    return(min(x[which(x > 0)], na.rm = T))
  }) / 2
  saveRDS(colmins, "data/processed/colmins.rds")
  
  ####### Climatic and glaciological
  data$gl_distance = log_transform(data$gl_distance)
  data$gl_area = log_transform(data$gl_area)
  data$gl_coverage = sin(sqrt(data$gl_coverage))
  data$clim_pr = log_transform(data$clim_pr)
  
  ####### Stream params and nutrients
  data$pc_water_temp = log_transform(data$pc_water_temp)
  data$pc_conductivity = log_transform(data$pc_conductivity)
  data$pc_turbidity = log_transform(data$pc_turbidity)
  data$nut_srp = log_transform(data$nut_srp)
  data$nut_din = log_transform(data$nut_din)
  
  ####### Biomass
  data$bacterial_abundance = log_transform(data$bacterial_abundance)
  data$chla = log_transform(data$chla)
  
  ####### Mineralogy
  data$min_calcite = log_transform(data$min_calcite)
  data$min_clays = asin(sqrt(data$min_clays))
  data$min_feldspar = log_transform(data$min_feldspar)
  data$min_quartz = log_transform(data$min_quartz)
  
  return(data)
}

ProcessMagData <- function() {
  mag_data = read.csv('data/raw/bacteria/MAGs_cov_norm.txt', sep = '\t')
  tree = read.tree('data/raw/bacteria/treeBacteria.tree')
  tree$tip.label = map_chr(tree$tip.label, function(x)
    strsplit(x, '.fa')[[1]][1])
  
  rownames(mag_data) = mag_data$MAGs
  mag_data$MAGs = NULL
  mag_data = mag_data[, !startsWith(colnames(mag_data), 'GLR')]
  mag_data = mag_data[, !startsWith(colnames(mag_data), 'X3')]
  colnames(mag_data)[colnames(mag_data) == 'GL140_1'] = 'GL140_UP'
  colnames(mag_data)[colnames(mag_data) == 'GL140_2'] = 'GL140_UP'
  colnames(mag_data)[colnames(mag_data) == 'GL140_3'] = 'GL140_UP'
  colnames(mag_data)[colnames(mag_data) == 'GL140_4'] = 'GL140_DN'
  colnames(mag_data)[colnames(mag_data) == 'GL140_5'] = 'GL140_DN'
  colnames(mag_data)[colnames(mag_data) == 'GL140_6'] = 'GL140_DN'
  colnames(mag_data) = gsub('DownB', 'DN', colnames(mag_data))
  colnames(mag_data) = gsub('Down', 'DN', colnames(mag_data))
  colnames(mag_data) = gsub('UpB', 'UP', colnames(mag_data))
  colnames(mag_data) = gsub('Up', 'UP', colnames(mag_data))
  mag_data = as.data.frame(do.call(cbind, by(
    t(mag_data), INDICES = names(mag_data), FUN = colMeans
  )))
  
  colsums_all = colSums(mag_data)
  mag_data = mag_data[rowSums(mag_data) > 0, ]
  mag_data = mag_data[rowMeans(mag_data > 10) > 0.2, ]
  colsums_filtered = colSums(mag_data)
  print(mean(colsums_filtered / colsums_all)) # 0.9851491 of counts are kept, 2333 / 2868 ASVs
  
  write.csv(mag_data, file = 'data/processed/MAG_cov_filtered.tsv', quote = F)
  return(list(tree=tree, mag_data=mag_data))}

AddRichness <- function(data, tree, mag_data) {
  data$Shannon = map_dbl(data$Sample, function(x)
    ifelse(x %in% colnames(mag_data),
           diversity(mag_data[, colnames(mag_data) == x]),-999))
  data$Pielou = map_dbl(data$Sample, function(x)
    ifelse(
      x %in% colnames(mag_data),
      diversity(mag_data[, colnames(mag_data) == x]) / log(sum(mag_data[, colnames(mag_data) == x] > 0)),-999
    ))
  
  tree = keep.tip(tree, rownames(mag_data)[rownames(mag_data) %in% tree$tip.label])
  tree = midpoint.root(tree)
  mag_data = mag_data[rownames(mag_data) %in% tree$tip.label, ]
  mag_data = mag_data[order(match(rownames(mag_data), tree$tip.label)), ]
  coph = cophenetic(tree)
  
  mntd_df =  mntd(t(mag_data), coph, abundance.weighted = T)
  data$mntd = map_dbl(data$Sample, function(x)
    ifelse(x %in% colnames(mag_data),
           mntd_df[colnames(mag_data) == x],-999))
  mpd_df = mpd(t(mag_data), coph, abundance.weighted = T)
  data$mpd = map_dbl(data$Sample, function(x)
    ifelse(x %in% colnames(mag_data),
           mpd_df[colnames(mag_data) == x],-999))
  pd_df = pd(t(mag_data), tree, include.root = F)
  data$pd = map_dbl(data$Sample, function(x)
    ifelse(x %in% colnames(mag_data),
           pd_df[colnames(mag_data) == x, 'PD'],-999))
  data$Shannon[data$Shannon == -999] = NA
  data$Pielou[data$Pielou == -999] = NA
  data$mntd[data$mntd == -999] = NA
  data$mpd[data$mpd == -999] = NA
  data$pd[data$pd == -999] = NA
  
  data$pd[data$pd < 292] = NA
  data$pd = exp(data$pd)
  data$mpd = exp(data$mpd)
  
  data$Shannon[data$Date == 'Future'] = NA
  data$Pielou[data$Date == 'Future'] = NA
  data$mntd[data$Date == 'Future'] = NA
  data$mpd[data$Date == 'Future'] = NA
  data$pd[data$Date == 'Future'] = NA
  return(data)
}

MainD <- function() {
  dir.create('plots')
  dir.create('stats')
  
  # Load data, perform the bioclimatic PCA
  all_data = read.csv('data/processed/all_current_data_3_ssps.csv')
  colnames(all_data)[colnames(all_data) == 'SSP'] = 'Scenario'
  colnames(all_data)[colnames(all_data) == 'gl_dist'] = 'gl_distance'
  
  all_data_with_pca = BioclimPCA(all_data)
  all_data_with_bioclim_pca = all_data_with_pca$all_data
  pca_variables = all_data_with_pca$pca_vars
  PlotBioclimPCA(all_data_with_bioclim_pca, pca_variables)
  
  # Add MAG and richness information
  processed_mag_data = ProcessMagData()
  all_data_with_bioclim_pca_and_richness = AddRichness(all_data_with_bioclim_pca, processed_mag_data$tree, processed_mag_data$mag_data)
  all_data_processed = TransformData(all_data_with_bioclim_pca_and_richness)
  
  write.csv(all_data_processed, 'data/processed/all_current_clean_3_ssps.csv', quote = F)}
