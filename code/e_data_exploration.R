library(ggplot2)
library(ggpubr)
library(maps)
library(dplyr)
library(mapproj)

PlotPresentFuture <- function(data, var, var_name, log=FALSE){
  plot_data = data.frame(sample = data %>% filter(Site == 'UP', Scenario == 370, Date == 'Present') %>% pull('Sample'), 
                         Mountain_range = data %>% filter(Site == 'UP', Scenario == 370, Date == 'Present') %>% pull('Mountain_range')) 
  
  plot_data$present_126 = data %>% filter(Site == 'UP', Scenario == 126, Date == 'Present') %>% pull(var)
  plot_data$future_126 = data %>% filter(Site == 'UP', Scenario == 126, Date == 'Future') %>% pull(var)
  plot_data$present_370 = data %>% filter(Site == 'UP', Scenario == 370, Date == 'Present') %>% pull(var)
  plot_data$future_370 = data %>% filter(Site == 'UP', Scenario == 370, Date == 'Future') %>% pull(var)
  plot_data$present_585 = data %>% filter(Site == 'UP', Scenario == 585, Date == 'Present') %>% pull(var)
  plot_data$future_585 = data %>% filter(Site == 'UP', Scenario == 585, Date == 'Future') %>% pull(var)
  
  plot_data$change_126 = (plot_data$future_126 - plot_data$present_126)
  plot_data$change_370 = (plot_data$future_370 - plot_data$present_370)
  plot_data$change_585 = (plot_data$future_585 - plot_data$present_585)

  stat_126 = wilcox.test(plot_data$change_126, mu = 0)
  stat_370 = wilcox.test(plot_data$change_370, mu = 0)
  stat_585 = wilcox.test(plot_data$change_585, mu = 0)

  stats_df = data.frame(Variable = c(var, var, var),
                        RCP=c('2.6', '4.5', '8.5'),
                        median_change = c(median(plot_data$change_126),
                                          median(plot_data$change_370),
                                          median(plot_data$change_585)),
                        relative_change = c(median(plot_data$change_126 / plot_data$present_126),
                                            median(plot_data$change_370 / plot_data$present_370),
                                            median(plot_data$change_585 / plot_data$present_585)),
                        p_val = c(stat_126$p.value, stat_370$p.value, stat_585$p.value),
                        q25_change = c(quantile(plot_data$change_126, probs = 0.25),
                                       quantile(plot_data$change_370, probs = 0.25),
                                       quantile(plot_data$change_585, probs = 0.25)),
                        q75_change = c(quantile(plot_data$change_126, probs = 0.75),
                                       quantile(plot_data$change_370, probs = 0.75),
                                       quantile(plot_data$change_585, probs = 0.75)))
  
  p <- ggplot(plot_data, aes(x=present_370, y=future_370, colour=Mountain_range)) + geom_point(size=2) +
    geom_errorbar(aes(ymin=future_126, ymax=future_585), alpha=0.2) +
    geom_errorbarh(aes(xmin=present_126, xmax=present_585), alpha=0.2) +
    scale_x_log10() + scale_y_log10() + ggtitle(var_name) +
    scale_colour_manual(values = c('#CC2F50','#1A5A61','#E3C78D','#4BC992','#806A2A','#177FCF','#982DA1','#A7BFE8', 'dimgrey', 'tomato', 'grey'), 
                        name='Mountain range') +
    geom_abline(intercept = 0, linetype = 'dashed') + 
    xlab('Sampling year (2019-2022)') + ylab('Average for the 2070-2100 period') + theme_linedraw() + theme(panel.grid.major = element_blank(), 
                                                                                                            panel.grid.minor = element_blank(),
                                                                                                            panel.background = element_blank())
  if (log == TRUE){p = p + scale_x_log10() + scale_y_log10()}
  return(list(plot=p, stats=stats_df))}

PlotMap <- function(data){
  WorldData <- ggplot2::map_data('world') %>% fortify
  count_meta = data %>% filter(Site == 'UP', Scenario == 370, Date == 'Present') %>% group_by(Mountain_range) %>% summarise(latitude = mean(latitude), longitude = mean(longitude), 'Sample n.' = n())

  p1 = ggplot() +
    geom_map(data = WorldData, map = WorldData, aes(long, lat, group = group, map_id = region), fill = "dimgrey", colour = "dimgrey", size=0.2) +
    coord_map("rectangular", lat0=0, xlim=c(-180,180), ylim=c(-90, 90)) + xlab('') + ylab('') + 
    geom_point(data=count_meta, aes(x=longitude, y=latitude, colour = Mountain_range, size=`Sample n.`), alpha=1, size=9) + 
    geom_text(data=count_meta, aes(label=`Sample n.`, x=longitude, y=latitude),color='white', size=5)+
    scale_colour_manual(values = c('#CC2F50','#1A5A61','#E3C78D','#4BC992','#806A2A','#177FCF','#982DA1','#A7BFE8', 'dimgrey', 'tomato', 'grey'),
                        name='Mountain range') + 
    theme_bw() + theme(legend.title = element_text(size=8.5), 
                       legend.text=element_text(size=7.5), 
                       axis.title=element_text(size=7.5), 
                       legend.position = c(0.105, 0.4),
                       axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(),
                       axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank(),
                       legend.margin=margin(t = 0, unit='cm')) +
    guides(colour = guide_legend(override.aes = list(alpha=1, size=3)), shape = guide_legend(override.aes = list(size=3)), size = FALSE) 
  
  gl_dist = PlotPresentFuture(data, 'gl_distance', 'Distance to glacier terminus (m)', log=FALSE)
  sp2 = gl_dist$plot
  
  gl_area = PlotPresentFuture(data, 'gl_area', 'Glacier area (km2)', log=FALSE)
  sp3 = gl_area$plot
  
  p2 = ggarrange(sp2, sp3, labels = c('B','C'), nrow = 1, ncol = 2, legend = 'none', hjust=-0.2)
  p = ggarrange(p1, p2, labels = c('A',''), nrow = 2, ncol = 1, hjust=-0.2)
  ggsave(p, filename = 'plots/Fig_1_map_glacial_changes.pdf', width = 8, height = 8)
  
  return(list(gl_dist=gl_dist$stats, gl_area=gl_area$stats))}

MainE <- function(){
    data = read.csv('data/processed/all_current_data_3_ssps.csv')
    colnames(data)[colnames(data) == 'SSP'] = 'Scenario'
    colnames(data)[colnames(data) == 'gl_dist'] = 'gl_distance'
    gl_stats = PlotMap(data)
    
    clim_scd = PlotPresentFuture(data, 'clim_scd', ""*Yearly~snow~cover~days~(N)*"", log=FALSE)
    p1 = clim_scd$plot
    
    clim_pr = PlotPresentFuture(data, 'clim_pr', ""*Monthly~precipitations~(kg/m2/month)*"", log=FALSE)
    p2 = clim_pr$plot
    
    clim_tas = PlotPresentFuture(data, 'clim_tas', ""*Average~temperature~(degree*C)*"", log=FALSE)
    p3 = clim_tas$plot
    
    clim_tasmin = PlotPresentFuture(data, 'clim_tasmin', ""*Average~daily~minimal~temperature(degree*C)*"", log=FALSE)
    p4 = clim_tasmin$plot
    
    clim_tasmax = PlotPresentFuture(data, 'clim_tasmax', ""*Average~daily~maximal~temperature~(degree*C)*"", log=FALSE)
    p5 = clim_tasmax$plot
    
    leg = get_legend(p1)
    p = ggarrange(p1 + theme(legend.position = 'none'), p2 + theme(legend.position = 'none'), p3 + theme(legend.position = 'none'), 
                  p4 + theme(legend.position = 'none'), p5 + theme(legend.position = 'none'), leg,
                  labels = c('A','B','C','D','E'), nrow = 2, ncol = 3, vjust = 2)
    ggsave(p, filename = 'plots/Fig_S1_climate_variables.pdf', width = 12, height = 8)
    
    stats_df = rbind(gl_stats$gl_dist$stats, gl_stats$gl_area$stats, clim_scd$stats, clim_pr$stats, clim_tas$stats, clim_tasmin$stats, clim_tasmax$stats)
    write.csv(stats_df, file = 'stats/gl_clim_median_changes.csv', quote = F, row.names = F)}

