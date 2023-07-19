library(ggplot2)
library(ggpubr)
library(dplyr)
library(mgcv)
library(reshape2)
library(matrixStats)
library(purrr)
library(glmnet)
library(doMC)
library(foreach)
library(ggridges)
source('code/f_stream_parameters_models.R')

variables =  c('pc_water_temp', 'pc_turbidity', 'pc_conductivity', 'pc_ph', 'nut_din', 'nut_srp', 
              'chla', 'bacterial_abundance', 'Shannon', 'Pielou', 'mntd', 'mpd')

########################################################################################################################
# Stats functions
ComputeChanges <- function(up_data, var, colmins){
  median_changes = data.frame()
  for(ssp in c(126, 370, 585)){
    future_vals = up_data %>% filter(Scenario == ssp, Date == 'Future') %>% pull(paste0(var, '_predicted', collapse = ''))
    present_vals = up_data %>% filter(Scenario == ssp, Date == 'Present') %>% pull(paste0(var, '_predicted', collapse = ''))
    if (!(var %in% c('pc_ph', 'bacterial_abundance', 'Pielou', 'Shannon', 'mntd', 'mpd'))){
      future_vals = vapply(future_vals, function(x) exp(x) - as.double(colmins[var]), FUN.VALUE = numeric(1))
      present_vals = vapply(present_vals, function(x) exp(x) - as.double(colmins[var]), FUN.VALUE = numeric(1))}
    if (var  == 'bacterial_abundance'){
      future_vals = vapply(future_vals, function(x) exp(x) - 1, FUN.VALUE = numeric(1))
      present_vals = vapply(present_vals, function(x) exp(x) - 1, FUN.VALUE = numeric(1))}
    
    median_change = median(future_vals - present_vals)
    q25_change = quantile(future_vals - present_vals, probs = 0.25)
    q75_change = quantile(future_vals - present_vals, probs = 0.75)
    
    median_relative_change =  median((future_vals - present_vals)/present_vals)
    q25_relative_change =  quantile((future_vals - present_vals)/present_vals, probs = 0.25)
    q75_relative_change =  quantile((future_vals - present_vals)/present_vals, probs = 0.75)
    
    diff = future_vals - present_vals
    p_val = wilcox.test(diff)$p.value
    print(p_val)
    median_changes = rbind(median_changes, data.frame(Variable=var, Scenario=ssp, median_change=round(median_change,3), 
                                                      q25_change=round(q25_change,3), q75_change=round(q75_change,3),
                                                      median_relative_change=round(median_relative_change,3),
                                                      q25_relative_change=round(q25_relative_change,3), 
                                                      q75_relative_change=round(q75_relative_change,3),
                                                      p_val=p_val))}
  return(median_changes)}

GetChangesResults <- function(data, variables){
  colmins <- readRDS("data/processed/colmins.rds")
  up_data = data %>% filter(Site == 'UP')
  median_changes = data.frame()
  for (var in variables){
    var_results = ComputeChanges(up_data, var, colmins)
    median_changes = rbind(median_changes, var_results)}
  write.csv(median_changes, 'stats/stream_median_changes.csv', quote = F, row.names = F)}

########################################################################################################################
# Plot distributions functions
plot_distribution <- function(data, resp_var, x_lab, strip = T){
    p = ggplot(data[data$variable == resp_var,], aes(x=value, y=Date, fill=Date)) + 
      geom_density_ridges(quantile_lines = TRUE, quantiles = 2, alpha=0.75) + theme_linedraw() + scale_fill_manual(values =c('#2F5494', '#FA726B')) + facet_grid(~Scenario) +
      ylab('') + xlab(x_lab) + theme(axis.text.y = element_blank(),
                                    legend.title = element_blank(),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),
                                    panel.spacing = unit(1.5, "lines"),
                                    plot.margin = unit(c(0.25,0.25,0.5,0.25), units = "lines" ))
    if (strip == F){p = p + theme(strip.text = element_blank(), 
                                  strip.background = element_blank())}
    return(p)}

PlotPresentFutureComparison <- function(data){
  # Future/Present distributions
  plot_env_data = data %>% group_by(Date, Scenario) %>% 
    select(paste0(variables, '_predicted')) %>% 
    melt(measure.vars = paste0(variables, '_predicted'))
  plot_env_data$Scenario[plot_env_data$Scenario == 126] = 'RCP 2.6'
  plot_env_data$Scenario[plot_env_data$Scenario == 370] = 'RCP 4.5'
  plot_env_data$Scenario[plot_env_data$Scenario == 585] = 'RCP 8.5'

  plot_env_data$Date[plot_env_data$Date == 'Future'] = 'Average for the 2070-2100 period'
  plot_env_data$Date[plot_env_data$Date == 'Present'] = 'Sampling year (2019-2022)'

  p1 = plot_distribution(plot_env_data, 'pc_water_temp_predicted', ""*Water~temperature~(ln~degree*C)*"", T)
  p2 = plot_distribution(plot_env_data, 'pc_turbidity_predicted', ""*Turbidity~(ln~NTU)*"", T)
  p3 = plot_distribution(plot_env_data, 'pc_conductivity_predicted', ""*Conductivity~(ln~mu*S~cm^-1)*"", T)
  p4 = plot_distribution(plot_env_data, 'pc_ph_predicted', "pH", T)
  p5 = plot_distribution(plot_env_data, 'nut_din_predicted', ""*Dissolved~inorganic~nitrogen~(ln~mu*g~L^-1)*"", T)
  p6 = plot_distribution(plot_env_data, 'nut_srp_predicted', ""*Soluble~reactive~phosphate~(ln~mu*g~L^-1)*"", T)
  p7 = plot_distribution(plot_env_data, 'chla_predicted', ""*Chlorophyll-italic(a)~(ln~mu*g~g^-1)*"", T)
  p8 = plot_distribution(plot_env_data, 'bacterial_abundance_predicted', ""*Bacterial~abundance~(ln~cells~g^-1)*"", T)
  p9 = plot_distribution(plot_env_data, 'Shannon_predicted', ""*Shannon~index*"", T)
  p10 = plot_distribution(plot_env_data, 'Pielou_predicted', ""*Pielou~evenness*"", T)
  p11 = plot_distribution(plot_env_data, 'mntd_predicted', ""*alpha~-~MNTD*"", T)
  p12 = plot_distribution(plot_env_data, 'mpd_predicted', ""*alpha~-~MPD*"", T)
  p = ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, common.legend = T, ncol = 1) + theme(plot.margin = margin(0.1,0.5,0.2,0, "cm")) 
  ggsave(p, filename = 'plots/Fig_S2_env_scenarios.pdf', width = 6, height = 15)

  # Future/Present plots
  plot_env_present = data %>% filter(Date == 'Present', Site == 'UP') %>% group_by(Scenario, Mountain_range) %>%
    select(pc_water_temp_predicted, pc_turbidity_predicted, pc_conductivity_predicted, pc_ph_predicted, 
          nut_din_predicted, nut_srp_predicted, Mountain_range, Scenario) %>% 
    melt(measure.vars = c('pc_water_temp_predicted', 'pc_turbidity_predicted', 'pc_conductivity_predicted', 
                          'pc_ph_predicted', 'nut_din_predicted', 'nut_srp_predicted'), value.name = 'Present', id.vars = c('Scenario', 'Mountain_range'))
  plot_env_future = data %>% filter(Date == 'Future', Site == 'UP') %>% group_by(Scenario, Mountain_range)  %>%
    select(pc_water_temp_predicted, pc_turbidity_predicted, pc_conductivity_predicted, pc_ph_predicted, 
          nut_din_predicted, nut_srp_predicted, Mountain_range, Scenario) %>% 
    melt(measure.vars = c('pc_water_temp_predicted', 'pc_turbidity_predicted', 'pc_conductivity_predicted', 
                          'pc_ph_predicted', 'nut_din_predicted', 'nut_srp_predicted'), value.name = 'Future', id.vars = c('Scenario', 'Mountain_range'))
  plot_env_preds = plot_env_present
  colnames(plot_env_preds) = c('Scenario', 'Mountain_range', 'variable', 'Present')
  plot_env_preds$Future = plot_env_future$Future

  plot_env_preds$Scenario[plot_env_preds$Scenario == 126] = 'RCP 2.6'
  plot_env_preds$Scenario[plot_env_preds$Scenario == 370] = 'RCP 4.5'
  plot_env_preds$Scenario[plot_env_preds$Scenario == 585] = 'RCP 8.5'

  labels_facets <- list(
    'pc_water_temp'=""*A.~~~Water~temperature~(ln~degree*C)*"",
    'pc_turbidity'= ""*B.~~~Turbidity~(ln~NTU)*"",
    'pc_conductivity'=""*C.~~~Conductivity~(ln~mu*S~cm^-1)*"",
    'pc_ph'="D.   pH",
    'nut_din'=""*E.~~~Dissolved~inorganic~nitrogen~(ln~mu*g~L^-1)*"",
    'nut_srp'=""*F.~~~Soluble~reactive~phosphate~(ln~mu*g~L^-1)*"")
  facets_labeller <- function(variable,value){
    return(labels_facets[value])}

  p = ggplot(plot_env_preds[plot_env_preds$Scenario == 'RCP 4.5',], aes(x=Present, y=Future)) + 
    facet_wrap(~variable,  ncol=2, nrow=3, scales = "free", labeller = facets_labeller) + geom_abline(intercept = 0, slope = 1, colour='black', linetype='dashed') +
    geom_point(aes(colour=Mountain_range), size = 2, alpha = 0.5) + 
    geom_quantile(quantiles = 0.5, colour = 'black', size=1) + 
    theme_linedraw() + theme(panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            strip.text.x = element_text(hjust = 0, vjust = 0.25),
                            legend.position = 'none') + 
    scale_colour_manual(values = c('#CC2F50','#1A5A61','#E3C78D','#4BC992','#806A2A','#177FCF','#982DA1','#A7BFE8', 'dimgrey', 'tomato', 'grey')) +
    xlab('Sampling year (2019-2022)') +
    ylab('Average for the 2070-2100 period')
  ggsave(plot = p, 'plots/Fig_2_env_changes.pdf', width = 5.75, height = 8.5)}

MainG <- function(){
  data = read.csv('data/processed/all_projections_3_ssps.csv')
  GetChangesResults(data, variables)
  PlotPresentFutureComparison(data)}
