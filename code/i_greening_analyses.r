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
source('code/f_stream_parameters_models.R')

########################################################################################################################
# Greening functions
DataGreening <- function(data){
  up_data = pred_data %>% filter(Site == 'UP')
  up_data_present = up_data %>% filter(Date == 'Present')
  up_data_future = up_data %>% filter(Date == 'Future')
  chla_turb = data.frame(Turbidity = c(up_data_present$pc_turbidity_predicted,
                                        up_data_future$pc_turbidity_predicted),
                          Bacterial_abundance = c(up_data_present$bacterial_abundance_predicted,
                                                  up_data_future$bacterial_abundance_predicted),
                          Chlorophyll = c(up_data_present$chla_predicted,
                                          up_data_future$chla_predicted),
                          Shannon = c(up_data_present$shannon_predicted,
                                      up_data_future$shannon_predicted),
                          mntd = c(up_data_present$mntd_predicted,
                                    up_data_future$mntd_predicted),
                          SSP = c(up_data_present$SSP,
                                  up_data_future$SSP),
                          Date = c(rep('Sampling year (2019-2022)', length(up_data_present$Date)),
                                    rep('Average for the 2070-2100 period', length(up_data_future$Date))),
                          Sample = c(up_data_present$Sample,
                                      up_data_future$Sample))
    chla_turb$Bacterial_abundance = log10(exp(chla_turb$Bacterial_abundance))
  return(chla_turb)}

PlotGreening <- function(greening_data){
  greening_data = greening_data %>% filter(SSP == 370)
  p1 = ggplot(greening_data, aes(x=Turbidity, y=Chlorophyll, colour=Date)) + 
    geom_line(aes(group=Sample), size=0.5, colour='black', alpha=0.1) + geom_point(alpha=0.7, size=3) + 
    theme_bw() + scale_colour_manual(values=c('#15356B', '#E04D55')) +
    labs(x = ""*Turbidity~(ln~NTU)*"",
        y = ""*Chlorophyll-italic(a)~(ln~mu*g~g^-1)*"") + theme(legend.title = element_blank(),
                                                                panel.grid.major = element_blank(),
                                                                panel.grid.minor = element_blank())

  p2 = ggplot(greening_data, aes(x=Bacterial_abundance, y=Chlorophyll, colour=Date)) + 
    geom_line(aes(group=Sample), size=0.5, colour='black', alpha=0.1) + geom_point(alpha=0.7, size=3) + 
    theme_bw() + scale_colour_manual(values=c('#15356B', '#E04D55')) +
    labs(x = ""*Bacterial~abundance~(log[10]~cells~g^-1)*"",
        y = "") + theme(legend.title = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank())

  p3 = ggplot(greening_data, aes(x=Shannon, y=Chlorophyll, colour=Date)) + 
    geom_line(aes(group=Sample), size=0.5, colour='black', alpha=0.1) + geom_point(alpha=0.7, size=3) + 
    theme_bw() + scale_colour_manual(values=c('#15356B', '#E04D55')) +
    labs(x = ""*Shannon~index*"",
        y = "") + theme(legend.title = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank())

  p4 = ggplot(greening_data, aes(x=mntd, y=Chlorophyll, colour=Date)) + 
    geom_line(aes(group=Sample), size=0.5, colour='black', alpha=0.1) + geom_point(alpha=0.7, size=3) + 
    theme_bw() + scale_colour_manual(values=c('#15356B', '#E04D55')) +
    labs(x = ""*Shannon~index*"",
        y = "") + theme(legend.title = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank())

  p = ggarrange(p1, p2, p3, p4, labels = c('A','B', 'C', 'D'), ncol = 2, nrow = 2, common.legend = T, align = 'h') 
  ggsave(p, filename = 'Plots/Fig_3_chla_prediction.pdf', width = 8, height = 8)}

GreeningCorrelations <- function(greening_data){
  cor_results = data.frame()
  for (scenario in c(126, 370 ,585)){
    scenario_data_present = greening_data %>% filter(SSP == scenario, Date == 'Sampling year (2019-2022)')
    scenario_data_future = greening_data %>% filter(SSP == scenario, Date == 'Average for the 2070-2100 period')
    turb_cor = cor.test(scenario_data_future$Chlorophyll - scenario_data_present$Chlorophyll,
                        scenario_data_future$Turbidity - scenario_data_present$Turbidity, method='spearman')
    baca_cor = cor.test(scenario_data_future$Chlorophyll - scenario_data_present$Chlorophyll,
                        scenario_data_future$Bacterial_abundance - scenario_data_present$Bacterial_abundance, method='spearman')
    shan_cor = cor.test(scenario_data_future$Chlorophyll - scenario_data_present$Chlorophyll,
                        scenario_data_future$Shannon - scenario_data_present$Shannon, method='spearman')
    mntd_cor = cor.test(scenario_data_future$Chlorophyll - scenario_data_present$Chlorophyll,
                        scenario_data_future$mntd - scenario_data_present$mntd, method='spearman')
    scenario_results = data.frame(Variable = c('Turbidity', 'Bacterial abundance', 'Shannon', 'MNTD'),
                                  Scenario = rep(scenario, 4), 
                                  cor = c(turb_cor$statistic, baca_cor$statistic, shan_cor$statistic, mntd_cor$statistic),
                                  p = c(turb_cor$p.value, baca_cor$p.value, shan_cor$p.value, mntd_cor$p.value))
    cor_results = rbind(cor_results, scenario_results)}
  return(cor_results)}

GreeningCorrelations <- function(greening_data){
  cor_results = data.frame()
  for (scenario in c(126, 370 ,585)){
    scenario_data_present = greening_data %>% filter(SSP == scenario, Date == 'Sampling year (2019-2022)')
    scenario_data_future = greening_data %>% filter(SSP == scenario, Date == 'Average for the 2070-2100 period')
    turb_cor = cor.test(scenario_data_future$Chlorophyll - scenario_data_present$Chlorophyll,
                        scenario_data_future$Turbidity - scenario_data_present$Turbidity, method='spearman')
    baca_cor = cor.test(scenario_data_future$Chlorophyll - scenario_data_present$Chlorophyll,
                        scenario_data_future$Bacterial_abundance - scenario_data_present$Bacterial_abundance, method='spearman')
    shan_cor = cor.test(scenario_data_future$Chlorophyll - scenario_data_present$Chlorophyll,
                        scenario_data_future$Shannon - scenario_data_present$Shannon, method='spearman')
    mntd_cor = cor.test(scenario_data_future$Chlorophyll - scenario_data_present$Chlorophyll,
                        scenario_data_future$mntd - scenario_data_present$mntd, method='spearman')
    scenario_results = data.frame(Variable = c('Turbidity', 'Bacterial abundance', 'Shannon', 'MNTD'),
                                  Scenario = rep(scenario, 4), 
                                  cor = c(turb_cor$statistic, baca_cor$statistic, shan_cor$statistic, mntd_cor$statistic),
                                  p = c(turb_cor$p.value, baca_cor$p.value, shan_cor$p.value, mntd_cor$p.value))
    cor_results = rbind(cor_results, scenario_results)}
  return(cor_results)}
