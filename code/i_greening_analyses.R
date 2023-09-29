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

########################################################################################################################
# Greening functions
DataGreening <- function(data){
  up_data = data %>% filter(Site == 'UP') %>% arrange(Glacier)
  up_data_present = up_data %>% filter(Date == 'Present')
  up_data_future = up_data %>% filter(Date == 'Future')
  chla_turb = data.frame(Turbidity = c(up_data_present$pc_turbidity_predicted,
                                       up_data_future$pc_turbidity_predicted),
                          Bacterial_abundance = c(up_data_present$bacterial_abundance_predicted,
                                                  up_data_future$bacterial_abundance_predicted),
                          Chlorophyll = c(up_data_present$chla_predicted,
                                          up_data_future$chla_predicted),
                          Shannon = c(up_data_present$Shannon_predicted,
                                      up_data_future$Shannon_predicted),
                          mntd = c(up_data_present$mntd_predicted,
                                    up_data_future$mntd_predicted),
                          Scenario = c(up_data_present$Scenario,
                                       up_data_future$Scenario),
                          Date = c(rep('Sampling year (2019-2022)', length(up_data_present$Date)),
                                    rep('Average for the 2070-2100 period', length(up_data_future$Date))),
                          Sample = c(up_data_present$Sample,
                                      up_data_future$Sample))
  chla_turb$Bacterial_abundance = log10(exp(chla_turb$Bacterial_abundance))
  return(chla_turb)}

PlotGreening <- function(greening_data){
  greening_data = greening_data %>% filter(Scenario == 370)
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
         y = ""*Chlorophyll-italic(a)~(ln~mu*g~g^-1)*"") + theme(legend.title = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank())

  p3 = ggplot(greening_data, aes(y=Shannon, x=Chlorophyll, colour=Date)) + 
    geom_line(aes(group=Sample), size=0.5, colour='black', alpha=0.1) + geom_point(alpha=0.7, size=3) + 
    theme_bw() + scale_colour_manual(values=c('#15356B', '#E04D55')) +
    labs(y = ""*Shannon~index*"",
         x = ""*Chlorophyll-italic(a)~(ln~mu*g~g^-1)*"") + theme(legend.title = element_blank(),
                                                                  panel.grid.major = element_blank(),
                                                                  panel.grid.minor = element_blank())

  p4 = ggplot(greening_data, aes(y=mntd, x=Chlorophyll, colour=Date)) + 
    geom_line(aes(group=Sample), size=0.5, colour='black', alpha=0.1) + geom_point(alpha=0.7, size=3) + 
    theme_bw() + scale_colour_manual(values=c('#15356B', '#E04D55')) +
    labs(y = bquote(""*alpha~-~MNTD*""),
         x = ""*Chlorophyll-italic(a)~(ln~mu*g~g^-1)*"") + theme(legend.title = element_blank(),
                                                                  panel.grid.major = element_blank(),
                                                                  panel.grid.minor = element_blank())

  p = ggarrange(p1, p2, p3, p4, labels = c('A','B', 'C', 'D'), ncol = 2, nrow = 2, common.legend = T, align = 'h') 
  ggsave(p, filename = 'plots/Fig_3_greening.pdf', width = 8, height = 8)}

GreeningCorrelations <- function(greening_data){
  cor_results = data.frame()
  for (scenario in c(126, 370 ,585)){
    scenario_data_present = greening_data %>% filter(Scenario == scenario, Date == 'Sampling year (2019-2022)')
    scenario_data_future = greening_data %>% filter(Scenario == scenario, Date == 'Average for the 2070-2100 period')
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
                                  cor = c(turb_cor$estimate, baca_cor$estimate, shan_cor$estimate, mntd_cor$estimate),
                                  p = c(turb_cor$p.value, baca_cor$p.value, shan_cor$p.value, mntd_cor$p.value))
    cor_results = rbind(cor_results, scenario_results)}
  return(cor_results)}

MainI <- function(){
  data = read.csv('data/processed/all_projections_3_ssps.csv')
  greening_data = DataGreening(data)
  PlotGreening(greening_data)
  cors = GreeningCorrelations(greening_data)
  write.csv(cors, 'stats/greening_cors.csv', quote=F, row.names=F)}
