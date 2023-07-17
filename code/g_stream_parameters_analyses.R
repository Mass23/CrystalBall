# Median changes
up_data = rbind(up_data_present, up_data_future)
median_changes = data.frame()
for (var in c('pc_water_temp', 'pc_turbidity', 'pc_conductivity', 'pc_ph', 'nut_din', 'nut_srp', 
              'chla', 'bacterial_abundance', 'shannon', 'pielou', 'mntd', 'mpd', 'pd')){
  print('------------------------')
  for(ssp in c(126, 370, 585)){
    future_vals = up_data[(up_data$SSP == ssp) & (up_data$Date == 'Future'), paste0(var, '_predicted', collapse = '')]
    present_vals = up_data[(up_data$SSP == ssp) & (up_data$Date == 'Present'), paste0(var, '_predicted', collapse = '')]
    if (!(var %in% c('pc_ph', 'bacterial_abundance', 'pielou', 'shannon', 'mntd', 'mpd', 'pd'))){
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
    median_changes = rbind(median_changes, data.frame(Variable=var, SSP=ssp, median_change=round(median_change,3), 
                                                      q25_change=round(q25_change,3), q75_change=round(q75_change,3),
                                                      median_relative_change=round(median_relative_change,3),
                                                      q25_relative_change=round(q25_relative_change,3), 
                                                      q75_relative_change=round(q75_relative_change,3),
                                                      p_val=p_val))}}
write.csv(median_changes, 'Statistics/stream_median_changes.csv', quote = F, row.names = F)

##################################################################################################################################
# PLOTTING
##################################################################################################################################
# Pres/future datasets
change_water_temp = data.frame(Present = pred_data$pc_water_temp_predicted[(pred_data$Date == 'Present') & (pred_data$SSP == 370) & (pred_data$Site == 'UP')],
                               Future = pred_data$pc_water_temp_predicted[(pred_data$Date == 'Future') & (pred_data$SSP == 370) & (pred_data$Site == 'UP')])
change_turbidity = data.frame(Present = pred_data$pc_turbidity_predicted[(pred_data$Date == 'Present') & (pred_data$SSP == 370) & (pred_data$Site == 'UP')],
                              Future = pred_data$pc_turbidity_predicted[(pred_data$Date == 'Future') & (pred_data$SSP == 370) & (pred_data$Site == 'UP')])
change_chla = data.frame(Present = pred_data$chla_predicted[(pred_data$Date == 'Present') & (pred_data$SSP == 370) & (pred_data$Site == 'UP')],
                         Future = pred_data$chla_predicted[(pred_data$Date == 'Future') & (pred_data$SSP == 370) & (pred_data$Site == 'UP')])

#################################################################
# 1a. Plots present/future for env params
plot_env_data = pred_data %>% group_by(Date, SSP) %>% 
  select(pc_water_temp_predicted, pc_turbidity_predicted, pc_conductivity_predicted,
         pc_ph_predicted, nut_din_predicted, nut_srp_predicted, chla_predicted, bacterial_abundance_predicted) %>% 
  melt(measure.vars = c('pc_water_temp_predicted', 'pc_turbidity_predicted', 'pc_conductivity_predicted',
                        'pc_ph_predicted', 'nut_din_predicted', 'nut_srp_predicted', 'chla_predicted', 'bacterial_abundance_predicted'))
plot_env_data$SSP[plot_env_data$SSP == 126] = 'RCP 2.6'
plot_env_data$SSP[plot_env_data$SSP == 370] = 'RCP 4.5'
plot_env_data$SSP[plot_env_data$SSP == 585] = 'RCP 8.5'

plot_env_data$Date[plot_env_data$Date == 'Future'] = 'Average for the 2070-2100 period'
plot_env_data$Date[plot_env_data$Date == 'Present'] = 'Sampling year (2019-2022)'

plot_distribution <- function(data, resp_var, x_lab, strip = T){
  p = ggplot(data[data$variable == resp_var,], aes(x=value, y=Date, fill=Date)) + 
    geom_density_ridges(quantile_lines = TRUE, quantiles = 2, alpha=0.75) + theme_linedraw() + scale_fill_manual(values =c('#2F5494', '#FA726B')) + facet_grid(~SSP) +
    ylab('') + xlab(x_lab) + theme(axis.text.y = element_blank(),
                                   legend.title = element_blank(),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   panel.spacing = unit(1.5, "lines"),
                                   plot.margin = unit(c(0.25,0.25,0.5,0.25), units = "lines" ))
  if (strip == F){p = p + theme(strip.text = element_blank(), 
                                strip.background = element_blank())}
  return(p)}

p1 = plot_distribution(plot_env_data, 'pc_water_temp_predicted', ""*Water~temperature~(ln~degree*C)*"", T)
p2 = plot_distribution(plot_env_data, 'pc_turbidity_predicted', ""*Turbidity~(ln~NTU)*"", T)
p3 = plot_distribution(plot_env_data, 'pc_conductivity_predicted', ""*Conductivity~(ln~mu*S~cm^-1)*"", T)
p4 = plot_distribution(plot_env_data, 'pc_ph_predicted', "pH", T)
p5 = plot_distribution(plot_env_data, 'nut_din_predicted', ""*Dissolved~inorganic~nitrogen~(ln~mu*g~L^-1)*"", T)
p6 = plot_distribution(plot_env_data, 'nut_srp_predicted', ""*Soluble~reactive~phosphate~(ln~mu*g~L^-1)*"", T)
p7 = plot_distribution(plot_env_data, 'chla_predicted', ""*Chlorophyll-italic(a)~(ln~mu*g~g^-1)*"", T)
p8 = plot_distribution(plot_env_data, 'bacterial_abundance_predicted', ""*Bacterial~abundance~(ln~cells~g^-1)*"", T)
p = ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, common.legend = T, ncol = 1) + theme(plot.margin = margin(0.1,0.5,0.2,0, "cm")) 
ggsave(p, filename = 'Plots/SFig_2_env_scenarios.pdf', width = 6, height = 11)


# 1b. Plots present/future for env params
plot_env_present = pred_data %>% filter(Date == 'Present', Site == 'UP') %>% group_by(SSP, Mountain_range) %>%
  select(pc_water_temp_predicted, pc_turbidity_predicted, pc_conductivity_predicted, pc_ph_predicted, 
         nut_din_predicted, nut_srp_predicted, Mountain_range, SSP) %>% 
  melt(measure.vars = c('pc_water_temp_predicted', 'pc_turbidity_predicted', 'pc_conductivity_predicted', 
                        'pc_ph_predicted', 'nut_din_predicted', 'nut_srp_predicted'), value.name = 'Present', id.vars = c('SSP', 'Mountain_range'))
plot_env_future = pred_data %>% filter(Date == 'Future', Site == 'UP') %>% group_by(SSP, Mountain_range)  %>%
  select(pc_water_temp_predicted, pc_turbidity_predicted, pc_conductivity_predicted, pc_ph_predicted, 
         nut_din_predicted, nut_srp_predicted, Mountain_range, SSP) %>% 
  melt(measure.vars = c('pc_water_temp_predicted', 'pc_turbidity_predicted', 'pc_conductivity_predicted', 
                        'pc_ph_predicted', 'nut_din_predicted', 'nut_srp_predicted'), value.name = 'Future', id.vars = c('SSP', 'Mountain_range'))
plot_env_preds = plot_env_present
colnames(plot_env_preds) = c('SSP', 'Mountain_range', 'variable', 'Present')
plot_env_preds$Future = plot_env_future$Future

plot_env_preds$SSP[plot_env_preds$SSP == 126] = 'RCP 2.6'
plot_env_preds$SSP[plot_env_preds$SSP == 370] = 'RCP 4.5'
plot_env_preds$SSP[plot_env_preds$SSP == 585] = 'RCP 8.5'

labels_facets <- list(
  'pc_water_temp'=""*A.~~~Water~temperature~(ln~degree*C)*"",
  'pc_turbidity'= ""*B.~~~Turbidity~(ln~NTU)*"",
  'pc_conductivity'=""*C.~~~Conductivity~(ln~mu*S~cm^-1)*"",
  'pc_ph'="D.   pH",
  'nut_din'=""*E.~~~Dissolved~inorganic~nitrogen~(ln~mu*g~L^-1)*"",
  'nut_srp'=""*F.~~~Soluble~reactive~phosphate~(ln~mu*g~L^-1)*"")
facets_labeller <- function(variable,value){
  return(labels_facets[value])}

p = ggplot(plot_env_preds[plot_env_preds$SSP == 'RCP 4.5',], aes(x=Present, y=Future)) + 
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
ggsave(plot = p, 'Plots/Fig_2_env_changes.pdf', width = 5.75, height = 8.5)

#################################################################
# 2. Plots for chlorophyll
chla_turb = data.frame(Turbidity = c(pred_data$pc_turbidity_predicted[(pred_data$Date == 'Present') & (pred_data$Site == 'UP')],
                                     pred_data$pc_turbidity_predicted[(pred_data$Date == 'Future') & (pred_data$Site == 'UP')]),
                       Bacterial_abundance = c(pred_data$bacterial_abundance_predicted[(pred_data$Date == 'Present') & (pred_data$Site == 'UP')],
                                               pred_data$bacterial_abundance_predicted[(pred_data$Date == 'Future') & (pred_data$Site == 'UP')]),
                       Chlorophyll = c(pred_data$chla_predicted[(pred_data$Date == 'Present') & (pred_data$Site == 'UP')],
                                       pred_data$chla_predicted[(pred_data$Date == 'Future') & (pred_data$Site == 'UP')]),
                       Shannon = c(pred_data$shannon_predicted[(pred_data$Date == 'Present') & (pred_data$Site == 'UP')],
                                   pred_data$shannon_predicted[(pred_data$Date == 'Future') & (pred_data$Site == 'UP')]),
                       SSP = c(pred_data$SSP[(pred_data$Date == 'Present') & (pred_data$Site == 'UP')],
                               pred_data$SSP[(pred_data$Date == 'Future') & (pred_data$Site == 'UP')]),
                       Date = c(rep('Sampling year (2019-2022)', length(pred_data$Date[(pred_data$Date == 'Present') & (pred_data$Site == 'UP')])),
                                rep('Average for the 2070-2100 period', length(pred_data$Date[(pred_data$Date == 'Future') & (pred_data$Site == 'UP')]))),
                       Sample = c(pred_data$Sample[(pred_data$Date == 'Present') & (pred_data$Site == 'UP')],
                                  pred_data$Sample[(pred_data$Date == 'Future') & (pred_data$Site == 'UP')]))
chla_turb$Bacterial_abundance = log10(exp(chla_turb$Bacterial_abundance))

p1 = ggplot(chla_turb[chla_turb$SSP == 370,], aes(x=Turbidity, y=Chlorophyll, colour=Date)) + 
  geom_line(aes(group=Sample), size=0.5, colour='black', alpha=0.1) + geom_point(alpha=0.7, size=3) + 
  #geom_smooth(chla_turb[(chla_turb$SSP == 370) & (chla_turb$Date == 'Sampling year (2019-2022)'),], mapping = aes(x=Turbidity, y=Chlorophyll), colour='#FA726B', size=1.5, method='lm', se=F) + 
  #geom_smooth(chla_turb[chla_turb$SSP == 370,], mapping = aes(x=Turbidity, y=Chlorophyll), colour='#2F5494', size=1.5, method='lm', se=F) + 
  theme_bw() + scale_colour_manual(values=c('#15356B', '#E04D55')) +
  labs(x = ""*Turbidity~(ln~NTU)*"",
       y = ""*Chlorophyll-italic(a)~(ln~mu*g~g^-1)*"") + theme(legend.title = element_blank(),
                                                               panel.grid.major = element_blank(),
                                                               panel.grid.minor = element_blank())

p2 = ggplot(chla_turb[chla_turb$SSP == 370,], aes(x=Bacterial_abundance, y=Chlorophyll, colour=Date)) + 
  geom_line(aes(group=Sample), size=0.5, colour='black', alpha=0.1) + geom_point(alpha=0.7, size=3) + 
  #geom_smooth(chla_turb[(chla_turb$SSP == 370) & (chla_turb$Date == 'Sampling year (2019-2022)'),], mapping = aes(x=Water_temp, y=Chlorophyll), colour='#FA726B', size=1.5, method='lm', se=F) + 
  #geom_smooth(chla_turb[chla_turb$SSP == 370,], mapping = aes(x=Water_temp, y=Chlorophyll), colour='#2F5494', size=1.5, method='lm', se=F) + 
  theme_bw() + scale_colour_manual(values=c('#15356B', '#E04D55')) +
  labs(x = ""*Bacterial~abundance~(log[10]~cells~g^-1)*"",
       y = "") + theme(legend.title = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())

p3 = ggplot(chla_turb[chla_turb$SSP == 370,], aes(x=Shannon, y=Chlorophyll, colour=Date)) + 
  geom_line(aes(group=Sample), size=0.5, colour='black', alpha=0.1) + geom_point(alpha=0.7, size=3) + 
  #geom_smooth(chla_turb[(chla_turb$SSP == 370) & (chla_turb$Date == 'Sampling year (2019-2022)'),], mapping = aes(x=Water_temp, y=Chlorophyll), colour='#FA726B', size=1.5, method='lm', se=F) + 
  #geom_smooth(chla_turb[chla_turb$SSP == 370,], mapping = aes(x=Water_temp, y=Chlorophyll), colour='#2F5494', size=1.5, method='lm', se=F) + 
  theme_bw() + scale_colour_manual(values=c('#15356B', '#E04D55')) +
  labs(x = ""*Shannon~index*"",
       y = "") + theme(legend.title = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())

p = ggarrange(p1, p2, p3, labels = c('A','B', 'C'), ncol = 3, nrow = 1, common.legend = T, align = 'h') 
ggsave(p, filename = 'Plots/Fig_3_chla_prediction.pdf', width = 10, height = 4)

# Turbidity
cor.test(chla_turb$Chlorophyll[(chla_turb$Date == 'Average for the 2070-2100 period') & (chla_turb$SSP == 126)] - chla_turb$Chlorophyll[(chla_turb$Date == 'Sampling year (2019-2022)') & (chla_turb$SSP == 126)],
         chla_turb$Turbidity[(chla_turb$Date == 'Average for the 2070-2100 period') & (chla_turb$SSP == 370)] - chla_turb$Turbidity[(chla_turb$Date == 'Sampling year (2019-2022)') & (chla_turb$SSP == 370)], method='spearman')
cor.test(chla_turb$Chlorophyll[(chla_turb$Date == 'Average for the 2070-2100 period') & (chla_turb$SSP == 585)] - chla_turb$Chlorophyll[(chla_turb$Date == 'Sampling year (2019-2022)') & (chla_turb$SSP == 585)],
         chla_turb$Turbidity[(chla_turb$Date == 'Average for the 2070-2100 period') & (chla_turb$SSP == 126)] - chla_turb$Turbidity[(chla_turb$Date == 'Sampling year (2019-2022)') & (chla_turb$SSP == 126)], method='spearman')
cor.test(chla_turb$Chlorophyll[(chla_turb$Date == 'Average for the 2070-2100 period') & (chla_turb$SSP == 370)] - chla_turb$Chlorophyll[(chla_turb$Date == 'Sampling year (2019-2022)') & (chla_turb$SSP == 370)],
         chla_turb$Turbidity[(chla_turb$Date == 'Average for the 2070-2100 period') & (chla_turb$SSP == 585)] - chla_turb$Turbidity[(chla_turb$Date == 'Sampling year (2019-2022)') & (chla_turb$SSP == 585)], method='spearman')
# 126: S = 1330870, p-value < 2.2e-16, -0.843913 
# 370: S = 1360302, p-value < 2.2e-16, -0.8846908 
# 585: S = 1386768, p-value < 2.2e-16, -0.9213593 

# Bacterial abundance
cor.test(chla_turb$Chlorophyll[(chla_turb$Date == 'Average for the 2070-2100 period') & (chla_turb$SSP == 126)] - chla_turb$Chlorophyll[(chla_turb$Date == 'Sampling year (2019-2022)') & (chla_turb$SSP == 126)],
         chla_turb$Bacterial_abundance[(chla_turb$Date == 'Average for the 2070-2100 period') & (chla_turb$SSP == 370)] - chla_turb$Bacterial_abundance[(chla_turb$Date == 'Sampling year (2019-2022)') & (chla_turb$SSP == 370)], method='spearman')
cor.test(chla_turb$Chlorophyll[(chla_turb$Date == 'Average for the 2070-2100 period') & (chla_turb$SSP == 585)] - chla_turb$Chlorophyll[(chla_turb$Date == 'Sampling year (2019-2022)') & (chla_turb$SSP == 585)],
         chla_turb$Bacterial_abundance[(chla_turb$Date == 'Average for the 2070-2100 period') & (chla_turb$SSP == 126)] - chla_turb$Bacterial_abundance[(chla_turb$Date == 'Sampling year (2019-2022)') & (chla_turb$SSP == 126)], method='spearman')
cor.test(chla_turb$Chlorophyll[(chla_turb$Date == 'Average for the 2070-2100 period') & (chla_turb$SSP == 370)] - chla_turb$Chlorophyll[(chla_turb$Date == 'Sampling year (2019-2022)') & (chla_turb$SSP == 370)],
         chla_turb$Bacterial_abundance[(chla_turb$Date == 'Average for the 2070-2100 period') & (chla_turb$SSP == 585)] - chla_turb$Bacterial_abundance[(chla_turb$Date == 'Sampling year (2019-2022)') & (chla_turb$SSP == 585)], method='spearman')
# 126: S = 105936, p-value < 2.2e-16, 0.8532263 
# 370: S = 143870, p-value < 2.2e-16, 0.8006689 
# 585: S = 66400 , p-value < 2.2e-16, 0.9080032

#################################################################
# 3. Plots for feature importances
models_list = list('pc_water_temp126'=water_temp_126,
                   'pc_water_temp370'=water_temp_370, 
                   'pc_water_temp585'=water_temp_585, 
                   'pc_turbidity126'=turbidity_126, 
                   'pc_turbidity370'=turbidity_370, 
                   'pc_turbidity585'=turbidity_585, 
                   'pc_conductivity126'=conductivity_126, 
                   'pc_conductivity370'=conductivity_370, 
                   'pc_conductivity585'=conductivity_585, 
                   'pc_ph126'=ph_126, 
                   'pc_ph370'=ph_370, 
                   'pc_ph585'=ph_585,
                   'nut_din126'=din_126, 
                   'nut_din370'=din_370, 
                   'nut_din585'=din_585,
                   'nut_srp126'=srp_126, 
                   'nut_srp370'=srp_370, 
                   'nut_srp585'=srp_585,
                   'chla126'=chla_126,
                   'chla370'=chla_370,
                   'chla585'=chla_585,
                   'bacterial_abundance126'=sba_126,
                   'bacterial_abundance370'=sba_370,
                   'bacterial_abundance585'=sba_585)

registerDoMC(8)
smooth_estimates_full = foreach(var = c('pc_water_temp', 'pc_turbidity', 'pc_conductivity', 
                                        'pc_ph', 'nut_din', 'nut_srp', 'chla', 'bacterial_abundance'), .combine = bind_rows) %dopar%{
                                          smooth_estimates_var = data.frame()
                                          print(var)
                                          for (scenario in c(126, 370, 585)){
                                            print(scenario)
                                            models = models_list[[paste0(var, scenario, collapse = '')]]$models
                                            
                                            for (i in 1:90){
                                              model = models[[i]]
                                              smooth_est = smooth_estimates(model)
                                              smooth_est = smooth_est[smooth_est$smooth != 's(latitude,longitude)',]
                                              smooth_est$Variable = var
                                              smooth_est$Scenario = scenario
                                              smooth_est$model_i = i
                                              smooth_estimates_var = bind_rows(smooth_estimates_var, smooth_est)}}
                                          return(smooth_estimates_var)}

smooth_estimates_plot = smooth_estimates_full %>% 
  melt(measure.vars = c('gl_dist', 'gl_area', 'gl_coverage', 'clim_tas', 'clim_scd', 'clim_pr', 
                        'min_calcite', 'min_clays', 'min_feldspar', 'min_quartz'),
       id.vars = c('smooth', 'est', 'se', 'Variable', 'Scenario','model_i'),
       variable.name = 'smooth_var', value.name = 'smooth_value') %>%
  select(-smooth) %>% filter(!is.na(smooth_value)) %>% 
  group_by(Variable, Scenario, smooth_var, smooth_value) %>% 
  summarise(est = median(est), se = median(se))

smooth_estimates_plot$Scenario[smooth_estimates_plot$Scenario == 126] = 'RCP 2.6'
smooth_estimates_plot$Scenario[smooth_estimates_plot$Scenario == 370] = 'RCP 4.5'
smooth_estimates_plot$Scenario[smooth_estimates_plot$Scenario == 585] = 'RCP 8.5'

p1 = ggplot(smooth_estimates_plot, aes(colour=Scenario)) + facet_grid(Variable~smooth_var, scales = 'free') + 
  geom_line(stat="smooth", aes(x=as.numeric(smooth_value), y=est), method='loess', se=F) + 
  geom_line(stat="smooth", aes(x=as.numeric(smooth_value), y=est-se), method='loess', linetype='dashed', se=F) + 
  geom_line(stat="smooth", aes(x=as.numeric(smooth_value), y=est+se), method='loess', linetype='dashed', se=F) +
  theme_linedraw() + theme(panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.text.y = element_text(size = 7)) +
  ylab('Estimate') + xlab('Smooth value') + scale_colour_manual(values = c('#AB79F5', '#1E81FC', '#F58965'))

smooth_estimates_n_feats = smooth_estimates_full %>% 
  melt(measure.vars = c('gl_dist', 'gl_area', 'gl_coverage', 'clim_tas', 'clim_scd', 'clim_pr', 
                        'min_calcite', 'min_clays', 'min_feldspar', 'min_quartz'),
       id.vars = c('smooth', 'est', 'se', 'Variable', 'Scenario','model_i'),
       variable.name = 'smooth_var', value.name = 'smooth_value') %>% filter(!is.na(smooth_value))%>%
  select(-smooth, -est, -se, -smooth_value)  %>% 
  group_by(Variable, Scenario, smooth_var) %>% distinct() %>%
  summarise(prop_selected = n()/90) %>% mutate(Scenario = as.factor(Scenario))
p2 = ggplot(smooth_estimates_n_feats, aes(y=Scenario, x='1', fill = Scenario, alpha=prop_selected, label=round(prop_selected, 3))) + 
  facet_grid(Variable~smooth_var, scales = 'free', space = 'free')+
  geom_tile() + geom_text() + xlab('') + ylab('') +
  theme_linedraw() + theme(panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.text.y = element_text(size = 7),
                           axis.ticks.x = element_blank(),
                           axis.text.x = element_blank(),
                           legend.position = 'none') + scale_fill_manual(values = c('#AB79F5', '#1E81FC', '#F58965'))
p = ggarrange(p1, p2, ncol = 1, nrow = 2, common.legend = T)
ggsave(plot = p, filename = 'Plots/SFig_3_models_variables.pdf', width = 10, height = 20)



