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

variables =  c('pc_water_temp', 'pc_turbidity', 'pc_conductivity', 'pc_ph', 'nut_din', 'nut_srp', 
              'chla', 'bacterial_abundance', 'Shannon', 'Pielou', 'mntd', 'mpd')

se <- function(x) sd(x) / sqrt(length(x))

PredictionsStack <- function(submodels, stackmodels, new_data, covariates, kfold){
  all_stack_preds = data.frame()
  indices = split(seq_along(submodels), ceiling(seq_along(submodels)/(kfold-1)))
  for (i in 1:kfold){
    i_models = submodels[indices[[i]]]
    i_models_preds = PredictModelList(i_models, new_data, i)
    stack_model_preds = predict(stackmodels[[i]], newx = i_models_preds, type = 'response', s = stackmodels[[i]]$lambda.min, se.fit=T)
    new_data$fit = as.vector(stack_model_preds)
    all_stack_preds = rbind(all_stack_preds, new_data)}
  return(all_stack_preds %>% group_by_at(c('latitude', 'longitude', 'covariate', 'covariate_val', 'covariate_N', covariates)) %>% summarise(pred = median(fit), se = se(fit), var_median = median(var_median)))}

CreatePredictionTable <- function(resp_var, data, sub_models, stack_models, kfold){
    present_up_data = data %>% filter(Date == 'Present', Site == 'UP')
    lat_lon_table = data %>% select(latitude, longitude) %>% distinct()
    indices = split(seq_along(sub_models), ceiling(seq_along(sub_models)/(kfold-1)))

    # List all covariates
    all_vars = c()
    for (i in 1:kfold){
        i_sub_models = sub_models[indices[[i]]]
        stack_model = stack_models[[i]]
        
        # List covariates
        i_vars = c()
        for (sub_model in 1:9){
            vars = colnames(i_sub_models[[sub_model]]$model)
            vars = vars[!(vars %in% c('latitude', 'longitude', resp_var))]
            i_vars = c(i_vars, vars)}
        all_vars = c(all_vars, i_vars)}
    
    # Create the data frame
    full_table = data.frame()
    for (var in unique(all_vars)){
        other_vars = unique(all_vars)[unique(all_vars) != var]
        var_table = data.frame()
        min_val = present_up_data %>% pull(var) %>% min(na.rm=T)
        max_val = present_up_data %>% pull(var) %>% max(na.rm=T)
        for (value in seq(min_val, max_val, by = (max_val - min_val)/20)){
            value_tab = lat_lon_table
            value_tab[, var] = value
            var_table = rbind(var_table, value_tab)}
        for (other_var in other_vars){var_table[,other_var] = present_up_data %>% pull(other_var) %>% median()}
        var_table$var_median = present_up_data %>% pull(var) %>% median()
        var_table$covariate = var
        var_table$covariate_val = var_table %>% pull(var)
        var_table$covariate_N = mean(all_vars == var) * 3
        full_table = rbind(full_table, var_table)}
    return(as.data.frame(full_table))}

PlotResponseCurves <- function(data, variables, scenario){
    all_data = data.frame()
    for (resp_var in variables){
        print(resp_var)
        for (scenario in c('ssp126', 'ssp370', 'ssp585')){
            new_data = expand.grid()
            models_file = readRDS(paste0('data/processed/',resp_var,'_',scenario,'.rds', collapse=''))

            new_data = CreatePredictionTable(resp_var, 
                                             data %>% filter(Scenario == as.integer(strsplit(scenario, 'ssp')[[1]][2])), 
                                             models_file$sub, models_file$stack, 10)
            predictions = PredictionsStack(models_file$sub, models_file$stack, new_data, unique(new_data$covariate), 10)
            predictions$Variable = resp_var
            predictions$Scenario = scenario
            all_data = rbind(all_data, predictions)}}

    all_data = all_data %>% select(Variable, covariate, covariate_val, Scenario, var_median, pred, se, covariate_N) %>% 
                            group_by(Variable, covariate, covariate_val, Scenario) %>% 
                            summarise(pred = median(pred), se = median(se), var_median = median(var_median), 
                                      covariate_val = median(covariate_val), covariate_N=mean(covariate_N))
    write.csv(all_data, 'data/processed/stream_response_curves.csv', quote = F, row.names = F)
    return(all_data)}

CreatePlots <- function(){
    all_data = read.csv('data/processed/stream_response_curves.csv')
    print(unique(all_data$covariate_N))

    p = ggplot(all_data, aes(colour=Scenario, alpha=covariate_N/4*3)) + facet_grid(Variable~covariate, scales = 'free') + 
                geom_line(aes(x=as.numeric(covariate_val), y=pred)) + 
                geom_line(aes(x=as.numeric(covariate_val), y=pred-se), linetype='dashed') + 
                geom_line(aes(x=as.numeric(covariate_val), y=pred+se), linetype='dashed') +
                theme_linedraw() + theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank(),
                                         strip.text.y = element_text(size = 7)) +
                ylab('Response') + xlab('Smooth value') + scale_colour_manual(values = c('#AB79F5', '#1E81FC', '#F58965')) + 
                theme(axis.text.x = element_text(angle =45, hjust = 1)) + labs(alpha='Proportion of models') + scale_alpha(breaks=c(0.1,0.5,1.0),limits = c(0, 1))
ggsave(plot = p, filename = paste0('plots/Fig_S4_stream_model_response_curves.pdf'), width = 10, height = 11)
}

MainH <- function(){
  data = read.csv('data/processed/all_projections_3_ssps.csv')
  data = data %>% filter(Site == 'UP', Date == 'Present')
  resp_data = PlotResponseCurves(data, variables)
  CreatePlots()}