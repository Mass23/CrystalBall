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

PredictionsStack <- function(submodels, stack_models, new_data, covariates, kfold){
  all_stack_preds = data.frame(row.names = as.character(1:nrow(new_data)))
  indices = split(seq_along(submodels), ceiling(seq_along(submodels)/(kfold-1)))
  for (i in 1:kfold){
    i_models = submodels[indices[[i]]]
    i_models_preds = PredictModelList(i_models, new_data, i)
    stack_model_preds = predict(stack_models[[i]], newx = i_models_preds, type = 'response', s = stack_models[[i]]$lambda.min, se.fit=T)
    stack_model_preds$Sample = paste0('S',1:nrow(stack_model_preds))
    all_stack_preds = rbind(all_stack_preds, stack_model_preds)}
  all_stack_preds = as.data.frame(all_stack_preds)
  return(all_stack_preds %>% group_by(all_of(c('Sample', covariates))) %>% summarise(pred = median(fit), se = median(se.fit)))}

CreatePredictionTable <- function(resp_var, data, sub_models, stack_models, kfold){
    present_up_data = data %>% filter(Date == 'Present', Site == 'UP')
    full_table = data.frame()
    for (i in 1:kfold){
        i_sub_models = sub_models[[seq((i-1)*10,i*10)]]
        stack_model = stack_models[[i]]
        
        # List variables
        i_vars = c()
        for (sub_model in 1:9){
            vars = colnames(i_sub_models[[sub_model]]$model)
            vars = vars[!(vars %in% c('latitude', 'longitude', resp_var))]
            i_vars = unique(c(i_vars, vars))}
        
        # Create data frame
        i_table = data.frame()
        for (var in i_vars){
            other_vars = i_vars[i_vars != var]
            i_var_table = data.frame(latitude=data$latitude, longitude=data$longitude)
            i_var_table[,var] = map_dbl(1:nrow(i_var_table), function(i) data %>% filter(latitude == i_var_table$latitude[i], 
                                                                                         longitude == i_var_table$longitude[i]) %>% pull(var))
            for (other_var %in% other_vars){i_var_table[,other_var] = data %>% pull(other_var) %>% median()}
            i_var_table$i_fold = i
            i_table = rbind(i_table, i_var_table)}
        full_table = rbind(full_table, i_table)}
    return(full_table)}

PlotResponseCurves <- function(data, variables, scenario){
    for (resp_var in variables){
        if (resp_var == 'pc_water_temp'){variables = c('gl_area', 'gl_distance','gl_coverage', 'clim_tas', 'clim_pr', 'clim_scd')}
        else {variables = c('gl_area', 'gl_distance','gl_coverage', 'clim_tas', 'clim_pr', 'clim_scd', 'min_calcite', 'min_clays', 'min_feldspar', 'min_quartz')}
        for (scenario in c('ssp126', 'ssp370', 'ssp585')){
            new_data = expand.grid()
            models_file = readRDS(paste0('data/processed/',resp_var,'_',scenario,'.rds', collapse=''))

            new_data = CreatePredictionTable(data %>% filter(Scenario == as.integer(strsplit(scenario, 'ssp')[[1]][2])), models_file$sub, models_file$stack, 10)
            predictions = PredictionsStack(models_file$sub, models_file$stack, new_data, 10)
            

        }


    }
}

MainG <- function(){
  data = read.csv('data/processed/all_projections_3_ssps.csv')
  PlotResponseCurves(data, variables)}