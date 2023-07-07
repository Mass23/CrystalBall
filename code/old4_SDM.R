##################### LOAD PACKAGES ##################### 
library(mgcv)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(purrr)
library(reshape2)
library(matrixStats)
library(foreach)
library(doMC)

setwd('~/Desktop/CrystalBall/code')
set.seed(23)

##################### FUNCTIONS ##################### 
PredictModelList <- function(models, new_data, seed){
  set.seed(seed) 
  preds = data.frame()
  for (i in 1:length(models)){
    pred = predict.gam(models[[i]], newdata = new_data, newdata.guaranteed = T, type = "response")
    preds = rbind(preds, pred)}
  preds = as.data.frame(preds)
  return(vapply(1:ncol(preds), function(i) median(preds[,i]), FUN.VALUE = numeric(1)))}

RemoveCorrelatedVars <- function(top_50_comb, data){
  all_vars = unique(unlist(unlist(top_50_comb %>% select(-median_log_p))))
  all_corrs = cor(data %>% select(all_vars))
  all_corrs[rownames(all_corrs) == colnames(all_corrs)] = 0
  
  idx_to_remove = c()
  for (i in 1:50){
    vars = top_50_comb[i, colnames(top_50_comb) %in% c('V1', 'V2', 'V3')] %>% unlist() %>% unname()
    max_corr = max(all_corrs[rownames(all_corrs) %in% vars, colnames(all_corrs) %in% vars])
    if (max_corr >= 0.7){idx_to_remove = c(idx_to_remove, i)}}
  if (length(idx_to_remove) > 0){top_50_comb = top_50_comb[-idx_to_remove,]}
  return(top_50_comb)}

SelectVars <- function(data, variables, seed){
  set.seed(seed) 
  vars_selection = data.frame(var=variables)
  vars_selection$p = map_dbl(vars_selection$var, function(x) summary(gam(data = data, 
                                                                         formula = eval(parse(text=paste0('log_abundance ~ s(', x,", k=3, bs='ts')"))), 
                                                                         family = gaussian()))$s.pv)
  
  vars_cor = as.data.frame(t(combn(vars_selection$var,3)))
  vars_cor$median_log_p = map_dbl(1:nrow(vars_cor), function(i) median(-log(vars_selection$p[vars_selection$var %in% vars_cor[i,]])))
  
  top_50_comb = vars_cor %>% top_n(50, median_log_p)
  top_combs = RemoveCorrelatedVars(top_50_comb, data)
  
  return(top_combs %>% top_n(1, median_log_p) %>% select(-median_log_p) %>% sample_n(1) %>% unlist())}

FitModel <- function(train_data, variables){
  formula = paste0('log_abundance ~ ', 's(', variables[1], ", k=3, bs='ts') + s(",
                                             variables[2], ", k=3, bs='ts') + s(",
                                             variables[3], ", k=3, bs='ts')", collapse = '')
  model = gam(data = train_data, formula = eval(parse(text=formula)), family = gaussian(), select = T)
  return(model)}

GAMCVKFoldFeatureSelection <- function(full_data, variables, kfold){
  set.seed(23) 
  t0 = Sys.time()
  cols_to_keep = c(variables, 'Glacier', 'log_abundance')
  present_data = full_data %>% filter(Date == 'Present') %>% na.omit()
  future_data = full_data %>% filter(Date == 'Future') %>% na.omit()
  present_data = present_data[sample(1:nrow(present_data)),]
  
  all_true_validate = c()
  all_true_train = c()
  all_pred_validate = c()
  all_pred_train = c()
  
  all_models = list()
  
  glaciers = unique(present_data$Glacier)
  groups = vapply(glaciers, function(x) sample(1:kfold,1), FUN.VALUE = numeric(1))
  present_data$fold = vapply(present_data$Glacier, function(x) groups[glaciers == x], FUN.VALUE = numeric(1))
  vars_selection_tab = c()
  
  for (kf in 1:kfold){
    print(kf)
    set.seed(kf) 
    
    train = present_data[present_data$fold != kf,] 
    validate = present_data[present_data$fold == kf,]
    train_folds = seq(1,kfold)[seq(1,kfold) != kf]
    
    vars_to_keep = SelectVars(train, variables, 23)
    vars_selection_tab = c(vars_selection_tab, vars_to_keep)
    
    model = FitModel(train, vars_to_keep)
    all_models[[length(all_models) + 1]] = model
    
    # Create predictions and gather true data
    k_true_validate = validate %>% pull(log_abundance) %>% as.vector()
    k_pred_validate = predict.gam(model, newdata = validate, newdata.guaranteed = T, type = "response")
    k_true_train = train %>% pull(log_abundance) %>% as.vector()
    k_pred_train = predict.gam(model, newdata = train, newdata.guaranteed = T, type = "response")
    
    # Append to the 
    all_true_validate = c(all_true_validate, k_true_validate)
    all_pred_validate = c(all_pred_validate, k_pred_validate)
    all_true_train = c(all_true_train, k_true_train)
    all_pred_train = c(all_pred_train, k_pred_train)}
  
  r2_train = cor(all_true_train, all_pred_train, use = "complete.obs", method = 'pearson')^2
  r2_validate = cor(all_true_validate, all_pred_validate, use = "complete.obs", method = 'pearson')^2
  varimp = table(vars_selection_tab)/sum(table(vars_selection_tab))*3
  
  print(paste0('Final model on train set r2: ', r2_train))
  print(paste0('Final model on validation set r2: ', r2_validate))
  print(varimp/3)
  print(Sys.time() - t0)
  return(list(models=all_models, r2_train=r2_train, r2_validate=r2_validate, var_imp=varimp))}

CreateResRow <- function(model, present_data, future_data, mag_data, mean_rel_ab, half_min_non_zero, mag_idx, scenario){
  # Make model predictions
  present_preds = PredictModelList(model$models, present_data[present_data$Site == 'UP',], 23)
  present_preds = exp(present_preds) - half_min_non_zero
  
  future_preds =  PredictModelList(model$models, future_data[future_data$Site == 'UP',], 23)
  future_preds = exp(future_preds) - half_min_non_zero
  
  # Append to data frame
  res = as.data.frame(model$var_imp)
  res$r2 = model$r2_validate
  res$MAG = rownames(mag_data)[mag_idx]
  res$mean_rel_ab = mean_rel_ab
  res$scenario = scenario
  
  res$median_present = median(present_preds)
  res$q25_present = quantile(present_preds, probs = 0.25)
  res$q75_present = quantile(present_preds, probs = 0.75)
  
  res$median_future = median(future_preds)
  res$q25_future = quantile(future_preds, probs = 0.25)
  res$q75_future = quantile(future_preds, probs = 0.75)
  
  res$median_change = median(future_preds - present_preds)
  res$q25_change = quantile(future_preds - present_preds, probs = 0.25)
  res$q75_change = quantile(future_preds - present_preds, probs = 0.75)
  
  res$median_rel_change = median((future_preds - present_preds)/present_preds)
  res$q25_rel_change = quantile(future_preds - present_preds, probs = 0.25)
  res$q75_rel_change = quantile(future_preds - present_preds, probs = 0.75)
  
  wilc_test = wilcox.test(future_preds - present_preds)
  res$wilcox_p = wilc_test$p.value
  res$wilcox_stat = wilc_test$statistic
  return(res)}

ProcessMag <- function(mag_idx, scenario, var_data_scaled, mag_data, variables){
  present_data = var_data_scaled %>% filter(SSP == scenario, Date == 'Present')
  future_data = var_data_scaled %>% filter(SSP == scenario, Date == 'Future')
  
  train_data = present_data[present_data$Sample %in% colnames(mag_data),]
  train_data$abundance = map_dbl(train_data$Sample, function(x) ifelse(x %in% colnames(mag_data), 
                                                                       floor(mag_data[mag_idx,x]), 
                                                                       -999))
  train_data$abundance[train_data$abundance == -999] = NA
  
  # Abundance data transformation
  half_min_non_zero = min(train_data$abundance[train_data$abundance > 0]) / 2
  train_data$log_abundance = log(train_data$abundance + half_min_non_zero)
  mean_rel_ab = mean(as.double(mag_data[mag_idx,]) / colSums(mag_data))
  
  # Create model
  set.seed(mag_idx+scenario)
  model = GAMCVKFoldFeatureSelection(train_data, variables, 10)
  res = CreateResRow(model, present_data, future_data, mag_data, mean_rel_ab, half_min_non_zero, mag_idx, scenario)
  return(res)}

Dataloader <- function(){
  var_data = read.csv('../data/processed/projection_data_3_ssps.csv')
  var_data = var_data %>% filter(Mountain_range != 'Alaska Range')
  mag_data = read.csv('../data/processed/MAGs_cov_filtered.tsv', row.names = 'X')
  selected_variables = c("bioclim_PC1","bioclim_PC2","bioclim_PC3","bioclim_PC4","bioclim_PC5","bioclim_PC6",
                         "clim_tas","clim_pr","clim_scd",
                         "pc_water_temp_predicted","pc_turbidity_predicted","pc_conductivity_predicted","pc_ph_predicted",
                         "nut_din_predicted", "nut_srp_predicted", "chla_predicted",
                         "gl_dist","gl_coverage","gl_area","gl_index",
                         "min_feldspar", "min_quartz", "min_calcite", "min_clays",
                         "elevation","latitude","longitude",'clim_tasmin','clim_tasmax','clim_fcf','clim_fgd')
  loaded_data = list(var_data = var_data, mag_data = mag_data, selected_variables = selected_variables)
  return(loaded_data)}

##################### MAIN ##################### 
Main4 <- function(loaded_data = NULL){
  if (is.null(loaded_data)){loaded_data = Dataloader()}
  
  var_table = loaded_data$var_data
  mag_table = loaded_data$mag_data
  variables = loaded_data$selected_variables
  
  var_table_scaled = var_table %>% mutate_at(variables, funs(c(scale(.))))
  #registerDoMC(cores = 10)
  
  #1:nrow(mag_table) .errorhandling = 'remove', 
  set.seed(23)
  model_out = foreach(mag = 1:nrow(mag_table), .combine = rbind) %:% # for real: nrow(mag_table)
    foreach(scenario = c(126, 370, 585), .combine = rbind) %do% 
    ProcessMag(mag, scenario, var_table_scaled, mag_table, variables)
  write.csv(model_out, '../data/results/model_res_all_scenarios_final.csv', quote = F, row.names = F)}
  
