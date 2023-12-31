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
library(glmnet)
source('code/f_stream_parameters_models.R')

set.seed(23)
variables = c("bioclim_PC1","bioclim_PC2","bioclim_PC3","bioclim_PC4","bioclim_PC5","bioclim_PC6",
                         "clim_tas","clim_pr","clim_scd",'clim_tasmax','clim_tasmin',
                         "pc_water_temp_predicted","pc_turbidity_predicted","pc_conductivity_predicted","pc_ph_predicted",
                         "nut_din_predicted", "nut_srp_predicted", "chla_predicted",
                         "gl_distance","gl_coverage","gl_area",
                         "min_feldspar", "min_quartz", "min_calcite", "min_clays",
                         "elevation","latitude","longitude","slope","abs_latitude")

##################### FUNCTIONS ##################### 
PredictModelListStrain <- function(models, new_data, seed){
  set.seed(seed) 
  preds = data.frame(row.names = as.character(1:nrow(new_data)))
  for (i in 1:length(models)){
    pred = predict.bam(models[[i]], newdata = new_data, newdata.guaranteed = T, type = 'response')
    preds = rbind(preds, pred)}
  preds = t(as.matrix(preds))
  return(preds)}
  
StackingGLMStrain <- function(pred_data, target_data){
  stack_model_ridge = cv.glmnet(x = pred_data, y=target_data, family = gaussian(), alpha=0, nfolds=5, type.measure="mse")
  {return(stack_model_ridge)}}
  
PredictStackingGLMStrain <- function(models, train_data, validate_data, train_target, seed){
  train_preds = PredictModelListStrain(models, train_data, seed)
  stack_glm = StackingGLMStrain(train_preds, train_target)

  test_preds = PredictModelListStrain(models, validate_data, seed)
  stack_preds = predict(stack_glm, newx = test_preds, type = 'response', s = stack_glm$lambda.min)
  return(list(predictions=as.vector(stack_preds[,1]), model=stack_glm))}
  
RemoveCorrelatedVars <- function(top_50_comb, data){
  all_vars = unique(unlist(unlist(top_50_comb %>% select(-median_log_p))))
  all_corrs = cor(data %>% select(all_vars))
  all_corrs[rownames(all_corrs) == colnames(all_corrs)] = 0
  
  idx_to_remove = c()
  for (i in 1:50){
    vars = top_50_comb[i, colnames(top_50_comb) %in% c('V1', 'V2', 'V3', 'V4')] %>% unlist() %>% unname()
    max_corr = max(all_corrs[rownames(all_corrs) %in% vars, colnames(all_corrs) %in% vars])
    if (max_corr >= 0.7){idx_to_remove = c(idx_to_remove, i)}}
  if (length(idx_to_remove) > 0){top_50_comb = top_50_comb[-idx_to_remove,]}
  return(top_50_comb)}

SelectVarsStrain <- function(data, variables, seed){
  set.seed(seed) 
  vars_selection = data.frame(var=variables)
  vars_selection$p = map_dbl(vars_selection$var, function(x) summary(bam(data = data, 
                                                                         formula = eval(parse(text=paste0('log_abundance ~ s(', x,", k=3, bs='ts')"))), 
                                                                         family = gaussian()))$s.pv)
  
  vars_cor = as.data.frame(t(combn(vars_selection$var,4)))
  vars_cor$median_log_p = map_dbl(1:nrow(vars_cor), function(i) median(-log(vars_selection$p[vars_selection$var %in% vars_cor[i,]])))
  
  top_50_comb = vars_cor %>% top_n(50, median_log_p)
  top_combs = RemoveCorrelatedVars(top_50_comb, data)
  
  return(top_combs %>% top_n(1, median_log_p) %>% select(-median_log_p) %>% sample_n(1) %>% unlist())}

FitModelStrain <- function(train, variables, n){
  set.seed(n)
  feat_selected = SelectVarsStrain(train, variables, n)
  form_final = paste0("log_abundance ~ s(", variables[1], ", k=3, bs='ts') + s(",
                                            variables[2], ", k=3, bs='ts') + s(",
                                            variables[3], ", k=3, bs='ts') + s(",
                                            variables[4], ", k=3, bs='ts')",  collapse = '')

  train_n = train[train$fold != n,]
  train_n = train_n[!duplicated(train_n$Glacier),]
  
  model_n = bam(data = train_n, formula = eval(parse(text=form_final)), select = T)
  return(model_n)}

FinalPredictionsStackStrain <- function(submodels, stack_models, new_data, kfold){
  all_stack_preds = data.frame(row.names = as.character(1:nrow(new_data)))
  indices = split(seq_along(submodels), ceiling(seq_along(submodels)/(kfold-1)))
  for (i in 1:kfold){
    i_models = submodels[indices[[i]]]
    i_models_preds = PredictModelList(i_models, new_data, i)
    stack_model_preds = as.vector(predict(stack_models[[i]], newx = i_models_preds, type = 'response', s = stack_models[[i]]$lambda.min))
    all_stack_preds = rbind(all_stack_preds, stack_model_preds)}
  all_stack_preds = as.data.frame(all_stack_preds)
  return(vapply(1:ncol(all_stack_preds), function(i) mean(all_stack_preds[,i]), FUN.VALUE = numeric(1)))}

StackedGAMCVKFoldFeatureSelectionStrain <- function(full_data, variables, kfold){
  set.seed(23) 
  t0 = Sys.time()
  cols_to_keep = c(variables, 'Glacier', 'log_abundance')
  present_data = full_data %>% filter(Date == 'Present') %>% na.omit()
  future_data = full_data %>% filter(Date == 'Future') %>% na.omit()
  
  true_y_validate = c()
  true_y_train = c()
  pred_y_validate = c()
  pred_y_train = c()
  
  all_models = list()
  stack_models = list()
  
  glaciers = unique(present_data$Glacier)
  groups = vapply(glaciers, function(x) sample(1:kfold,1), FUN.VALUE = numeric(1))
  present_data$fold = vapply(present_data$Glacier, function(x) groups[glaciers == x], FUN.VALUE = numeric(1))
  vars_selection_tab = c()
  
  for (kf in 1:kfold){
    set.seed(kf)
    present_data = present_data[sample(1:nrow(present_data)),]
    
    train = present_data[present_data$fold != kf,]
    validate = present_data[present_data$fold == kf,]
    train_folds = seq(1,kfold)[seq(1,kfold) != kf]
    
    vars_to_keep = SelectVarsStrain(train, variables, 23)
    vars_selection_tab = c(vars_selection_tab, vars_to_keep)
    
    # Build the models
    i_models = foreach(n=1:(kfold-1)) %do% FitModelStrain(train, vars_to_keep, n)
    for (i in 1:length(i_models)){all_models[[length(all_models) + 1]] = i_models[[i]]}

    # Keep predicted/observed values for validation and training sets
    validation = PredictStackingGLMStrain(i_models, train, validate, as.vector(train %>% pull(log_abundance)), kf)
    pred_validate = validation$predictions
    stack_models[[length(stack_models) + 1]] = validation$model
    true_y_validate = c(true_y_validate, as.vector(validate %>% pull(log_abundance)))
    pred_y_validate = c(pred_y_validate, pred_validate)
    
    pred_train = PredictStackingGLMStrain(i_models, train, train, as.vector(train %>% pull(log_abundance)), kf)$predictions
    true_y_train = c(true_y_train, as.vector(train %>% pull(log_abundance)))
    pred_y_train = c(pred_y_train, pred_train)}
  
  r2_train = cor(true_y_train, pred_y_train, use = "complete.obs", method = 'pearson')^2
  r2_validate = cor(true_y_validate, pred_y_validate, use = "complete.obs", method = 'pearson')^2
  varimp = table(vars_selection_tab)/sum(table(vars_selection_tab))*3
  
  print(paste0('Final model on train set r2: ', r2_train))
  print(paste0('Final model on validation set r2: ', r2_validate))
  print(Sys.time() - t0)
  return(list(sub_models=all_models, stack_models=stack_models, r2_train=r2_train, r2_validate=r2_validate, var_imp=varimp))}

CreateResRow <- function(model, present_data, future_data, mag_data, mean_rel_ab, half_min_non_zero, mag_idx, scenario){
  # Make model predictions
  present_preds = FinalPredictionsStackStrain(model$sub_models, model$stack_models, present_data[present_data$Site == 'UP',], 10)
  present_preds = exp(present_preds) - half_min_non_zero
  
  future_preds =  FinalPredictionsStackStrain(model$sub_models, model$stack_models, future_data[future_data$Site == 'UP',], 10)
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
  present_data = var_data_scaled %>% filter(Scenario == scenario, Date == 'Present', Site == 'UP')
  future_data = var_data_scaled %>% filter(Scenario == scenario, Date == 'Future', Site == 'UP')
  
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
  model = StackedGAMCVKFoldFeatureSelectionStrain(train_data, variables, 10)
  res = CreateResRow(model, present_data, future_data, mag_data, mean_rel_ab, half_min_non_zero, mag_idx, scenario)
  return(res)}

Dataloader <- function(){
  var_data = read.csv('data/processed/all_projections_3_ssps.csv')
  var_data = var_data %>% filter(Mountain_range != 'Alaska Range')
  var_data$abs_latitude = abs(var_data$latitude)
  mag_data = read.csv('data/processed/MAG_cov_filtered.tsv', row.names = 'X')
  loaded_data = list(var_data = var_data, mag_data = mag_data, selected_variables = variables)
  return(loaded_data)}

##################### MAIN ##################### 
MainJ <- function(loaded_data = NULL){
  if (is.null(loaded_data)){loaded_data = Dataloader()}
  
  var_table = loaded_data$var_data
  mag_table = loaded_data$mag_data
  variables = loaded_data$selected_variables
  
  registerDoMC(cores = 42)
  
  set.seed(23)
  model_out = foreach(mag = 1:nrow(mag_table), .combine = rbind) %:%
    foreach(scenario = c(126, 370, 585), .combine = rbind) %dopar% 
    ProcessMag(mag, scenario, var_table, mag_table, variables)
  write.csv(model_out, 'stats/model_res_all_scenarios_final.csv', quote = F, row.names = F)
  print(quantile(model_out %>% select(MAG, r2) %>% distinct() %>% pull(r2)))}
  
