library(dplyr)
library(mgcv)
library(reshape2)
library(matrixStats)
library(purrr)
library(glmnet)
library(doMC)
library(foreach)

####### FUNCTIONS ################################### 
PredictModelList <- function(models, new_data, seed){
  set.seed(seed) 
  preds = data.frame(row.names = as.character(1:nrow(new_data)))
  for (i in 1:length(models)){
    pred = predict.gam(models[[i]], newdata = new_data, newdata.guaranteed = T, type = 'response')
    preds = rbind(preds, pred)}
  preds = t(as.matrix(preds))
  return(preds)}
  
StackingGLM <- function(pred_data, target_data){
  stack_model_lasso = cv.glmnet(x = pred_data, y=target_data, family = gaussian(), alpha=1, nfolds=5, type.measure="mse")
  stack_model_elnet = cv.glmnet(x = pred_data, y=target_data, family = gaussian(), alpha=0.5, nfolds=5, type.measure="mse")
  stack_model_ridge = cv.glmnet(x = pred_data, y=target_data, family = gaussian(), alpha=0.5, nfolds=5, type.measure="mse")
  error_lasso = stack_model_lasso$cvm[stack_model_lasso$lambda == stack_model_lasso$lambda.min]
  error_elnet = stack_model_elnet$cvm[stack_model_elnet$lambda == stack_model_elnet$lambda.min]
  error_ridge = stack_model_ridge$cvm[stack_model_ridge$lambda == stack_model_ridge$lambda.min]
  errors = c(error_lasso, error_elnet, error_ridge)
  if (min(errors) == errors[1]) {return(stack_model_lasso)}
  if (min(errors) == errors[2]) {return(stack_model_elnet)}
  if (min(errors) == errors[3]) {return(stack_model_ridge)}}
  
PredictStackingGLM <- function(models, train_data, validate_data, train_target, seed){
  train_preds = PredictModelList(models, train_data, seed)
  stack_glm = StackingGLM(train_preds, train_target)

  test_preds = PredictModelList(models, validate_data, seed)
  stack_preds = predict(stack_glm, newx = test_preds, type = 'response', s = stack_glm$lambda.min)
  return(list(predictions=as.vector(stack_preds[,1]), model=stack_glm))}

FinalPredictionsStack <- function(submodels, stack_models, new_data, kfold){
  all_stack_preds = data.frame(row.names = as.character(1:nrow(new_data)))
  indices = split(seq_along(submodels), ceiling(seq_along(submodels)/(kfold-1)))
  for (i in 1:kfold){
    i_models = submodels[indices[[i]]]
    i_models_preds = PredictModelList(i_models, new_data, i)
    stack_model_preds = as.vector(predict(stack_models[[i]], newx = i_models_preds, type = 'response', s = stack_models[[i]]$lambda.min))
    all_stack_preds = rbind(all_stack_preds, stack_model_preds)}
  all_stack_preds = as.data.frame(all_stack_preds)
  return(vapply(1:ncol(all_stack_preds), function(i) mean(all_stack_preds[,i]), FUN.VALUE = numeric(1)))}

SelectFeatures <- function(train, features, resp_var, seed){
  set.seed(seed)
  vars_selection = data.frame()
  for (feat in features){
    form_feats = paste0(resp_var, " ~ s(latitude, longitude, bs='sos', m=1, k=-1) + s(", feat,", k=3, bs='ts')", collapse = '')
    model_selection = bam(data = train, formula = eval(parse(text=form_feats)))
    summary_p = summary(model_selection)$s.pv[2]
    vars_selection = rbind(vars_selection, data.frame(Feature=feat, log_p=-log(summary_p)))}
  feat_selected = vars_selection %>% top_n(3, log_p) %>% pull(Feature)
  return(feat_selected)}

CreateModel <- function(train, features, resp_var, n){
  set.seed(n)
  
  feat_selected = SelectFeatures(train, features, resp_var, n)

  form_final = paste0(resp_var," ~ s(latitude, longitude, bs='sos', m=1, k=-1) + s(",
                      feat_selected[1], ", k=3, bs='ts') + s(", 
                      feat_selected[2], ", k=3, bs='ts') + s(", 
                      feat_selected[3], ", k=3, bs='ts')",  collapse = '')
  
  # Randomly select one sample for each GFS
  train_n = train[train$fold != n,]
  train_n = train_n[!duplicated(train_n$Glacier),]
  
  model_n = bam(data = train_n, formula = eval(parse(text=form_final)), select = T)
  return(model_n)}

GAMCVKFoldFeatureSelection <- function(full_data, resp_var, features, kfold){
  set.seed(23)
  registerDoMC(10)
  
  cols_to_keep = c(resp_var, 'Date', 'Glacier', 'latitude', 'longitude', features)
  full_data = full_data %>% select(all_of(cols_to_keep)) %>% na.omit()
  present_data = full_data %>% filter(Date == 'Present')
  
  true_y_validate = c()
  true_y_train = c()
  pred_y_validate = c()
  pred_y_train = c()
  
  all_models = list()
  stack_models = list()
  all_r2_weights = c()
  
  glaciers = unique(present_data$Glacier)
  groups = vapply(glaciers, function(x) sample(1:kfold,1), FUN.VALUE = numeric(1))
  present_data$fold = vapply(present_data$Glacier, function(x) groups[glaciers == x], FUN.VALUE = numeric(1))
  
  for (i in 1:kfold){
    set.seed(i)
    print(i)
    present_data = present_data[sample(1:nrow(present_data)),]
    
    train = present_data[present_data$fold != i,]
    validate = present_data[present_data$fold == i,]
    train_folds = seq(1,kfold)[seq(1,kfold) != i]

    # Build the models
    i_models = foreach(n=train_folds) %do% CreateModel(train, features, resp_var, n)
    for (i in 1:length(i_models)){all_models[[length(all_models) + 1]] = i_models[[i]]}

    # Keep predicted/observed values for validation and training sets
    #pred_validate = PredictModelList(i_models, resp_var, train, validate, i)
    validation = PredictStackingGLM(i_models, train, validate, as.vector(train[resp_var] %>% pull()), i)
    pred_validate = validation$predictions
    stack_models[[length(stack_models) + 1]] = validation$model
    true_y_validate = c(true_y_validate, as.vector(validate[resp_var] %>% pull()))
    pred_y_validate = c(pred_y_validate, pred_validate)
    
    #pred_train = PredictModelList(i_models, resp_var, train, train, i)
    pred_train = PredictStackingGLM(i_models, train, train, as.vector(train[resp_var] %>% pull()), i)$predictions
    true_y_train = c(true_y_train, as.vector(train[resp_var] %>% pull()))
    pred_y_train = c(pred_y_train, pred_train)}
  
  r2_train = cor(true_y_train, pred_y_train, use = "complete.obs", method = 'pearson')^2
  r2_validate = cor(true_y_validate, pred_y_validate, use = "complete.obs", method = 'pearson')^2
  
  print(paste0('Train r2: ', r2_train))
  print(paste0('Validation r2: ', r2_validate))
  
  return(list(all_models=all_models, stack_models=stack_models, r2_train=r2_train, r2_validate=r2_validate, n_glaciers=length(unique(full_data$Glacier))))}

ProcessVariable <- function(data_126, data_370, data_585, resp_var){
  if (resp_var == 'pc_water_temp'){variables = c('gl_area', 'gl_distance','gl_coverage', 'clim_tas', 'clim_pr', 'clim_scd')}
  else {variables = c('gl_area', 'gl_distance','gl_coverage', 'clim_tas', 'clim_pr', 'clim_scd', 'min_calcite', 'min_clays', 'min_feldspar', 'min_quartz')}
  water_temp_126 = GAMCVKFoldFeatureSelection(data_126, resp_var, variables, 10)
  water_temp_370 = GAMCVKFoldFeatureSelection(data_370, resp_var, variables, 10)
  water_temp_585 = GAMCVKFoldFeatureSelection(data_585, resp_var, variables, 10)
  return(list(ssp126=water_temp_126, ssp370=water_temp_370, ssp585=water_temp_585))}

CreatePredDatasets <- function(data){
  up_data_present = data %>% filter(Site == 'UP', Date == 'Present') %>% distinct() %>% 
    select(-pc_water_temp, -pc_turbidity, -pc_conductivity, -pc_ph, -nut_din, -nut_srp, 
           -chla, -bacterial_abundance, -Shannon, -Pielou, -mntd, -mpd, -pd) %>% na.omit()
  up_data_future = data %>% filter(Site == 'UP', Date == 'Future') %>% distinct() %>% 
    select(-pc_water_temp, -pc_turbidity, -pc_conductivity, -pc_ph, -nut_din, -nut_srp, 
           -chla, -bacterial_abundance, -Shannon, -Pielou, -mntd, -mpd, -pd) %>% na.omit()
  return(rbind(up_data_present, up_data_future))}

MainF <- function(){
  all_data = read.csv('data/processed/all_current_clean_3_ssps.csv')
  data_126 = all_data[all_data$Scenario == 126,]
  data_370 = all_data[all_data$Scenario == 370,]
  data_585 = all_data[all_data$Scenario == 585,]
  
  # Others
  variables_to_model = c('pc_water_temp', 'pc_turbidity', 'pc_conductivity', 'pc_ph', 'nut_din', 'nut_srp', 'chla',
                         'bacterial_abundance', 'Shannon', 'Pielou', 'mntd', 'mpd')
  
  pred_data = CreatePredDatasets(all_data)
  
  stats_df = data.frame()
  for (variable in variables_to_model){
    print(variable)
    processed_var = ProcessVariable(data_126, data_370, data_585, variable)
    pred_data[,paste0(variable, '_predicted', collapse = '')] = NA
    
    for (scenario in c('ssp126', 'ssp370', 'ssp585')){
      submodels = processed_var[[scenario]]$all_models
      stack_models = processed_var[[scenario]]$stack_models
      
      # Gather stats
      r2_train = processed_var[[scenario]]$r2_train
      r2_validation = processed_var[[scenario]]$r2_validate
      n_glaciers = processed_var[[scenario]]$n_glaciers
      stats_var = data.frame(Variable=variable, Scenario=scenario, Train_r2=r2_train, Validation_r2=r2_validation, N=n_glaciers)
      stats_df = rbind(stats_df, stats_var)
      
      # Add predictions to the data
      predictions = FinalPredictionsStack(submodels, stack_models, pred_data, 10)
      scenario_var_models = list(sub=submodels, stack=stack_models, scenario=scenario, variable=variable)
      saveRDS(scenario_var_models, paste0('data/processed/',variable,'_',scenario,'.rds',collapse=''))
      pred_data[paste0('ssp', pred_data$Scenario) == scenario, paste0(variable, '_predicted')] = predictions[paste0('ssp', pred_data$Scenario) == scenario]}}
  
  pred_data = pred_data %>% arrange(Sample)
  write.csv(pred_data, 'data/processed/all_projections_3_ssps.csv', quote = F, row.names = F)
  write.csv(stats_df, 'stats/stream_parameters_performance.csv', quote = F, row.names = F)}

