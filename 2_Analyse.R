source('code/d_data_processing.R')
source('code/e_data_exploration.R')
source('code/f_stream_parameters_models.R')
#source('code/g_stream_parameters_analyses.R')
#source('code/h_strain_models.R')
#source('code/i_strain_analyses.R')
#source('code/j_functional_analyses.R')

Main <- function(){
  print('Running d_data_processing.R ...')
  #MainD()
  
  print('Running e_data_exploration.R')
  #MainE()
  
  print('Running f_stream_parameters_models.R ...')
  MainF()
}

Main()
