source('code/d_data_processing.R')
source('code/e_data_exploration.R')
source('code/f_stream_parameters_models.R')
source('code/g_stream_parameters_analyses.R')


Main <- function(){
  #print('Running d_data_processing.R ...')
  #MainD()
  
  print('Running e_data_exploration.R')
  #MainE()
  
  print('Running f_stream_parameters_models.R ...')
  #MainF()

  print('Running g_stream_parameters_analyses.R ...')
  MainG()

  print('Running h_stream_response_curves.R ...')
  #MainH()

  print('Running i_greening_analyses.R ...')
  #MainI()


}

Main()
