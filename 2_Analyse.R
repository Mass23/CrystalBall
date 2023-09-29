Sys.setenv("LANGUAGE"="EN")

Main <- function(){
  #source('code/d_data_processing.R')
  #print('Running d_data_processing.R ...')
  #MainD()
  
  #source('code/e_data_exploration.R')
  #print('Running e_data_exploration.R')
  #MainE()
  
  #source('code/f_stream_parameters_models.R')
  #print('Running f_stream_parameters_models.R ...')
  #MainF()

  #source('code/g_stream_parameters_analyses.R')
  #print('Running g_stream_parameters_analyses.R ...')
  #MainG()

  source('code/h_stream_response_curves.R')
  print('Running h_stream_response_curves.R ...')
  MainH()

  #source('code/i_greening_analyses.R')
  #print('Running i_greening_analyses.R ...')
  #MainI()

  #source('code/j_strain_models.R')
  #print('Running j_strain_models.R ...')
  #MainJ()

  #source('code/k_strain_changes_analyses.R')
  #print('Running k_strain_changes_analyses.R ...')
  #MainK()

  #source('code/l_strain_drivers_analyses.R')
  #print('Running l_strain_drivers_analyses.R ...')
  #MainL()

  #source('code/m_functional_analysis.R')
  #print('Running m_functional_analysis.R ...')
  #MainM()
}

Main()
