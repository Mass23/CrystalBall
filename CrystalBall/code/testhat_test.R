setwd('~/Desktop/CrystalBall/code')
source('4_SDM.R')

library(testthat)
local_edition(3)

##################### TEST ##################### 
test_that("Test", {
  set.seed(23)
  test_var_data = data.frame(var1 = rep(sample(1:100, 100, replace = T), 6),
                             var2 = rep(sample(1:100, 100, replace = T), 6),
                             var3 = rep(seq(100,1,-1), 6),
                             var4 = rep(1:100, 6),
                             SSP = c(rep(126, 200), rep(370, 200), rep(585, 200)),
                             Date = rep(c(rep('Present', 100), rep('Future', 100)), 3),
                             Sample = rep(as.character(1:100), 6),
                             Glacier = rep(as.character(1:100), 6),
                             Site = rep('UP', 600))
  test_mag_data = matrix(ncol = 100, nrow = 0)
  test_mag_data = rbind(test_mag_data, seq(1,100,1))
  test_mag_data = rbind(test_mag_data, seq(100,1,-1))
  test_mag_data = rbind(test_mag_data, sample(1:100, 100, replace = T))
  test_mag_data = rbind(test_mag_data, sample(1:100, 100, replace = T))
  test_mag_data = as.data.frame(test_mag_data)
  colnames(test_mag_data) = as.character(1:100)
  
  test_variables = c('var1','var2','var3','var4')
  loaded_data = list(var_data = test_var_data, mag_data = test_mag_data, selected_variables = test_variables)
  expect_snapshot(Main(loaded_data))})
################################################ 

