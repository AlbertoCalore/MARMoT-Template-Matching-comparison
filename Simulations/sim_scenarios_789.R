
# Libraries ---------------------------------------------------------------


rm(list=ls())
setwd("~/Server/Scenarios/")

if(system.file(package="nnet") == "") install.packages("nnet")
library(nnet)
source("../Functions/data_generation.R")


# Data --------------------------------------------------------------------


dati = readRDS("../Data preprocessing/Dati/dati_FINAL.rds")


# Parameters --------------------------------------------------------------


seed = 1234
n_sim = 500000
n_treat = 500
noise = 1
rm_left = ceiling(0.1 * n_treat)
rm_right = ceiling(0.05 * n_treat) 

treatment = "PRVDR_NUM"


# Other functions ---------------------------------------------------------


rm.gc = function(obj){
  rm(obj)
  gc()
}

rm.left = function(scenario, rm_left_trt){
  scenario_left = scenario
  scenario_left[["sim_data"]] = scenario_left[["sim_data"]][!scenario_left[["sim_data"]][, treatment] %in% rm_left_trt, ]
  scenario_left[["sim_data"]][, treatment] = droplevels(scenario_left[["sim_data"]][, treatment])
  rm_left_trt = paste0(treatment, rm_left_trt)
  scenario_left[["true_effects"]] = scenario_left[["true_effects"]][!names(scenario_left[["true_effects"]]) %in% rm_left_trt]
  return(scenario_left)
}

rm.right = function(scenario, rm_right_trt){
  scenario_right = scenario
  scenario_right[["sim_data"]] = scenario_right[["sim_data"]][!scenario_right[["sim_data"]][, treatment] %in% rm_right_trt, ]
  scenario_right[["sim_data"]][, treatment] = droplevels(scenario_right[["sim_data"]][, treatment])
  rm_right_trt = paste0(treatment, rm_right_trt)
  scenario_right[["true_effects"]] = scenario_right[["true_effects"]][!names(scenario_right[["true_effects"]]) %in% rm_right_trt]
  return(scenario_right)
}


# Simulations -------------------------------------------------------------


for (i in 1:10) {
  scenario = gen.data(dati, n_sim, n_treat, noise, i)
  path = paste0("../Scenarios/Scenario_7/seed", sep = as.character(i), ".rda") #<--------
  save(scenario, file = path)
  
  tab_trt = table(scenario[["sim_data"]][, treatment])
  rm_left_trt = names(sort(tab_trt))[1:rm_left]
  rm_right_trt = names(sort(tab_trt, decreasing = T))[1:rm_right]
  
  scenario1 = scenario
  
  scenario = rm.left(scenario1, rm_left_trt)
  path_left = paste0("../Scenarios/Scenario_8/seed", sep = as.character(i), ".rda") #<--------
  save(scenario, file = path_left)
  rm.gc(scenario)
  
  scenario = rm.right(scenario1, rm_right_trt)
  path_right = paste0("../Scenarios/Scenario_9/seed", sep = as.character(i), ".rda") #<--------
  save(scenario, file = path_right)
  rm.gc(scenario)
  
  rm.gc(scenario1)
}


# Note --------------------------------------------------------------------










