
# Setup -------------------------------------------------------------------


rm(list=ls())
setwd("~/Server/Test/")


# Libraries ---------------------------------------------------------------


if(system.file(package="ggplot2") == "") install.packages("ggplot2")
library(ggplot2)
if(system.file(package="Matching") == "") install.packages("Matching")
if(system.file(package="MARMoT") == "") install.packages("MARMoT")
library(MARMoT)
source("../Functions/Template_matching.R")
source("../Functions/all_for_a_run.R")


# Common parameters -------------------------------------------------------


treatment = "PRVDR_NUM"
outcome = "NEW_30DAYREADMISSION"
id = "ID"
confounders = c("BENE_SEX_IDENT_CD", "SP_CHRNKIDN", "SP_COPD", "SP_DIABETES", 
                "NEW_PROCED_PCI", "NEW_HYPERTENSION", "NEW_PERIPHVASC", 
                "NEW_RACE", "NEW_AGE_CAT5")
reduced_confounders = confounders[!confounders %in% c("NEW_PERIPHVASC")]

origin = "../Scenarios/Scenario_6/"                         #<----------
path = paste0(origin, "seed")
replications = 10
save_path = paste0(origin, "out")                                   


# MARMoT ------------------------------------------------------------------
# -------------------------------------------------------------------------


exp1 = expression(MARMoT(scenario$sim_data, confounders, treatment, 
                         n.cores = 6, reference = "median", caliper = 0.25))
set.seed(1234)
out1 = run(path, replications, exp1, confounders, treatment, reduced_confounders, outcome, id)
save(out1, file = paste0(save_path, sep = 1, ".rda"))


# -------------------------------------------------------------------------


exp2 = expression(MARMoT(scenario$sim_data, confounders, treatment, 
                         n.cores = 6, reference = "mean0", caliper = 0.25))
set.seed(1234)
out2 = run(path, replications, exp2, confounders, treatment, reduced_confounders, outcome, id)
save(out2, file = paste0(save_path, sep = 2, ".rda"))


# -------------------------------------------------------------------------


exp3 = expression(MARMoT(scenario$sim_data, confounders, treatment, 
                         n.cores = 6, reference = "max", caliper = 0.25))
set.seed(1234)
out3 = run(path, replications, exp3, confounders, treatment, reduced_confounders, outcome, id)
save(out3, file = paste0(save_path, sep = 3, ".rda"))


# -------------------------------------------------------------------------


exp4 = expression(MARMoT(scenario$sim_data, confounders, treatment, 
                         n.cores = 6, reference = "median", caliper = 0.1))
set.seed(1234)
out4 = run(path, replications, exp4, confounders, treatment, reduced_confounders, outcome, id)
save(out4, file = paste0(save_path, sep = 4, ".rda"))


# -------------------------------------------------------------------------


exp5 = expression(MARMoT(scenario$sim_data, confounders, treatment, 
                         n.cores = 6, reference = "mean0", caliper = 0.1))
set.seed(1234)
out5 = run(path, replications, exp5, confounders, treatment, reduced_confounders, outcome, id)
save(out5, file = paste0(save_path, sep = 5, ".rda"))


# -------------------------------------------------------------------------


exp6 = expression(MARMoT(scenario$sim_data, confounders, treatment, 
                         n.cores = 6, reference = "max", caliper = 0.1))
set.seed(1234)
out6 = run(path, replications, exp6, confounders, treatment, reduced_confounders, outcome, id)
save(out6, file = paste0(save_path, sep = 6, ".rda"))


# -------------------------------------------------------------------------


exp7 = expression(MARMoT(scenario$sim_data, reduced_confounders, treatment, 
                         n.cores = 6, reference = "median", caliper = 0.25))
set.seed(1234)
out7 = run(path, replications, exp7, confounders, treatment, reduced_confounders, outcome, id)
save(out7, file = paste0(save_path, sep = 7, ".rda"))


# -------------------------------------------------------------------------


exp8 = expression(MARMoT(scenario$sim_data, reduced_confounders, treatment, 
                         n.cores = 6, reference = "mean0", caliper = 0.25))
set.seed(1234)
out8 = run(path, replications, exp8, confounders, treatment, reduced_confounders, outcome, id)
save(out8, file = paste0(save_path, sep = 8, ".rda"))


# -------------------------------------------------------------------------


exp9 = expression(MARMoT(scenario$sim_data, reduced_confounders, treatment, 
                         n.cores = 6, reference = "max", caliper = 0.25))
set.seed(1234)
out9 = run(path, replications, exp9, confounders, treatment, reduced_confounders, outcome, id)
save(out9, file = paste0(save_path, sep = 9, ".rda"))


# -------------------------------------------------------------------------


exp10 = expression(MARMoT(scenario$sim_data, reduced_confounders, treatment, 
                         n.cores = 6, reference = "median", caliper = 0.1))
set.seed(1234)
out10 = run(path, replications, exp10, confounders, treatment, reduced_confounders, outcome, id)
save(out10, file = paste0(save_path, sep = 10, ".rda"))


# -------------------------------------------------------------------------


exp11 = expression(MARMoT(scenario$sim_data, reduced_confounders, treatment, 
                         n.cores = 6, reference = "mean0", caliper = 0.1))
set.seed(1234)
out11 = run(path, replications, exp11, confounders, treatment, reduced_confounders, outcome, id)
save(out11, file = paste0(save_path, sep = 11, ".rda"))


# -------------------------------------------------------------------------


exp12 = expression(MARMoT(scenario$sim_data, reduced_confounders, treatment, 
                         n.cores = 6, reference = "max", caliper = 0.1))
set.seed(1234)
out12 = run(path, replications, exp12, confounders, treatment, reduced_confounders, outcome, id)
save(out12, file = paste0(save_path, sep = 12, ".rda"))


# -------------------------------------------------------------------------


# -------------------------------------------------------------------------


# Template matching -------------------------------------------------------
# -------------------------------------------------------------------------


exp20 = expression(Template.matching(scenario$sim_data, confounders, treatment, 3/5 , 10000, "cem"))
set.seed(1234)
out20 = run(path, replications, exp20, confounders, treatment, reduced_confounders, outcome, id)
save(out20, file = paste0(save_path, sep = 20, ".rda"))


# -------------------------------------------------------------------------


exp21 = expression(Template.matching(scenario$sim_data, confounders, treatment, 3/5 , 10000, "Match"))
set.seed(1234)
out21 = run(path, replications, exp21, confounders, treatment, reduced_confounders, outcome, id)
save(out21, file = paste0(save_path, sep = 21, ".rda"))


# -------------------------------------------------------------------------


exp22 = expression(Template.matching(scenario$sim_data, reduced_confounders, treatment, 3/5 , 10000, "cem"))
set.seed(1234)
out22 = run(path, replications, exp22, confounders, treatment, reduced_confounders, outcome, id)
save(out22, file = paste0(save_path, sep = 22, ".rda"))


# -------------------------------------------------------------------------


exp23 = expression(Template.matching(scenario$sim_data, reduced_confounders, treatment, 3/5 , 10000, "Match"))
set.seed(1234)
out23 = run(path, replications, exp23, confounders, treatment, reduced_confounders, outcome, id)
save(out23, file = paste0(save_path, sep = 23, ".rda"))


# -------------------------------------------------------------------------



















