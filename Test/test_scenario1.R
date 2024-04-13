
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

origin = "../Scenarios/Scenario_1/"                         #<----------
path = paste0(origin, "seed")
replications = 10
save_path = paste0(origin, "out")                                   


# MARMoT ------------------------------------------------------------------
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



















