

# Run ---------------------------------------------------------------------


run = function(path, replications, exp, confounders, 
               treatment, reduced_confounders, outcome, id){
  
  out = multiple.seed.run(path, replications, exp, confounders, 
                          treatment, reduced_confounders, outcome, id)
  mean_out = mean.out(out)
  grap1 = graph.confounders(mean_out$og_matrix, mean_out$bal_matrix)
  grap2 = graph.scatter(out)
  
  stats_results = list(mean_out$trt_stats, mean_out$og_stats, mean_out$bal_stats, 
                       mean_out$OR_stats)
  names(stats_results) = c("Trt stats", "ASB Pre", "ASB Post", "OR stats")
  
  graphs_results = list(grap1, grap2)
  names(graphs_results) = c("Pre-Post", "Scatterplot")
  
  call = as.character(exp)
  names(call) = "call"
  call_list = list(call)
  names(call_list) = "call"
  
  final_results = c(call_list, stats_results, graphs_results)
  
  return(final_results)
}


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------


# Function necessities ----------------------------------------------------


if(system.file(package="ggplot2") == "") install.packages("ggplot2")
library(ggplot2)
if(system.file(package="ggrepel") == "") install.packages("ggrepel")
library(ggrepel)
if(system.file(package="epitools") == "") install.packages("epitools")
library(epitools)


# To each single ----------------------------------------------------------


# Diagnostic stats --------------------------------------------------------


diagnostic.stat = function(out, data, id, treatment, confounders){
  
  diag_names = c("Min", "1st quartile", "Median", "3rd quartile", "Max", "Mean", 
                 "Over 5%", "Over 10%", "% of inclusion", "Rep. mean", "Rep. max", 
                 "Trt. size", "Pop. size", "% profiles", "Time (m)")
  diagnostic = c(rep(NaN, length(diag_names)))
  names(diagnostic) = diag_names
  
  ASB_stats = ASB(out$balanced_data, confounders, treatment)
  diagnostic[1:length(out$ASB_post$Stat)] = out$ASB_post$Stat
  
  diagnostic["Time (m)"] = out$Time
  
  diagnostic["Pop. size"] = nrow(out$balanced_data)
  
  id_table = sort(table(out$balanced_data[, id]), decreasing = T)
  
  diagnostic["Rep. mean"] = mean(id_table)
  
  diagnostic["Rep. max"] = id_table[1]
  
  id_pre = length(unique(data[, id]))
  id_post = length(unique(out$balanced_data[, id]))
  
  diagnostic["% of inclusion"] = (1 - (id_pre - id_post) / id_pre) * 100
  
  diagnostic["Trt. size"] = mean(table(out$balanced_data[, treatment]))
  
  diagnostic["% profiles"] = perc.profiles(data[, confounders], 
                                           out$balanced_data[, confounders])
  
  return(diagnostic)
}


# Diagnostic pop distribution  --------------------------------------------


diagnostic.pop.dist = function(data, treatment, confounders){
  
  tab = table(data[, treatment])
  stat = rep(NA, 9)
  names(stat) = c("Num", "Min", "1st quartile", "Median", "3rd quartile", "Max", "Mean", "Sd", "Profiles")
  stat["Num"] = nrow(data)
  stat["Min"] = min(tab)
  stat["1st quartile"] = quantile(tab, 0.25)
  stat["Median"] = median(tab)
  stat["3rd quartile"] = quantile(tab, 0.75)
  stat["Max"] = max(tab)
  stat["Mean"] = mean(tab)
  stat["Sd"] = stats::sd(tab)
  stat["Profiles"] = n.profiles(data[, confounders])
  
  return(stat)
}


# OR ----------------------------------------------------------------------


diagnostic.OR = function(beta, bal_data, treatment, outcome){
  
  trt_beta_names = names(beta)[grepl("^PRVDR_NUMTRT_", names(beta))]
  or_true = exp(beta[trt_beta_names])
  if(length(or_true) == length(levels(bal_data[, treatment]))){
     or_true = or_true / or_true[1]
     or_true = or_true[-1]
  }
  or_true = c(1, or_true)
  names(or_true) = NULL
  
  compute.or = function(data, treatment, outcome){
    tab = table(data[, treatment], data[, outcome])
    or = oddsratio.fisher(tab)$measure
    or = or[, "estimate"]
    return(or)
  }
  or_bal = compute.or(bal_data, treatment, outcome)
  or_df = rbind(or_true, or_bal)
  rownames(or_df) = c("True_OR", "Balanced_OR")
  or_df = t(or_df)
  or_df = data.frame(or_df)
  
  Bias = abs(or_df[, "True_OR"] - or_df[, "Balanced_OR"])/or_df[, "True_OR"]
  or_df = cbind(or_df, Bias)
  
  bias_mean = mean(Bias)
  bias_sd = stats::sd(Bias)
  stat = c(bias_mean, bias_sd)
  stat = stat
  names(stat) = c("Relative bias mean", "Relative bias sd")
  
  return(stat)
}


# Profiling ---------------------------------------------------------------


profiling.single = function(row){
  out = paste0(row, collapse = "_")
  return(out)
}

profiling = function(dati){
  out = apply(dati, 1, profiling.single)
  return(out)
}

n.profiles = function(data){
  prof = profiling(data)
  tab = table(prof)
  n = length(tab)
  return(n)
}

perc.profiles = function(data, bal_data){
  og = n.profiles(data)
  bal = n.profiles(bal_data)
  out = bal / og * 100
  return(out)
}

# Single seed run diagnostic ----------------------------------------------


single.seed.run.diagnostic = function(data, treatment, confounders, beta, bal_output, outcome, id){
  
  trt_stats = diagnostic.pop.dist(data, treatment, confounders)
  OR_stats = diagnostic.OR(beta, bal_output$balanced_data, treatment, outcome)
  bal_stats = diagnostic.stat(bal_output, data, id, treatment, confounders)
  #sink("/dev/null")
  bal_matrix = ASB(bal_output$balanced_data, confounders, treatment, verbose = F)$Matrix
  og_ASB = ASB(data, confounders, treatment, verbose = F)
  #sink()
  og_stats = og_ASB$Stat
  og_matrix = og_ASB$Matrix
  og_trt_table = table(data[, treatment])
  
  out = list(trt_stats, OR_stats, bal_stats, bal_matrix, og_stats, og_matrix, og_trt_table)
  names(out) = c("trt_stats", "OR_stats", "bal_stats", "bal_matrix", "og_stats", "og_matrix", "og_trt_table")
  return(out)
}


# Single seed run ---------------------------------------------------------


single.seed.run = function(exp, scenario, confounders, treatment, reduced_confounders, outcome, id){
  
  bal_output = eval(exp)
  out = single.seed.run.diagnostic(scenario$sim_data, treatment, confounders, 
                                   scenario$true_effects, bal_output, outcome, id)
  
  return(out)
}


# Multiple seed runs ------------------------------------------------------


multiple.seed.run = function(path, replications, exp, confounders, treatment, reduced_confounders, outcome, id){
  out = list()
  for (i in 1:replications) {
    print(paste0("Replication ", i))
    complete_path = paste0(path, sep = as.character(i), ".rda")
    load(complete_path)
    #sink("/dev/null")
    tmp = single.seed.run(exp, scenario, confounders, treatment, reduced_confounders, outcome, id)
    #sink()
    out[[i]] = tmp
  }
  return(out)
}


# Mean of the results -----------------------------------------------------


div = function(x, n){
  out = x / n
  return(out)
}

mean.out = function(out){
  sum_list = mapply("+", out[[1]], out[[2]], SIMPLIFY = FALSE)
  n = length(out)
  for (h in 3:n) {
    sum_list = mapply("+", sum_list, out[[h]], SIMPLIFY = FALSE)
  }
  out = lapply(sum_list, div, n)
  return(out)
}


# On the mean -------------------------------------------------------------


# Diagnostic graphs alternative -------------------------------------------


graph.confounders = function(og_matrix, bal_matrix){
  
  pre = og_matrix
  post = bal_matrix
  
  tmp_bef = data.frame(pre[, 0])
  tmp_aft = data.frame(pre[, 0])
  
  tmp_bef$value = rowMeans(pre)
  tmp_bef$group = rownames(tmp_bef)
  tmp_bef$time = "Before"
  
  tmp_aft$value = rowMeans(post)
  tmp_aft$group = rownames(tmp_aft)
  tmp_aft$time = "After"
  
  tmp = rbind(tmp_bef, tmp_aft) 
  tmp$time = factor(tmp$time, levels = c("Before", "After"), ordered = T)
  
  bef_aft = ggplot(tmp, aes(x = time, y = value, 
                            group = group, label = ifelse(time == "Before", group, ""))) + 
    geom_point(aes(color = group)) + geom_line(aes(color = group)) + 
    labs(x = NULL, y = "Mean ASB") + 
    geom_label_repel(size = 3, max.overlaps = Inf) + guides(color = "none") + theme_minimal()
  
  return(bef_aft)
}


# Other -------------------------------------------------------------------


cbind.out = function(out){
  
  cbind_list = mapply(cbind, out[[1]], out[[2]], SIMPLIFY = FALSE)
  n = length(out)
  for (h in 3:n) {
    cbind_list = mapply(cbind, cbind_list, out[[h]], SIMPLIFY = FALSE)
  }
  
  return(cbind_list)
}


# Scatterplot  ------------------------------------------------------------


graph.scatter = function(out){
  
  cbind_out = cbind.out(out)
  x = c(cbind_out$og_trt_table)
  y = colMeans(cbind_out$bal_matrix)
  df = data.frame(cbind(x, y))
  
  bal_dim = ggplot(df, aes(x = x, y = y)) + geom_point() + 
    geom_smooth(method = "loess", se = FALSE, span = 0.75, color = "#A2C8FF", formula = y ~ x) + 
    labs(x = "Treatment population before balancing", y = "Mean ASB") + 
    theme_minimal()
  
  return(bal_dim)
}
