Template.matching = function(data, confounders, treatment, dim_proportion, nsample, matching){
  
  start_time = Sys.time()
  
  # Iperparametri (Dimensione template e numero di campioni) ----------------
  
  
  dim_template = min(table(data[, treatment]))*dim_proportion
  n = nsample
  
  
  # Vettore delle caratteristiche -------------------------------------------
  
  
  chain_tables = function(data, confounders){
    chain = c()
    for (c in confounders) {
      freq = table(data[,c])
      chain = c(chain, freq)
    }
    chain = chain/nrow(data)*100
    return(chain)
  }
  
  vec_template = chain_tables(data, confounders)
  
  
  # Vettore dei livelli -----------------------------------------------------
  
  #NOTA!!! In ASB ho tolto solo il complemento ad 1 delle dicotomiche, qui tolgo tutti gli ultimi livelli
  last_level = c()
  for (l in confounders) {
    nlev = length(levels(data[, l]))
    last_level = c(last_level, rep(1, (nlev-1)))
    last_level = c(last_level, 0)
  }
  last_level = last_level == 1
  
  # Estrazione di n campioni ---------------------------------------------
  
  
  set.seed(1234)
  
  campione = matrix(data = NA, nrow = dim_template, ncol = n)
  var_camp = matrix(data = NA, nrow = n, ncol = length(vec_template))
  for (i in 1:n) {
    campione[, i] = sample(1:nrow(data), dim_template, replace = F)
    camp = data[campione[, i], ]
    var_camp[i, ] = chain_tables(camp, confounders)
  }
  
  
  
  # Distanza di Mahalanobis da media popolazione generale -------------------
  
  
  dummy = matrix(data = NA, ncol = length(vec_template), nrow = nrow(data))
  j = 0
  for (c in confounders) {
    for (l in levels(data[, c])) {
      j = j + 1
      dummy[, j] = ifelse(data[, c] == l, 1, 0)
    }
  }
  dummy = as.data.frame(dummy)
  
  correl = cor(model.matrix(~.-1, data = dummy[, last_level]))
  distanze = mahalanobis(var_camp[,last_level], vec_template[last_level], correl)
  
  
  # Seleziono il template ---------------------------------------------------
  
  
  min = which(distanze == min(distanze))
  template = data[campione[, min], ]
  
  
  # Rendere i fattori numerici ----------------------------------------------
  
  
  as.double.factor = function(data) {
    tmp = data
    for (i in colnames(tmp)) {
      if (is.factor(tmp[, i])) {
        lev = levels(tmp[, i])
        for (g in lev) {
          tmp[, i] = ifelse(tmp[, i] == g, (match(g, lev) - 1), tmp[, i])
        }
      }
    }
    return(tmp)
  }
  
  
  # Match -------------------------------------------------------------------
  
  
  trt = unique(data[, treatment])
  template[, "identifier"] = 1
  matched = data[0,]
  
  
  for (t in trt) {
    tratt = data[data[, treatment] == t, ]
    tratt[, "identifier"] = 0
    step = rbind(tratt, template)
    
    if(matching == "Match"){
      num.step = as.double.factor(step) #devo trasformare tutto in numerico, ci sarà un modo migliore
      
      mtch = Matching::Match(Y = NULL, Tr = num.step[, "identifier"], X = num.step[, confounders], 
                             exact = TRUE, estimand = "ATT", M=1, ties=FALSE, distance.tolerance = 1)
      selected.index = mtch$index.control
      selected = step[selected.index, ]
    }
    
    if(matching == "cem"){
      cem_data = step[, c(confounders, "identifier")]
      sink("/dev/null")
      cem = cem::cem(treatment = "identifier", data = cem_data, baseline.group = "1")
      sink()
      
      strata = cem[["strata"]]
      tmp = cbind(step, strata)
      
      sel = function(s, data){
        pop = data[data[, "identifier"] == 0, ]
        s_pop = pop[pop[, "strata"] == s, ]
        selected = s_pop[sample(nrow(s_pop), 1), ]
        return(selected)
      }
      
      unique_strata = unique(strata)
      out = lapply(unique_strata, sel, tmp)
      
      selected = out[[1]]
      for (i in 2:length(out)) {
        selected = rbind(selected, out[[i]])
      }
      selected[, "strata"] = NULL
    }
    
    matched = rbind(matched, selected)
  }
  
  matched[, "identifier"] = NULL
  
  # Output ------------------------------------------------------------------
  
  output.maker = function(balanced_data, ASB_pre, ASB_post, time){
    output = list()
    output$balanced_data = balanced_data
    output$ ASB_pre = ASB_pre
    output$ASB_post = ASB_post
    output$Time = time
    return(output)
  }
  
  ASB_pre = MARMoT::ASB(data, confounders, treatment)
  ASB_post = MARMoT::ASB(matched, confounders, treatment)
  
  time = difftime(Sys.time(), start_time, units = "mins")

  output = output.maker(matched, ASB_pre, ASB_post, time)
  
  return(output)
}