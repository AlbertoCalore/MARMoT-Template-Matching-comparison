

gen.data = function(dati, n_sim, n_treat, noise, seed){
  
  set.seed(seed)
  
  gen.data.no.seed = function(dati, n_sim, n_treat, noise, seed){
    
    
    # Recurrent functions -----------------------------------------------------
    
    
    noise.maker = function(v, noise){
      noise = noise * sd(v)
      nor = rnorm(length(v), 0, noise)
      v = v + nor
      v = round(v, digits = 2)
      v[v > 1] = 1
      v[v < 0] = 0
      return(v)
    }
    
    noise.maker.multi = function(m, noise){
      for (h in 1:ncol(m)) {
        m[, h] = noise.maker(m[, h], noise)
      }
      m = m / rowSums(m)
      return(m)
    }
    
    draw.bin.single = function(prob_single){
      out = rbinom(1, 1, prob_single)
      return(out)
    }
    
    draw.bin = function(prob){
      out = sapply(prob, draw.bin.single)
      return(out)
    }
    
    logistic.sim = function(dati, sim_data, formula_sim, noise){
      glm_bin = glm(formula_sim, data = dati, family = binomial)
      glm_bin_pred = predict(glm_bin, newdata = sim_data, type = "response")
      glm_bin_pred = noise.maker(glm_bin_pred, noise)
      pred = data.frame(glm_bin_pred)
      pred$sim = NA
      for (h in unique(pred$glm_bin_pred)) {
        n = nrow(pred[pred["glm_bin_pred"] == h, ])
        pred[pred["glm_bin_pred"] == h, "sim"] = rbinom(n, 1, h)
      }
      out = factor(pred[, "sim"])
      return(out)
    }
    
    draw.multi.single = function(prob_row){
      out_vec = rmultinom(1, 1, prob_row)
      out = which(out_vec == 1)
      names(out) = c()
      return(out)
    }
    
    draw.multi = function(prob){
      out = apply(prob, 1, draw.multi.single)
      return(out)
    }
    
    multilogit.sim = function(dati, sim_data, formula_sim, noise){
      words = unlist(strsplit(as.character(formula_sim), split = "~"))
      Y = dati[, words[2]]
      formula_sim = as.formula(paste0("Y", sep = " ~ ", words[3]))
      glm_multi = multinom(formula_sim, dati, maxit = 500)
      pred = predict(glm_multi, newdata = sim_data, type = "probs")
      pred = noise.maker.multi(pred, noise)
      resp = draw.multi(pred)
      return(resp)
    }
    
    
    # Data cleaning -----------------------------------------------------------
    
    
    treatment = "PRVDR_NUM"
    
    dati$NEW_AGE_CAT5 = droplevels(dati$NEW_AGE_CAT5)
    dati$NEW_ADMISSION_DIAGNOSIS = NULL
    dati$NEW_SUM_COMORBIDITIES = NULL
    dati$NEW_COMORBIDITIES = NULL
    dati$DESYNPUF_ID = NULL
    dati$NEW_AGE_CAT10 = NULL
    dati$NEW_LENGTH_OF_STAY = NULL
    dati$NEW_COST_X_DAY = NULL
    
    to_factor = c("SP_CHRNKIDN", "SP_COPD", "SP_DIABETES", "NEW_HYPERTENSION", "NEW_PERIPHVASC", "NEW_90DAYSDEATH", "NEW_30DAYREADMISSION")
    for (i in to_factor) {
      dati[, i] = factor(dati[, i])
    }
    
    dati$NEW_90DAYSDEATH = NULL
    
    sim_data = data.frame(matrix(NA, ncol = ncol(dati), nrow = n_sim))
    colnames(sim_data) = colnames(dati)
    
    
    # Some lists --------------------------------------------------------------
    
    
    outcomes = c("NEW_90DAYSDEATH", "NEW_LENGTH_OF_STAY", "NEW_30DAYREADMISSION", "NEW_COST_X_DAY")
    confounders = c(colnames(dati)[!colnames(dati) %in% outcomes])
    confounders = confounders[confounders != treatment]
    
    quantitatives = c()
    qualitatives = c()
    for (i in colnames(dati)) {
      if(is.numeric(dati[, i]))
      {quantitatives = c(quantitatives, i)} else{
        qualitatives = c(qualitatives, i)
      }
    }
    
    table(colnames(dati) %in% c(quantitatives, qualitatives))
    
    outcomes_qual = outcomes[outcomes %in% qualitatives]
    outcomes_quant = outcomes[outcomes %in% quantitatives]
    
    confounders_qual = confounders[confounders %in% qualitatives]
    confounders_quant = confounders[confounders %in% quantitatives] 
    
    
    # Basic characteristics ---------------------------------------------------
    
    
    #sex
    prob_sex1 = prop.table(table(dati[, "BENE_SEX_IDENT_CD"]))[2]
    sim_data[, "BENE_SEX_IDENT_CD"] = rbinom(n_sim, 1, prob_sex1)
    sim_data[, "BENE_SEX_IDENT_CD"] = factor(sim_data[, "BENE_SEX_IDENT_CD"])
    
    #race
    formula_sim = as.formula("NEW_RACE ~ BENE_SEX_IDENT_CD")
    sim_data[, "NEW_RACE"] = logistic.sim(dati, sim_data, formula_sim, noise)
    
    #age
    formula_sim = as.formula("NEW_AGE_CAT5 ~ BENE_SEX_IDENT_CD + NEW_RACE")
    sim_data[, "NEW_AGE_CAT5"] = multilogit.sim(dati, sim_data, formula_sim, noise)
    sim_data[, "NEW_AGE_CAT5"] = factor(sim_data[, "NEW_AGE_CAT5"], 
                                        levels = sort(unique(sim_data[, "NEW_AGE_CAT5"])), ordered = T)
    levels(sim_data[, "NEW_AGE_CAT5"]) = levels(dati[, "NEW_AGE_CAT5"])
    
    
    # Diseases ----------------------------------------------------------------
    
    
    # diabetes
    formula_sim = as.formula("SP_DIABETES ~ BENE_SEX_IDENT_CD + NEW_AGE_CAT5 + NEW_RACE")
    sim_data[, "SP_DIABETES"] = logistic.sim(dati, sim_data, formula_sim, noise)
    
    #periperal vascular disease
    formula_sim = as.formula("NEW_PERIPHVASC ~ BENE_SEX_IDENT_CD + NEW_AGE_CAT5 + NEW_RACE")
    sim_data[, "NEW_PERIPHVASC"] = logistic.sim(dati, sim_data, formula_sim, noise)
    
    #Hypertension
    formula_sim = as.formula("NEW_HYPERTENSION ~ BENE_SEX_IDENT_CD + NEW_AGE_CAT5 + NEW_RACE")
    sim_data[, "NEW_HYPERTENSION"] = logistic.sim(dati, sim_data, formula_sim, noise)
    
    #COPD
    formula_sim = as.formula("SP_COPD ~ BENE_SEX_IDENT_CD + NEW_AGE_CAT5 + NEW_RACE + SP_DIABETES")
    sim_data[, "SP_COPD"] = logistic.sim(dati, sim_data, formula_sim, noise)
    
    #Cronic kidney disease
    formula_sim = as.formula("SP_CHRNKIDN ~ BENE_SEX_IDENT_CD + NEW_AGE_CAT5 + NEW_RACE + SP_DIABETES + SP_COPD + NEW_HYPERTENSION")
    sim_data[, "SP_CHRNKIDN"] = logistic.sim(dati, sim_data, formula_sim, noise)
    
    
    # Procedure ---------------------------------------------------------------
    
    
    formula_sim = as.formula("NEW_PROCED_PCI ~ BENE_SEX_IDENT_CD + NEW_AGE_CAT5 + SP_DIABETES + 
                         SP_COPD + NEW_HYPERTENSION + SP_CHRNKIDN + NEW_PERIPHVASC")
    sim_data[, "NEW_PROCED_PCI"] = logistic.sim(dati, sim_data, formula_sim, noise)
    
    
    # Check -------------------------------------------------------------------
    
    #print("Sanity check")
    #par(mfrow = c(1,2))
    #for(i in confounders_qual){
    #  print(i)
    #  plot(prop.table(table(dati[, i])))
    #  print(prop.table(table(dati[, i])))
    #  plot(prop.table(table(sim_data[, i])))
    #  print(prop.table(table(sim_data[, i])))
    #}
    #par(mfrow = c(1,1))
    
    
    # Treatments --------------------------------------------------------------

    
    Y = class.ind(dati[, "PRVDR_NUM"])
    formula_sim = as.formula("Y ~ BENE_SEX_IDENT_CD + SP_CHRNKIDN + SP_COPD + 
                                SP_DIABETES + NEW_PROCED_PCI + NEW_HYPERTENSION + 
                                NEW_PERIPHVASC + NEW_RACE + NEW_AGE_CAT5")
    trt_model = multinom(formula_sim, data = dati, maxit = 500)
    
    beta = coef(trt_model)
    x = model.matrix(as.formula(" ~ BENE_SEX_IDENT_CD + SP_CHRNKIDN + SP_COPD + 
                               SP_DIABETES + NEW_PROCED_PCI + NEW_HYPERTENSION + 
                               NEW_PERIPHVASC + NEW_RACE + NEW_AGE_CAT5"), sim_data)
    
    #new treatments
    n_fake_trt = n_treat - length(unique(dati[, "PRVDR_NUM"]))
    fake_beta = matrix(0, nrow = n_fake_trt, ncol = ncol(beta)) 
    colnames(fake_beta) = colnames(beta)
    
    for (g in colnames(fake_beta)) {
      mn = quantile(beta[, g], 0.15)
      mx = quantile(beta[, g], 0.75)
      fake_col = runif(nrow(fake_beta), min = mn, max = mx)
      fake_beta[, g] = fake_col
    }
    
    new_beta = rbind(beta, fake_beta)
    rownames(new_beta) = paste0("TRT_0", c(2:(nrow(new_beta)+1)))
    pred = x %*% t(new_beta)
    rm(x)
    gc()
    resp = exp(pred) / (1 + rowSums(exp(pred)))
    TRT_01 = 1 - rowSums(resp)
    resp = cbind(TRT_01, resp)
    
    sim_data[, treatment] = draw.multi(resp)
    sim_data[, treatment] = factor(sim_data[, treatment], levels = c(1:n_treat))
    levels(sim_data[, treatment]) = colnames(resp)
    
    #sim_data[, treatment] = trt_sim
    colnames(sim_data)[colnames(sim_data) == "trt_sim"] = treatment
    
    
    # Outcome -----------------------------------------------------------------
    

    formula_out = as.formula("NEW_30DAYREADMISSION ~ BENE_SEX_IDENT_CD + 
                         SP_CHRNKIDN + SP_COPD + SP_DIABETES + NEW_PROCED_PCI +
                         NEW_HYPERTENSION + NEW_PERIPHVASC + NEW_RACE + 
                         NEW_AGE_CAT5 + PRVDR_NUM")
    glm_bin = glm(formula_out, data = dati, family = binomial)
    b = coef(glm_bin)
    
    
    x = model.matrix(as.formula(" ~ BENE_SEX_IDENT_CD + 
                         SP_CHRNKIDN + SP_COPD + SP_DIABETES + NEW_PROCED_PCI +
                         NEW_HYPERTENSION + NEW_PERIPHVASC + NEW_RACE + 
                         NEW_AGE_CAT5 + PRVDR_NUM"), sim_data)
    
    b_trt = c(b[grepl(paste0("^", treatment), names(b))])
    fake_b_trt = runif(n_fake_trt, quantile(b_trt, 0.1), quantile(b_trt, 0.9))
    b = c(b, fake_b_trt)
    names(b)[(length(b) - n_treat + 2):length(b)] = paste0(treatment, paste0("TRT", sep = "_0", c(2:n_treat)))
    
    prob = plogis(x %*% b)
    rm(x)
    gc()
    prob = noise.maker(prob, noise)
    y = draw.bin(prob)
    sim_data[, "NEW_30DAYREADMISSION"] = factor(y)
    
    
    # ID ----------------------------------------------------------------------
    
    
    sim_data$ID = paste0("ID", c(1:n_sim))
    
    
    # Output ------------------------------------------------------------------
    
    
    scenario = c(n_sim, n_treat, noise, seed)
    names(scenario) = c("n_sim", "n_treat", "noise", "seed")
    
    #out = list(sim_data, b, new_beta, scenario)
    #names(out) = c("sim_data", "true_effects", "treatment_beta", "scenario")
    
    out = list(sim_data, b, scenario)
    names(out) = c("sim_data", "true_effects", "scenario")
    return(out)
    
  }
  
  out = gen.data.no.seed(dati, n_sim, n_treat, noise, seed)
  return(out)
}












