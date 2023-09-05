# Purpose: To profile across a range of values for a given parameter
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date 9/4/23

like_prof = function(em_inputs, sim, assessment_name, share_M = FALSE) {
  
  # mle_vals = mle_vals # Define maximum likelihood estimate value
  parameter_names = c("ln_fish_selpars", "ln_srv_selpars", "ln_M", "RecPars")
  # OM values
  om_vals = list(log(c(fish_len_slope, fish_len_midpoint)), log(c(srv_len_slope, srv_len_midpoint)), 
                 log(M), log(r0))

  # set up inputs
  data = em_inputs$data # data
  like_prof_all = data.frame()
  
  for(p in 1:length(parameter_names)) {
    
    # get parameter length
    par_len = length(parameters[names(parameters) == parameter_names[p]][[1]])
    if(parameter_names[p] == "RecPars") par_len = 1 # only profile across r0
    if(share_M == TRUE) par_len = 1 # only profile across a shared value of M
    
    for(l in 1:par_len) {
      map = em_inputs$map # initialize mapping
      seq_values = 1:par_len
      # Create a vector with NA at position p and sequential numbers elsewhere
      if(l == 1 & parameter_names[p] != "RecPars" & share_M == FALSE) {
        map_vals = ifelse(seq_values == l, NA, seq_values-1)
        # set up mapping
        map = c(map, list(factor(map_vals))) # combine factors 
        names(map)[length(map)] = parameter_names[p] # name factor
      }
      if(l > 1 & share_M == FALSE) {
        map_vals = ifelse(seq_values == l, NA, seq_values)
        map_vals = ifelse(map_vals %in% c(NA, 1), map_vals, seq_values - 1)
        # set up mapping
        map = c(map, list(factor(map_vals))) # combine factors 
        names(map)[length(map)] = parameter_names[p] # name factor
      } # end if statement
      
      # get profile values
      prof_vals = seq(om_vals[[p]][l] - 1.5, om_vals[[p]][l] + 1.5, 0.1)
      
      if(parameter_names[p] == "RecPars") map$RecPars = factor(c(NA, NA))
      if(share_M == TRUE) { # if we are sharing the M
        map$ln_M = factor(c(NA, NA)) 
        prof_vals = seq(mean(om_vals[[p]]) - 1.5, mean(om_vals[[p]]) + 1.5, 0.1)
      } # end if statement for sharing m
      
      # setup data frame for storage
      like_prof_tmp = data.frame(jnLL = NA, Catch = NA, Fish_Index = NA, Fish_Age = NA, Fish_Len = NA,
                             Srv_Index = NA, Srv_Age = NA, Srv_Len = NA, rec_pen = NA, prof_val = NA, 
                             om_val = NA, sim = NA)
      
      for(v in 1:length(prof_vals)) {
        parameters = em_inputs$parameters # reset parameters
        # which value to fix parameter at
        parameters[names(parameters) == parameter_names[p]][[1]][l] = prof_vals[v]
        # different way of profiling for single M (repeating 2 to make sure the same value gets mapped to both sexes)
        if(share_M == TRUE) parameters[names(parameters) == parameter_names[p]][[1]] = rep(prof_vals[v], 2)
        tryCatch({ # tryCatch keeps the loop going despite erros 
          
          # make ad obj
          model_fxn = TMB::MakeADFun(data, parameters, map, random = NULL,
                                     DLL= "Sex_Str_EM", silent = TRUE,
                                     checkParameterOrder = TRUE, tracepar = TRUE)
          
          # Optimize model here w/ nlminb
          mle_optim <- stats::nlminb(model_fxn$par, model_fxn$fn, model_fxn$gr,
                                     control = list(iter.max = 1e5, eval.max = 1e5))
          Report <- model_fxn$report(model_fxn$env$last.par.best) # get regular report
          sd_rep <- TMB::sdreport(model_fxn) # get sdreport
          
          like_prof_tmp[v,] = c(Report$jnLL, sum(Report$catch_nLL), 
                                                    sum(Report$fish_index_nLL), sum(Report$fish_age_comp_nLL),
                                                    sum(Report$fish_len_comp_nLL), sum(Report$srv_index_nLL), sum(Report$srv_age_comp_nLL),
                                                    sum(Report$srv_len_comp_nLL), Report$rec_nLL, prof_vals[v], om_vals[[p]][l],
                                                    sim)

        }, error = function(error) {cat("ERROR :",conditionMessage(error), "\n")}) # end try catch statement
      
        } # end v loop
      
      like_prof_tmp$par_name = paste(parameter_names[p], l, sep = "_") # rename parameter name
      like_prof_tmp$assessment_name = assessment_name # differentiate assessments
      like_prof_all = rbind(like_prof_tmp, like_prof_all) # bind values
      print(l)
    } # end l loop for parameter ordering
  } # end p loop for parameter names
  
  return(like_prof_all)
} # end function