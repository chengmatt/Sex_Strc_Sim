# Purpose: Utility functions 
# Author: Matthew LH. Cheng (UAF-CFOS)
# Date: 8/3/23

#' Title Get Age-Length Transition Matrix for use in Age-Structured Models
#' 
#' @param age_bins vector of ages
#' @param len_bins vector of lengths
#' @param mean_length mean length at age 
#' @param cv cv for LAA
#'
#' @return
#' @export
#'
#' @examples
get_al_trans_matrix = function(age_bins, len_bins, mean_length, sd) {
  
  # Construct age length matrix
  age_length = matrix(0.0, nrow = length(age_bins), ncol = length(len_bins))

  for(a in 1:length(age_bins)) {
    for(l in 1:length(len_bins)) {
      
      if (l == 1) { # Probability of being between 1st and 2nd length bin given age a
        age_length[a, l] = pnorm(len_bins[2], mean_length[a], sd[a])
      } else if (l == length(len_bins)) { # Probability of being larger than the last length bin given age a
        age_length[a, l] = 1 - pnorm(len_bins[length(len_bins)], mean_length[a], sd[a])
      } else { # a of being in between length bins given age a
        age_length[a, l] = pnorm(len_bins[l+1], mean_length[a], sd[a]) -  
                               pnorm(len_bins[l], mean_length[a], sd[a])
      }
    } # end l loop
  } # end a loop
  return(age_length)
} # end function

#' Title Generalized Logistic Function
#'
#' @param slope Slope of logistic function
#' @param bins Number of bins
#' @param midpoint Midpoint of logistic function
#'
#' @return
#' @export
#'
#' @examples
logist = function(slope, bins, midpoint) {
  return(  1 / (1 + exp(-slope * (bins - midpoint)) ))
} # end function

#' Title Logistic Function with a95 parameterization
#'
#' @param a95 age at 95% selex
#' @param bins Number of bins
#' @param a50 Midpoint of logistic function
#'
#' @return
#' @export
#'
#' @examples

logist_19 = function(a50, bins, a95) {
  return(1 / (1 + exp(-log(19) * ((bins - a50) / (a95 - a50)) )))
} # end function

#' Title Von Bertalannfy Growth Function
#'
#' @param age_bins vector of ages
#' @param k growth rate parameter
#' @param L_inf Asymptotic size
#' @param t0 Size at age0
#' @param cv coefficient of variation for generating LAA samples
#'
#' @return
#' @export
#'
#' @examples
vonB_est = function(age_bins, k , L_inf, t0, sd) {
  return(L_inf * (1-exp(-k*(age_bins -t0)) ) + rnorm(1, 0, sd))
} #end function


#' Title Take additional newton steps with TMB model
#'
#' @param n.newton number of additional newton steps we want to take
#' @param ad_model MakeADFUN model object
#' @param mle_optim Optimized model object
#'
#' @return
#' @export
#'
#' @examples
add_newton = function(n.newton, ad_model, mle_optim) {
  
  tryCatch(expr = for(i in 1:n.newton) {
    g = as.numeric(ad_model$gr(mle_optim$par))
    h = optimHess(mle_optim$par, fn = ad_model$fn, gr = ad_model$gr)
    mle_optim$par = mle_optim$par - solve(h,g)
    mle_optim$objective = ad_model$fn(mle_optim$par)
  }, error = function(e){e})
  
}



#' Title Get biological information for sex-aggregated and sex-specific WAA, LAA, and age-length transition matrix
#'
#' @param n_sexes Number of sexes
#' @param n_ages Number of ages
#' @param age_bins Vector of ages
#' @param len_bins Vector of length 
#' @param LAA Dataframe with length-at-age data
#' @param LW Dataframe with weight-length data
#' @param sim Simulation number
#'
#' @return
#' @export
#'
#' @examples
get_biologicals = function(n_sexes, n_ages, age_bins, len_bins, LAA, LW, sim) {
  
  # Get WAA values from data
  waa_sex = matrix(0, ncol = n_sexes, nrow = n_ages)
  waa_nosex = matrix(0, ncol = 1, nrow = n_ages)
  laa_sex = matrix(0, ncol = n_sexes, nrow = n_ages)
  laa_nosex = matrix(0, ncol = 1, nrow = n_ages)
  al_matrix_sexsp = array(0, dim = c(c(length(age_bins), length(len_bins), n_sexes)))
  al_matrix_sexagg = array(0, dim = c(c(length(age_bins), length(len_bins), 1)))
  iter = sim
  
  for(s in 1:n_sexes) {
    
    LAA_obs = LAA[LAA$sim == iter & LAA$sex == s, ]
    WL_obs = LW[LW$sim == iter & LW$sex == s, ]

    # Get sex-specific WAA
    waa_sex_sp = get_WAA(LAA_obs_age = LAA_obs$ages, 
                         LAA_obs_len = LAA_obs$lens,
                         WL_obs_len = WL_obs$lens,
                         WL_obs_wt = WL_obs$wts,
                         ages = age_bins)
    
    waa_sex[,s] = waa_sex_sp[[1]] # waa
    laa_sex[,s] = waa_sex_sp[[2]] # laa
    al_matrix_sexsp[,,s] = get_al_trans_matrix(age_bins, len_bins, waa_sex_sp[[2]], waa_sex_sp[[6]]) # get al transition matrix
    
    if(s == n_sexes) { # sex-aggregated waa
      LAA_obs = LAA[LAA$sim == iter, ]
      WL_obs = LW[LW$sim == iter, ]
      waa_sex_agg = get_WAA(LAA_obs_age = LAA_obs$ages, 
                            LAA_obs_len = LAA_obs$lens,
                            WL_obs_len = WL_obs$lens,
                            WL_obs_wt = WL_obs$wts,
                            ages = age_bins)
      waa_nosex[,1] = waa_sex_agg[[1]] # get waa
      laa_nosex[,1] = waa_sex_agg[[2]] # laa
      al_matrix_sexagg[,,1] = get_al_trans_matrix(age_bins, len_bins, waa_sex_agg[[2]], waa_sex_agg[[6]])
    } # if sex-aggregated waa
  } # end s loop
  
  return(list(waa_nosex = waa_nosex, al_matrix_sexagg = al_matrix_sexagg, 
              waa_sex = waa_sex, al_matrix_sexsp = al_matrix_sexsp,
              laa_sex = laa_sex, laa_nosex = laa_nosex,
              sd_waa_sex_agg =  waa_sex_agg[[6]], 
              sd_waa_sex_sp_m =  waa_sex_sp[[6]]))
  
} # end function

#' Title Take additional newton steps with TMB model
#'
#' @param n.newton number of additional newton steps we want to take
#' @param ad_model MakeADFUN model object
#' @param mle_optim Optimized model object
#'
#' @return
#' @export
#'
#' @examples

add_newton <- function(n.newton, ad_model, mle_optim) {
  
  tryCatch(expr = for(i in 1:n.newton) {
    g = as.numeric(ad_model$gr(mle_optim$par))
    h = optimHess(mle_optim$par, fn = ad_model$fn, gr = ad_model$gr)
    mle_optim$par = mle_optim$par - solve(h,g)
    mle_optim$objective = ad_model$fn(mle_optim$par)
  }, error = function(e){e})
  
}

#' Title Run Estimation Model
#'
#' @param data list of data inputs
#' @param parameters list of parameter random starting values
#' @param map map to fix parameters
#' @param n.newton number of additional newton steps to take
#' @param iter.max number of iterations for nlminb to run (default = 1e5)
#' @param eval.max number of nlminb evaluations (default = 1e5)
#' @param silent whether or not to print stuff out
#' @param getsdrep whether to return a sdreport
#' @param DLL DLL character
#' @return
#' @export
#'
#' @examples

run_EM <- function(data, parameters, map, n.newton, random = NULL, DLL,
                   iter.max = 1e6, eval.max = 1e6, 
                   silent = TRUE, getsdrep = TRUE) {
  
  # Make AD Function here
  model_fxn <- TMB::MakeADFun(data, parameters, map, random = random,
                              DLL= DLL, silent = silent, 
                              checkParameterOrder = TRUE, tracepar = TRUE)
  
  # Optimize model here w/ nlminb
  mle_optim <- stats::nlminb(model_fxn$par, model_fxn$fn, model_fxn$gr, 
                             control = list(iter.max = iter.max, eval.max = eval.max))
  
  # Take additional newton steps
  add_newton(n.newton = n.newton, ad_model = model_fxn, mle_optim = mle_optim)
  
  if(getsdrep == TRUE) {
    # Get report with mle optimized parameters
    model_fxn$rep <- model_fxn$report(model_fxn$env$last.par.best) # Need to pass both fixed and random effects!!!
    # Get sd report here from TMB
    sd_rep <- TMB::sdreport(model_fxn)
  } # if get sdrep = TRUE
  
  return(list(model_fxn = model_fxn, mle_optim = mle_optim, sd_rep = sd_rep))
  
}
  
#' Title Run Model
#'
#' @param data Data in list format
#' @param parameters Parameters in list format
#' @param map Map parameters in list format 
#' @param DLL DLL
#' @param iter.max iterations for optimization 
#' @param eval.max evaluations for optimization
#' @param n.newton additional newton steps to take
#' @param silent Whether or not we want to output iterations
#' @param random Parameter for random effects 
#'
#' @return
#' @export
#'
#' @examples
run_model = function(data, parameters, map, DLL = "Sex_Str_EM", iter.max = 3e5, eval.max = 3e5, n.newton = 3,
                     silent = TRUE, random = NULL) {
  # make ad object
  model_fxn = TMB::MakeADFun(data, parameters, map, random = random, 
                             DLL= DLL, silent = silent,  checkParameterOrder = TRUE, tracepar = TRUE)
  
  # Optimize model here w/ nlminb
  mle_optim <- stats::nlminb(model_fxn$par, model_fxn$fn, model_fxn$gr, control = list(iter.max = iter.max, eval.max = eval.max))
  add_newton(n.newton = n.newton, ad_model = model_fxn, mle_optim = mle_optim) # take extra newton steps if needed
  model_fxn$sd_rep <- TMB::sdreport(model_fxn) # get standard errors from inverse hessian
  model_fxn$rep <- model_fxn$report(model_fxn$env$last.par.best) # Need to pass both fixed and random effects!!!
  return(model_fxn)
}


