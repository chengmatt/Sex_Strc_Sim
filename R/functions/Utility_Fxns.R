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
  
  # Get midpoint of length bins to make sure there are still probabilities left 
  len_mids = len_bins[1:(length(len_bins) - 1)] + diff(len_bins) / 2
  # Construct age length matrix
  age_length = matrix(0.0, nrow = length(age_bins), ncol = length(len_mids))

  for(a in 1:length(age_bins)) {
    for(l in 2:length(len_bins)) {
      
      if (l == 2) { # Probability of being between 1st and 2nd length bin given age a
        age_length[a, l - 1] = pnorm(len_bins[2], mean_length[a], sd)
      } else if (l == length(len_bins)) { # Probability of being larger than the last length bin given age a
        age_length[a, l - 1] = 1 - pnorm(len_bins[length(len_bins) - 1], mean_length[a], sd)
      } else { # a of being in between length bins given age a
        age_length[a, l - 1] = pnorm(len_bins[l], mean_length[a], sd) -  
                               pnorm(len_bins[l - 1], mean_length[a], sd)
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
  return(1 / (1 + exp(-slope * (bins - midpoint)) ))
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
  return(1 / (1 + 19 ^ ((a50 - bins) / a95) ))
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
vonB = function(age_bins, k , L_inf, t0, sd) {
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
#' @param len_mids Vector of length midpoints
#' @param LAA Dataframe with length-at-age data
#' @param LW Dataframe with weight-length data
#' @param sim Simulation number
#'
#' @return
#' @export
#'
#' @examples
get_biologicals = function(n_sexes, n_ages, age_bins, len_mids, LAA, LW, sim) {
  
  # Get WAA values from data
  waa_sex = matrix(0, ncol = n_sexes, nrow = n_ages)
  waa_nosex = matrix(0, ncol = 1, nrow = n_ages)
  al_matrix_sexsp = array(0, dim = c(c(length(age_bins), length(len_mids), n_sexes)))
  al_matrix_sexagg = array(0, dim = c(c(length(age_bins), length(len_mids), 1)))
  
  for(s in 1:n_sexes) {
    
    # Get sex-specific WAA
    waa_sex_sp = get_WAA(LAA_obs_age = LAA$ages[LAA$sim == sim & LAA$sex == s], 
                         LAA_obs_len = LAA$lens[LAA$sim == sim & LAA$sex == s],
                         WL_obs_len = LW$lens[LW$sim == sim & LAA$sex == s],
                         WL_obs_wt = LW$wts[LW$sim == sim & LAA$sex == s],
                         ages = age_bins)
    
    waa_sex[,s] = waa_sex_sp[[1]] # waa
    al_matrix_sexsp[,,s] = get_al_trans_matrix(age_bins, len_bins, waa_sex_sp[[2]], waa_sex_sp[[6]]) # get al transition matrix
    
    if(s == n_sexes) { # sex-aggregated waa
      waa_sex_agg = get_WAA(LAA_obs_age = LAA$ages[LAA$sim == sim], 
                            LAA_obs_len = LAA$lens[LAA$sim == sim],
                            WL_obs_len = LW$lens[LW$sim == sim],
                            WL_obs_wt = LW$wts[LW$sim == sim],
                            ages = age_bins)
      waa_nosex[,1] = waa_sex_agg[[1]] # get waa
      al_matrix_sexagg[,,1] = get_al_trans_matrix(age_bins, len_bins, waa_sex_agg[[2]], waa_sex_agg[[6]])
    } # if sex-aggregated waa
  } # end s loop
  
  return(list(waa_nosex = waa_nosex, al_matrix_sexagg = al_matrix_sexagg, 
              waa_sex = waa_sex, al_matrix_sexsp = al_matrix_sexsp))
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
  
#' Title Scale to zero and one
#'
#' @param x Vector of values you want to rescale
#'
#' @return
#' @export
#'
#' @examples
scale_zero_one = function(x) {
  return((x - min(x))/(max(x) - min(x)))
}