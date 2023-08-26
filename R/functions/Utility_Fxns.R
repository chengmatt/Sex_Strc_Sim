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
