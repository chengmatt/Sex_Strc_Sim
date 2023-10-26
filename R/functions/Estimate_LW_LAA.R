# Purpose: To estimate length weight and length-at-age relationships
# Author: Matthew LH. Cheng (UAF-CFOS)
# Date 8/23/23


#' Title Predict weights based on observed lengths
#'
#' @param length Vector of lengths
#' @param alpha condition factor/average tissue density
#' @param beta allometric scaling parameter
#'
#' @return
#' @export
#'
#' @examples
pred_wl = function(length, alpha, beta) {
  weight = alpha * length^beta
  return(weight)
} # end function


#' Title Predict lengths based on observed ages
#'
#' @param age Vector of ages
#' @param linf Asymptotic length
#' @param k Brody growth coefficient
#' @param t0 Theoretical length at age0
#'
#' @return
#' @export
#'
#' @examples
pred_laa = function(age, linf, k, t0) {
  len = linf * (1 - exp(-k * (age - t0)))
  return(len)
} # end function


#' Title Minimze to get parameters for weight length relationship
#'
#' @param obs_len Vector of observed lengths
#' @param obs_wt Vector of observed weights
#' @param ln_alpha log alpha
#' @param ln_beta log beta
#' @param ln_sigma log sigma
#'
#' @return
#' @export
#'
#' @examples
nLL_WL = function(obs_len, obs_wt, ln_alpha, ln_beta, ln_sigma) {
  
  # Exponentiate parameters
  alpha = exp(ln_alpha)
  beta = exp(ln_beta)
  sigma = exp(ln_sigma)
  
  # Get predicted weights
  pred_wts = pred_wl(length = obs_len, alpha = alpha, beta = beta)
  nLL = -1 * sum(dnorm(obs_wt, pred_wts, sd = sigma, log = TRUE)) # Get nLL
  return(nLL)
} # end function

#' Title Get length-at-age relationship and minimze parameters
#'
#' @param obs_lens vector of observed lengths
#' @param obs_age vector of observed ages
#' @param ln_linf asymptotic length in log
#' @param ln_k log brody coeff
#' @param t0 theoretical length at age 0
#' @param ln_sigma log sigma
#'
#' @return
#' @export
#'
#' @examples
nLL_LAA = function(obs_lens, obs_age, ln_linf, ln_k, t0, ln_sigma) {
  
  # Exponentiate parameters
  linf = exp(ln_linf)
  k = exp(ln_k)
  sigma = exp(ln_sigma)
  
  # Get predicted lengths
  pred_lens = pred_laa(age = obs_age, linf = linf, k = k, t0 = t0)
  nLL = -1 * sum(dnorm(obs_lens, pred_lens, sd = sigma, log = TRUE)) # Get nLL
  return(nLL)
} # end function

#' Title Function to estimate WL and LAA relationship, and then turn into vector of WAA
#'
#' @param LAA_obs_age LAA observed ages
#' @param LAA_obs_len LAA observed lengths
#' @param WL_obs_len WL observed lengths
#' @param WL_obs_wt WL observed weights
#' @param ages vector of ages
#'
#' @return
#' @export
#'
#' @examples
get_WAA = function(LAA_obs_age, LAA_obs_len, WL_obs_len, WL_obs_wt, ages) {
  
  # Get weight-length relationship
  WL_rel = bbmle::mle2(nLL_WL,
                        start = list(ln_alpha = log(1e-05),
                                     ln_sigma = log(3)),
                        data = list(obs_len = WL_obs_len,
                                    obs_wt = WL_obs_wt,
                                    ln_beta = log(beta_wl[1])), # fix beta
                        method="Nelder-Mead",
                        optimizer="nlminb",
                        control=list(maxit=5e6))
  
  
  # Get length at age relationship
  LAA_rel = bbmle::mle2(nLL_LAA,
                        start = list(ln_linf = log(80),
                                     ln_k = log(0.2),
                                     ln_sigma = log(2),
                                     t0 = 0),
                        data = list(obs_age = LAA_obs_age,
                                    obs_lens = LAA_obs_len),
                        method="Nelder-Mead",
                        optimizer="nlminb",
                        control=list(maxit=5e6))
  
  # Extract parameters
  linf = exp(LAA_rel@coef[1])
  k = exp(LAA_rel@coef[2])
  t0 = LAA_rel@coef[3]
  laa_sd = exp(LAA_rel@coef[4])
  alpha = exp(WL_rel@coef[1])
  beta = beta_wl[1] # fixing beta here
  
  # Get weight at age now
  winf = (alpha * linf^beta)
  waa = winf * (1 - exp(-k * (ages - t0)))^beta
  laa = linf * (1 - exp(-k * (ages - t0)))
  return(list(waa, laa, linf, k, t0, laa_sd, alpha, beta))
} # end function
