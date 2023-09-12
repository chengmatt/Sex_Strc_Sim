# Purpose: Collection of functions to get the true fmsy from the OM
# Creator: Matthew LH. Cheng
# Date: 8/30/23

#' Title Get spawning biomass per recruit
#'
#' @param M  Female M
#' @param selex Female fishery selectivity
#' @param Trial_F Trial F
#' @param waa female weight at age
#' @param mat_at_age maturity at age
#' @param ages vector of ages
#'
#' @return
#' @export
#'
#' @examples

get_SBPR = function(M, selex, Trial_F, waa, mat_at_age, ages) {
  # Set up
  Na = vector(length = length(ages) * 4)
  SBPR = vector(length = length(ages) * 4)
  Za = M + (selex * Trial_F) # get total mortality
  Srva = exp(-Za) # get survival
  Na[1] = 1 # initialize population
  SBPR[1] = Na[1] * Srva[1] * waa[1] * mat_at_age[1]
  
  for(a in 2:(length(ages) * 4)) {
    if(a >= length(ages)) {  # if we are at thge max age of the iteartion, use the max age data to do calcs
      Na[a] = Na[a-1] * Srva[length(ages)]
      SBPR[a] = Na[a] * waa[length(ages)] * mat_at_age[length(ages)]
    } else{
      Na[a] = Na[a-1] * Srva[a]
      SBPR[a] = Na[a] * waa[a] * mat_at_age[a]
    } # end if else
  } # end a loop
  return(list(SBPR_Na = Na, SBPR_sum = sum(SBPR), SBPR_vec = SBPR))
} # end function



#' Title Get yield per recruit quantities
#'
#' @param M  Female M
#' @param selex Female fishery selectivity
#' @param Trial_F Trial F
#' @param waa female weight at age
#' @param ages vector of ages
#'
#' @return
#' @export
#'
#' @examples
get_YPR = function(M, selex, Trial_F, waa, ages) {
  
  # initialize population
  Na = vector(length = length(ages) * 4)
  YPR = vector(length = length(ages) * 4)
  Fa = selex * Trial_F # fishing mortality
  Za = M + (selex * Trial_F) # get total mortality
  Srva = exp(-Za) # get survival
  Na[1] = 1 # initialize population
  YPR[1] = (Na[1] * (1 - Srva[1]) * waa[1]) * (Fa[1] / Za[1]) # baranov
  
  for(a in 2:(length(ages) * 4)) {
    if(a >= length(ages)) {  # if we are at thge max age of the iteartion, use the max age data to do calcs
      Na[a] = Na[a-1] * Srva[length(ages)]
      YPR[a] = Na[a] * (1 - Srva[length(ages)]) * waa[length(ages)] * (Fa[length(ages)] / Za[length(ages)])
    } else{
      Na[a] = Na[a-1] * Srva[a]
      YPR[a] = Na[a] * (1 - Srva[a]) * waa[a] * (Fa[a] / Za[a])
    } # end if else
  } # end a loop
  return(list(YPR_Na = Na, YPR_sum = sum(YPR), YPR_vec = YPR))
} # end function

#' Title Get equilibrium recruitment
#'
#' @param SBPR_Fmsy SBPR Fmsy
#' @param SBPR_0 Unfished SBPR
#' @param waa weight at age
#' @param mat_at_age maturity at age
#' @param ages age bins
#'
#' @return
#' @export
#'
#' @examples
get_Req = function(SBPR_Fmsy, waa, mat_at_age, ages) {
  # Get unfished SBPR
  Z0 <- c(0, rep(M[1], length(ages)-1)) # Get Natural mortality
  N0 <- 1 * exp(-cumsum(Z0)) # Calculate Numbers over the lifespan of cohort
  SBPR_0 <- sum(N0 * waa * mat_at_age) # Get SSB in biomass units
  Req = (r0 * (4*h*SBPR_Fmsy - (1-h) * SBPR_0)) / ((5*h - 1) * SBPR_Fmsy) # get equilibrium recruitment
  return(Req)
}

#' Title Minimize to get Fmsy
#'
#' @param M  Female M
#' @param selex Female fishery selectivity
#' @param ln_Fmsy log Trial F
#' @param waa female weight at age
#' @param mat_at_age maturity at age
#' @param ages vector of ages
#' @param Init_N vector of inital ages
#' @return
#' @export
#'
#' @examples
get_Fmsy_nLL = function(ln_Fmsy, M, selex, waa, mat_at_age, ages, Init_N) {
  
  # set up optimization
  Fmsy = exp(ln_Fmsy) # exponentiate
  SBPR = get_SBPR(M = M, selex = selex, Trial_F = Fmsy, waa = waa, 
                  mat_at_age = mat_at_age, ages = ages) # get sbpr at fmsy
  YPR = get_YPR(M = M, selex = selex, Trial_F = Fmsy, waa = waa, ages = ages) # get YPR
  Req = get_Req(SBPR_Fmsy = SBPR$SBPR_sum, waa = waa, 
                mat_at_age = mat_at_age, ages = ages) # get equilibrium recruitment

  # set up optimization criteria
  catch_MSY = YPR$YPR_sum * Req 
  nLL = -1 * log(catch_MSY)
  
  return(nLL)
} # end function

#' Title Minimize to get Fmsy
#'
#' @param M  Female M
#' @param selex Female fishery selectivity
#' @param ln_Fmsy log Trial F
#' @param waa female weight at age
#' @param mat_at_age maturity at age
#' @param ages vector of ages
#' @param Init_N vector of inital ages
#' @return
#' @export
#'
#' @examples
get_Fmsy = function(ln_Fmsy, M, selex, waa, mat_at_age, ages, Init_N) {
  
  # Minimization to get Fmsy
  Fmsy = bbmle::mle2(get_Fmsy_nLL,
                     start = list(ln_Fmsy = ln_Fmsy),
                     data = list(M = M,
                                 selex = selex,
                                 waa = waa,
                                 mat_at_age = mat_at_age,
                                 ages = ages,
                                 Init_N = Init_N),
                     method="Nelder-Mead",
                     optimizer="nlminb",
                     control=list(maxit=1e6))
  
  return(list(Fmsy = exp(Fmsy@coef), obj = Fmsy@details$objective))
} # end function
