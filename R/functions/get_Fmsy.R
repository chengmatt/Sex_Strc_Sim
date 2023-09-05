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

#' Title Get equilibrium SSB
#'
#' @param M  Female M
#' @param selex Female fishery selectivity
#' @param Trial_F Trial F
#' @param waa female weight at age
#' @param mat_at_age maturity at age
#' @param ages vector of ages
#' @param Init_N vector of inital ages
#' @return
#' @export
#'
#' @examples
get_SSBe = function(M, selex, Trial_F, waa, mat_at_age, ages, Init_N) {
  
  # set up
  ssbe = vector(length = length(ages))
  N_equilibrium = matrix(nrow = length(ages), ncol = length(ages) * 2)
  N_equilibrium[,1] = Init_N # initial numbers at N
  Fa = Trial_F * selex
  Za = Fa + M
  Sa = exp(-Za)
  
  for(y in 2:(length(ages) * 2)) { # population projeciton
    # Initial recruitment
    N_equilibrium[1, y] = Init_N[1]
    for(a in 2:length(ages)) N_equilibrium[a,y] = N_equilibrium[a-1,y-1] * Sa[a-1] # not plus group
    N_equilibrium[length(ages),y] = N_equilibrium[length(ages)-1,y-1] * Sa[a-1] +
                                    N_equilibrium[length(ages),y-1] * Sa[a] # plus group
  } # end y loop
  
  # Get equilibrium SSB (subject to fishing for 1 year and calculate ssb)
  for(a in 1:n_ages) ssbe[a] = N_equilibrium[a, (length(ages) * 2)] * Sa[a] * waa[a] * mat_at_age[a]  

  return(list(N_equilibrium = N_equilibrium, ssbe_sum = sum(ssbe), ssbe_vec = ssbe))
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
get_Fmsy_nLL = function(ln_Fmsy, M, selex, waa, mat_at_age, ages, Init_N) {
  Fmsy = exp(ln_Fmsy) # exponentiate
  # Get quantities to compute Fmsy
  SBPR = get_SBPR(M = M, selex = selex, Trial_F = Fmsy, waa = waa, mat_at_age = mat_at_age, ages = ages)
  YPR = get_YPR(M = M, selex = selex, Trial_F = Fmsy, waa = waa, ages = ages)
  Bmsy = get_SSBe(Init_N = Init_N,M = M, selex = selex, Trial_F = Fmsy, waa = waa, mat_at_age = mat_at_age, ages = ages)
  nLL = -1 * log((YPR$YPR_sum * Bmsy$ssbe_sum) / SBPR$SBPR_sum)
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

