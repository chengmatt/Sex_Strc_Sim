# Purpose: To predict new recruits following a beverton holt stock
# recruit relationship
# Creator: Matthew LH. Cheng
# Date: 10/29/22

#' @param ssb Spawning stock biomass 
#' @param h steepness (how steep you want the BH curve to be)
#' @param r0 Virgin recruitment (where recruitment asymptotes at)
#' @param M Natural Mortality for Females

beverton_holt_recruit <- function(ssb, h, r0, M) {
  
  # Get SPR
  SPR_N <- vector()
  SPR_SSB0 <- vector()
  
  for(a in 1:(n_ages * 2)) { # project 2x the number of ages in assessment
    if(a == 1) SPR_N[a] = 1
    if(a > 1) SPR_N[a] = SPR_N[a-1] * exp(-M) 
    if(a <= n_ages) SPR_SSB0[a] = SPR_N[a] * waa[a,1] * mat_at_age[a,1]
    if(a > n_ages) SPR_SSB0[a] = SPR_N[a] * waa[n_ages,1] * mat_at_age[n_ages,1]
  } # end a loop
  
  # Now, get SPR rate
  SPR0 <- sum(SPR_SSB0)
  ssb0 <<- SPR0 * r0
  
  # Output to environemnt
  SPR_SSB0 <<- SPR_SSB0
  
  # Get BH components
  BH_first_part <- 4 * h * r0 * ssb
  BH_sec_part <- (ssb0 * (1 - h)) + ssb * ((5*h) - 1)
  
  # Now calculate BH parameterization
  rec <- BH_first_part / BH_sec_part
  
  return( rec )
  
}
