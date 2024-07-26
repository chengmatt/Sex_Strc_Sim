#' Title Get projected catch from Fmsy and HCR (Tier 1 - 3 NPFMC control rule)
#'
#' @param n_ages Number of ages
#' @param n_sexes Number of sexes
#' @param fmsy Value for fmsy
#' @param bmsy_val Value for bmsy_val
#' @param term_NAA Array of terminal year NAA
#' @param term_SSB terminal SSB
#' @param term_F_Slx Array of terminal Fishery selex
#' @param term_F Vector of terminal year Fs
#' @param M_s Vector of natural mortality rates 
#' @param WAA Array of weight at age values
#' @param MatAA Vector of maturity at age
#' @param r0 Virgin recruitment
#' @return Numeric value of ABC
#' @export
#'
#' @examples
get_proj_catch <- function(fmsy_val,
                          bmsy_val,
                          sex_ratio,
                          n_ages, 
                          n_sexes, 
                          term_NAA,
                          term_SSB,
                          term_F_Slx,
                          term_F,
                          M_s,
                          WAA,
                          r0,
                          MatAA
                          ) {
  
  # Empty matrices to store projections in
  N_proj <- array(data = 0, dim = c(n_ages, n_sexes))
  F_proj <- array(data = 0, dim = c(n_ages, n_sexes))
  Z_proj <- array(data = 0, dim = c(n_ages, n_sexes))
  CAA_proj <- array(data = 0, dim = c(n_ages, n_sexes))
  Catch_proj <- array(data = 0, dim = c(n_ages, n_sexes))
  term_F_Slx <- array(term_F_Slx, dim = c(n_ages, n_sexes)) # preserve data dimensions
  term_NAA <- array(term_NAA, dim = c(n_ages, n_sexes)) # preserve data dimensions
  
  # project population forward for a year
  for(s in 1:n_sexes) {
    for(a in 1:n_ages) {
      # Input into N_proj - deterministic recruitment
      N_proj[1,] = beverton_holt_recruit(ssb = term_SSB, h = h, r0 = r0, M = M_s[1]) * sex_ratio[s]
      # Get terminal total mortality here (constant M)
      term_ZAA = sum(term_F * term_F_Slx[a,s], na.rm = TRUE) + M_s[s]
      if(a < n_ages) { # not plus group
        N_proj[a + 1, s] = term_NAA[a, s] * exp(-term_ZAA)
      } else{ 
        N_proj[a, s] = N_proj[a, s] + (term_NAA[a, s] * exp(-term_ZAA))
      } # else = plus group
    } # end a (age) loop
  } # end s (sex) loop
  
  if(n_sexes == 1) proj_SSB = sum(N_proj[,1] * 0.5 * WAA[,1] * MatAA[,1])
  if(n_sexes != 1) proj_SSB = sum(N_proj[,1] * WAA[,1] * MatAA[,1])
  
  # Harvest control rule here
  bmsy_val_ratio = proj_SSB / bmsy_val # get b/bmsy_val
  if(bmsy_val_ratio < 1) Fval = ((proj_SSB / bmsy_val  - 0.05) / (1 - 0.05)) * fmsy_val
  if(bmsy_val_ratio >= 1) Fval = fmsy_val
  if(bmsy_val_ratio <= 0.1) Fval = fmsy_val
  
  # Get f projections for age and sex here 
  for(a in 1:n_ages) {
    for(s in 1:n_sexes) {
      F_proj[a,s] = Fval * term_F_Slx[a,s]
      # Get mortality and catch projections now
      Z_proj[a,s] <- F_proj[a,s] + M_s[s]
      # Now, get catch projections
      CAA_proj[a,s] <- N_proj[a,s] * (1 - exp(-Z_proj[a,s])) * (F_proj[a,s] / Z_proj[a,s]) 
      Catch_proj[a,s] <- CAA_proj[a,s] * WAA[a,s] # Turn to biomass units
      ABC <- sum(Catch_proj) # sum to get abc
    } # s loop
  } # a loop
  
  
  
  
  return(c(ABC))
  
} # end function