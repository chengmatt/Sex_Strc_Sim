# Purpose: To specify fishing mortality scenarios (contrast, increase, constant)

#' Title
#'
#' @param F_pattern Contrast, Increase, or Constant for fishing mortality pattern
#' @param n_years Number of years
#' @param fmsy Fmsy value
#'
#' @return
#' @export
#'
#' @examples
f_pattern_scenarios = function(F_pattern, n_years, fmsy) {
  # Set up fishing mortality scenarios
  if(F_pattern == "Contrast") {
    midpoint = ((n_years) / 2) # get midpoint
    Fmort1 = seq(fmsy * 0.1, fmsy * 1, length.out = length(1:midpoint))
    Fmort2 = seq(fmsy * 0.98, fmsy * 0.5, length.out = length((midpoint+1):n_years))
    Fmort = c(Fmort1, Fmort2)
  } # end contrast fishing mortality pattern
  
  if(F_pattern == "Constant") {
    Fmort = rep(fmsy, n_years)
  } # end contrast fishing mortality pattern
  
  if(F_pattern == "Increase") {
    Fmort = seq(fmsy * 0.1, fmsy, length.out = n_years)
  } # end contrast fishing mortality pattern
  return(Fmort)
}