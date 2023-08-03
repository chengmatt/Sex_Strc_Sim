# Purpose: To simulate datasets for use in sex-structured simulations
# Author: Matthew LH. Cheng (UAF - CFOS)
# Date: 8/2/23

library(here)
library(tidyverse)

# General Arguments
n_sims = 5
n_years = 50
age_bins = 1:30
len_bins = seq(10, 85, by = 1)
len_mids = len_bins[1:(length(len_bins) - 1)] + diff(len_bins) / 2 # Get midpoint of lengths
n_sexes = 2
n_ages = length(age_bins)
n_lens = length(len_bins)

# Selex
len_slope = 0.2
len_midpoint = 36.5

# Von B (F, M)
k = c(0.2, 0.1) 
L_inf = c(80, 70)
t0 = c(-1.31, -1.31)
vonB_cv = 0.2

# Recruitment
sigmaRec = 0.8
# Sex ratio
sexRatio = c(1, 1)
# Mortality
M = 0.2

# Load in all functions from the functions folder
fxn_path <- here("R", "functions")
files <- list.files(fxn_path)
for(i in 1:length(files)) source(here(fxn_path, files[i]))

# Read and create OM objects
read_params_create_OM_objects(spreadsheet_path = here("input", "Sablefish_Inputs.xlsx"), n_years = 50)

# Construct vonB LAA 
# Females
vonB_Female = vonB(age_bins = age_bins, k = k[1], L_inf = L_inf[1], t0 = t0[1], cv = 0)
# Males
vonB_Male = vonB(age_bins = age_bins, k = k[2], L_inf = L_inf[2], t0 = t0[2], cv = 0)
plot(age_bins, vonB_Female, col = "red", type = "l")
lines(age_bins, vonB_Male, col = "blue", type = "l")

# Get age-length transition matrix
# Female
al_matrix_Female = get_al_trans_matrix(age_bins = age_bins, len_bins = len_bins,
                          mean_length = vonB_Female, cv = vonB_cv)
# Male
al_matrix_Male = get_al_trans_matrix(age_bins = age_bins, len_bins = len_bins,
                                       mean_length = vonB_Male, cv = vonB_cv)

# Combine matrices
al_matrix = array(c(al_matrix_Female, al_matrix_Male), 
                  dim = c(length(age_bins), length(len_mids), n_sex))

# Construct selectivity
# Length-based selectivity
len_selex = logist(slope = len_slope, bins = len_mids, midpoint = len_midpoint)
plot(len_mids, len_selex)

# Age-based selectivity converted from length-based selectivity
age_selex_Female = al_matrix[,,1] %*% len_selex
age_selex_Male = al_matrix[,,2] %*% len_selex
plot(age_bins, age_selex_Female, type = "l", col = "red", ylim = c(0,1))
lines(age_bins, age_selex_Male, type = "l", col = "blue")

NAA = array(0, dim = c(n_years, length(age_bins), n_sex, n_sims)) # Numbers at age
NAL = array(0, dim = c(n_years, length(len_mids), n_sex, n_sims)) # Numbers at length
SSB = array(0, dim = c(n_years, n_sims)) # Spawning stock biomass

# Start simulation
for(sim in 1:n_sims) {
  
  # Generate recruitment deviates and deviations from equilibrium
  Rec_Dev = exp(rnorm(n_years, mean = -sigma_rec^2/2, sd = sigma_rec))
  Init_Dev = exp(rnorm(length(age_bins), mean = -sigma_rec^2/2, sd = sigma_rec))
  
  # Initialize population first
  for(s in 1:n_sexes) {
    # Get numbers at age
    NAA[1,,s,sim] = r0 * exp(-M * 0:(length(age_bins)-1)) * Init_Dev * sexRatio[s]
    # Convert NAA to numbers at length
    NAL[1,,s,sim] = t(al_matrix[,,s]) %*% NAA[1,,s,1]
  } # end first sex loop
  
  # Calculate SSB
  SSB[1,sim] = sum(NAA[1,,1,sim] * wt_at_age[1,,1,sim] * mat_at_age[1,,1,sim])
  
  for(y in 2:n_years) {
    for(a in 1:n_ages) {
      for(s in 1:n_sexes) {
        
        # Recruitment Age
        if(a == 1) NAA[y,1,s,sim] = beverton_holt_recruit(ssb = SSB[y - 1, sim], 
                                        h = h,  r0 = r0, M = M[1]) * Rec_Dev[y]
        # Project Population Forward
        if(a > 1 & a < n_ages) NAA[y,a,s,sim] = NAA[y - 1,a - 1,s,sim] * exp(-M)
        if(a == n_ages) { # Calculate abundance for plus group (Mortality of MaxAge-1 + Mortality of MaxAge)
          NAA[y,a,s,sim] = (NAA[y - 1,a - 1,s,sim] * exp(-M)) + (NAA[y - 1,a,s,sim] * exp(-M))
          # Convert NAA to numbers at length
          NAL[y,,s,sim] = t(al_matrix[,,s]) %*% NAA[y,,s,1]
        } # end if for plus group
        
      } # end second sex loop
    } # end age loop
    
    # Calculate SSB for recruitment next year
    SSB[y,sim] = sum(NAA[y,,1,sim] * wt_at_age[y,,1,sim] * mat_at_age[y,,1,sim])
    
    
  } # end year loop
} # end sim loop


plot(NAA[50,,1,1], type = "l")
plot(NAL[50,,1,1], type = "l")

