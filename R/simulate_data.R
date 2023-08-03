# Purpose: To simulate datasets for use in sex-structured simulations
# Author: Matthew LH. Cheng (UAF - CFOS)
# Date: 8/2/23

library(here)
library(tidyverse)

# General Arguments
n_sims = 5
n_years = 51
age_bins = 1:30
len_bins = seq(10, 85, by = 1)
len_mids = len_bins[1:(length(len_bins) - 1)] + diff(len_bins) / 2 # Get midpoint of lengths
n_sexes = 2
n_ages = length(age_bins)
n_lens = length(len_bins)
n_fish_fleets = 2
n_srv_fleets = 2

# Selex
fish_len_slope = 0.7
fish_len_midpoint = 50.5
srv_len_slope = 0.8
srv_len_midpoint = 40.5

# Von B (F, M)
k = c(0.3, 0.15) 
L_inf = c(75, 75)
t0 = c(-1.31, -1.31)
vonB_cv = 0.15

# Recruitment
sigmaRec = 0.8
# Sex ratio
sexRatio = c(1, 1)
# Mortality
M = c(0.2, 0.2)

# Load in all functions from the functions folder
fxn_path <- here("R", "functions")
files <- list.files(fxn_path)
for(i in 1:length(files)) source(here(fxn_path, files[i]))

# Read and create OM objects
read_params_create_OM_objects(spreadsheet_path = here("input", "Sablefish_Inputs.xlsx"), n_years = n_years)

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
                  dim = c(length(age_bins), length(len_mids), n_sexes))

# Construct selectivity
# Length-based selectivity - Fishery
fish_len_selex = logist(slope = fish_len_slope, bins = len_mids, midpoint = fish_len_midpoint)
srv_len_selex = logist(slope = srv_len_slope, bins = len_mids, midpoint = srv_len_midpoint)
plot(len_mids, fish_len_selex)
lines(len_mids, srv_len_selex, col = "blue")

# Age-based selectivity converted from length-based selectivity (fishery)
fish_age_selex_Female = al_matrix[,,1] %*% fish_len_selex
fish_age_selex_Male = al_matrix[,,2] %*% fish_len_selex
plot(age_bins, fish_age_selex_Female, type = "l", col = "red", ylim = c(0,1))
lines(age_bins, fish_age_selex_Male, type = "l", col = "blue")

# Age-based selectivity converted from length-based selectivity (survey)
srv_age_selex_Female = al_matrix[,,1] %*% srv_len_selex
srv_age_selex_Male = al_matrix[,,2] %*% srv_len_selex
plot(age_bins, srv_age_selex_Female, type = "l", col = "red", ylim = c(0,1))
lines(age_bins, srv_age_selex_Male, type = "l", col = "blue")

# Containers
NAA = array(0, dim = c(n_years, length(age_bins), n_sexes, n_sims)) # Numbers at age
NAL = array(0, dim = c(n_years, length(len_mids), n_sexes, n_sims)) # Numbers at length
FAA = array(0, dim = c(n_years, length(age_bins), n_sexes, n_fish_fleets, n_sims)) # Fishery Mortality at Age
ZAA = array(0, dim = c(n_years, length(age_bins), n_sexes, n_sims)) # Total Mortality at Age
CAA = array(0, dim = c(n_years, length(age_bins), n_sexes, n_fish_fleets, n_sims)) # Catch at Age
CAL = array(0, dim = c(n_years, length(len_mids), n_sexes, n_fish_fleets, n_sims)) # Catch at Length
Total_Catch_Sex = array(0, dim = c(n_years, n_sexes, n_fish_fleets, n_sims)) # Sex-Specific Catch
Total_Catch = array(0, dim = c(n_years, n_fish_fleets, n_sims)) # Aggregated Catch
SSB = array(0, dim = c(n_years, n_sims)) # Spawning stock biomass

# Fishing Mortality
Fmort = array(cumsum(rep(0.001, n_years)), dim = c(n_years, n_fish_fleets, n_sims))
# Fishery Age Selectivity
FishAge_Selex = array(c(fish_age_selex_Female, fish_age_selex_Male),dim = c(length(age_bins), n_sexes, n_fish_fleets))
# Fishery Length Selectivity
FishLen_Selex = array(c(fish_len_selex, fish_len_selex),dim = c(length(len_mids), n_sexes, n_fish_fleets))
# Fishery Age Comps
Fish_AgeComps = array(0, dim = c(n_years, length(age_bins), n_sexes, n_fish_fleets, n_sims)) 
Fish_Neff_Age = array(100, dim = c(n_years, n_fish_fleets)) # Age Effective Sample Size
# Fishery Length Comps
Fish_LenComps = array(0, dim = c(n_years, length(len_mids), n_sexes, n_fish_fleets, n_sims)) 
Fish_Neff_Len= array(300, dim = c(n_years, n_lens)) # Length Effective Sample Size
# Fishery Index
Fish_Index = array(0, dim = c(n_years, n_fish_fleets, n_sims))
q_Fish = c(0.03, 0.05)
cv_Fish_Index = c(0.1, 0.1)

# Survey Selex
SrvLen_Selex = array(c(srv_len_selex, srv_len_selex),dim = c(length(len_mids), n_sexes, n_srv_fleets))
SrvAge_Selex = array(c(srv_age_selex_Female, srv_age_selex_Male),dim = c(length(age_bins), n_sexes, n_srv_fleets))
# Survey Age Comps
Srv_AgeComps = array(0, dim = c(n_years, length(age_bins), n_sexes, n_srv_fleets, n_sims)) 
Srv_Neff_Age = array(100, dim = c(n_years, n_srv_fleets)) # Age Effective Sample Size
# Survey Length Comps
Srv_LenComps = array(0, dim = c(n_years, length(len_mids), n_sexes, n_srv_fleets, n_sims)) 
Srv_Neff_Len= array(300, dim = c(n_years, n_lens)) # Length Effective Sample Size
# Survey Index
Srv_Index = array(0, dim = c(n_years, n_srv_fleets, n_sims))
q_Srv = c(0.03, 0.05)
cv_Srv_Index = c(0.1, 0.1)

# Start simulation
for(sim in 1:n_sims) {
  
  # Initialize Population --------------------------------------------------

  # Generate recruitment deviates and deviations from equilibrium
  Rec_Dev = exp(rnorm(n_years, mean = -sigma_rec^2/2, sd = sigma_rec))
  Init_Dev = exp(rnorm(length(age_bins), mean = -sigma_rec^2/2, sd = sigma_rec))
  
  for(s in 1:n_sexes) {
    # Get numbers at age
    NAA[1,,s,sim] = r0 * exp(-M[s] * 0:(length(age_bins)-1)) * Init_Dev * sexRatio[s]
    # Convert NAA to numbers at length
    NAL[1,,s,sim] = t(al_matrix[,,s]) %*% NAA[1,,s,1]
  } # end first sex loop
  
  # Calculate SSB at time t = 1
  SSB[1,sim] = sum(NAA[1,,1,sim] * wt_at_age[1,,1,sim] * mat_at_age[1,,1,sim])
  
  for(y in 2:n_years) {
    # Calculate Deaths from Fishery
    FAA[y - 1,,,,sim] = Fmort[y-1,,sim] * FishAge_Selex[,,] # Fishing Mortality at Age
    for(a in 1:n_ages) {
      for(s in 1:n_sexes) {
        # Calculate Total Mortality
        ZAA[y - 1,a,s,sim] = M[s] + sum(FAA[y - 1,a,s,,sim])

# Project Population Forward ----------------------------------------------

        # Recruitment
        if(a == 1) {
          NAA[y,1,s,sim] = beverton_holt_recruit(ssb = SSB[y - 1, sim], 
                                                 h = h,  r0 = r0, M = M[1]) * Rec_Dev[y]
        } # end if recruitment age
        
        if(a > 1 & a < n_ages) {
          NAA[y,a,s,sim] = NAA[y - 1,a - 1,s,sim] * exp(-ZAA[y - 1,a - 1,s,sim])
        } # end if between recruitment age and plus group
        
        # Calculate abundance for plus group (Mortality of MaxAge-1 + Mortality of MaxAge)
        if(a == n_ages) { 
          NAA[y,a,s,sim] = (NAA[y - 1,a - 1,s,sim] * exp(-ZAA[y - 1,a - 1,s,sim])) + 
                           (NAA[y - 1,a,s,sim] * exp(-ZAA[y - 1,a,s,sim]))
          
          # Convert NAA to numbers at length
          NAL[y,,s,sim] = t(al_matrix[,,s]) %*% NAA[y,,s,1]
        } # end if for plus group
      } # end second sex loop
    } # end age loop
    
    # Calculate SSB here
    SSB[y,sim] = sum(NAA[y,,1,sim] * wt_at_age[y,,1,sim] * mat_at_age[y,,1,sim])
    
# Observation Model (Fishery) -------------------------------------------------------
    for(f in 1:n_fish_fleets) {
      for(s in 1:n_sexes) {
        
        # Get catch at age
        CAA[y-1,,s,f,sim] = (FAA[y-1,,s,f,sim] / ZAA[y - 1,,s,sim]) * 
                             (NAA[y-1,,s,sim] * (1 - exp(-ZAA[y - 1,,s,sim])))
        # Get Catch at Length
        CAL[y-1,,s,f,sim] = t(al_matrix[,,s]) %*% CAA[y-1,,s,f,sim]
        # Get Total Catch by Sex
        Total_Catch_Sex[y-1,s,f,sim] = CAA[y-1,,s,f,sim] %*% wt_at_age[y-1,,s,sim]
        
        ### Simulate Fishery Compositions ---------------------------------------
        
        # Fishery Age Compositions
        Fish_AgeComps[y-1,,s,f,sim] = rmultinom(1, size = Fish_Neff_Age[y,f], 
                                                prob = CAA[y-1,,s,f,sim]/sum(CAA[y-1,,s,f,sim]))
        # Fishery Length Compositions
        Fish_LenComps[y-1,,s,f,sim] = rmultinom(1, size = Fish_Neff_Len[y,f], 
                                                prob = CAL[y-1,,s,f,sim]/sum(CAL[y-1,,s,f,sim]))
        
      } # end third sex loop

      ### Simulate Fishery Index --------------------------------------------------
      # Calculate sd for deviations
      Fish_Index_sd = sqrt(log(cv_Fish_Index[f]^2+1))
      # Sample Fishery Index with lognormal error
      Fish_Index[y-1,f,sim] = sum(q_Fish[f] * FishAge_Selex[,,f] * NAA[y-1,,,sim]) * 
                              exp(rnorm(1, -Fish_Index_sd^2/2, Fish_Index_sd))

      # Get Total Catch
      Total_Catch[y-1,f,sim] = sum(Total_Catch_Sex[y-1,,f,sim])
      
    } # end fish fleet loop
    
# Observation Model (Survey) ----------------------------------------------
    for(sf in 1:n_srv_fleets) {
      for(s in 1:n_sexes) {
        
        ### Simulate Survey Compositions ---------------------------------------
        # Compute probability of age-selection
        Prob_Age_Selex = NAA[y-1,,s,sim] * SrvAge_Selex[,s,sf] / sum(NAA[y-1,,s,sim] * SrvAge_Selex[,s,sf])
        Prob_Len_Selex = NAL[y-1,,s,sim] * SrvLen_Selex[,s,sf] / sum(NAL[y-1,,s,sim] * SrvLen_Selex[,s,sf])
        
        # Survey Age Compositions
        Srv_AgeComps[y-1,,s,sf,sim] = rmultinom(1, size = Srv_Neff_Age[y,sf], 
                                                prob = Prob_Age_Selex)
        # Survey Length Compositions
        Srv_LenComps[y-1,,s,sf,sim] = rmultinom(1, size = Srv_Neff_Len[y,sf], 
                                                prob = Prob_Len_Selex)
         
      } # end fourth sex loop
      
      ### Simulate Survey Index --------------------------------------------------
      # Calculate sd for deviations
      Srv_Index_sd = sqrt(log(cv_Srv_Index[sf]^2+1))
      # Sample Fishery Index with lognormal error
      Srv_Index[y-1,sf,sim] = sum(q_Srv[sf] * SrvAge_Selex[,,sf] * NAA[y-1,,,sim]) * 
                              exp(rnorm(1, -Srv_Index_sd^2/2, Srv_Index_sd))
      
    } # end survey fleet loop
  } # end year loop
} # end sim loop



plot(SSB[,1], type = "l")
plot(CAL[50,,1,1,1], type = 'l')
plot(Total_Catch[,2,1])
plot(Fish_Index[-51,1,1], type = "l")
plot(Srv_Index[-51,1,1], type = "l")
plot(Fish_LenComps[5,,1,1,1], type = "l")
plot(Srv_LenComps[1,,1,1,1], type = "l")
