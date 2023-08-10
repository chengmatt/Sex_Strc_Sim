# Purpose: To simulate datasets for use in sex-structured simulations
# Author: Matthew LH. Cheng (UAF - CFOS)
# Date: 8/2/23

library(here)
library(tidyverse)

# General Arguments
n_sims = 5
n_years = 51
age_bins = 1:30
len_bins = seq(20, 85, by = 1)
len_mids = len_bins[1:(length(len_bins) - 1)] + diff(len_bins) / 2 # Get midpoint of lengths
n_sexes = 2
n_ages = length(age_bins)
n_lens = length(len_mids)
n_fish_fleets = 2
n_srv_fleets = 2

# Selex
fish_len_slope = 0.7
fish_len_midpoint = 30.5
srv_len_slope = 0.8
srv_len_midpoint = 25.5

# Von B (F, M)
k = c(0.15, 0.125) 
L_inf = c(80, 75)
t0 = c(-1.31, -1.31)
vonB_cv = 0.15

# Length-weight
alpha = c(0.00002, 0.00001)
beta = c(3, 3)
lw_cv = 0.15

# Recruitment
sigmaRec = 0.5
# Sex ratio
sexRatio = c(0.5, 0.5)
# Mortality
M = c(0.15, 0.15)

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
vonB = matrix(c(vonB_Female, vonB_Male), ncol = 2)
plot(age_bins, vonB_Female, col = "red", type = "l")
lines(age_bins, vonB_Male, col = "blue", type = "l")

# Construct LW 
lw_female = alpha[1] * len_bins^beta[1]
lw_male = alpha[2] * len_bins^beta[2]
lw = matrix(c(lw_female, lw_male), ncol = 2)
plot(len_bins, lw_female, col = "red", type = "l")
lines(len_bins, lw_male, col = "blue", type = "l")

# Construct WAA relationship
waa_female = alpha[1] * vonB_Female^beta[1]
waa_male = alpha[2] * vonB_Female^beta[2]
waa = matrix(c(waa_female, waa_male), ncol = 2)
plot(age_bins, waa_female, col = "red", type = "l")
lines(age_bins, waa_male, col = "blue", type = "l")


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

# Age-Lenght matrix unsexed
al_matrix_unsexed = get_al_trans_matrix(age_bins = age_bins, len_bins = len_bins,
                                        mean_length = rowMeans(matrix(c(vonB_Female, vonB_Male), ncol = 2)), 
                                        cv = vonB_cv)

# Construct selectivity
# Length-based selectivity - Fishery
fish_len_selex = logist(slope = fish_len_slope, bins = len_mids, midpoint = fish_len_midpoint)
srv_len_selex = logist(slope = srv_len_slope, bins = len_mids, midpoint = srv_len_midpoint)
plot(len_mids, fish_len_selex)
lines(len_mids, srv_len_selex, col = "blue")

# Age-based selectivity converted from length-based selectivity (fishery)
fish_age_selex_Female = al_matrix[,,1] %*% fish_len_selex
fish_age_selex_Male = al_matrix[,,2] %*% fish_len_selex
# fish_age_selex_Female = logist(slope = 0.7, bins = age_bins, midpoint = 4)
# fish_age_selex_Male = logist(slope = 0.3, bins = age_bins, midpoint = 8)
plot(age_bins, fish_age_selex_Female, type = "l", col = "red", ylim = c(0,1))
lines(age_bins, fish_age_selex_Male, type = "l", col = "blue")

# Age-based selectivity converted from length-based selectivity (survey)
srv_age_selex_Female = al_matrix[,,1] %*% srv_len_selex
srv_age_selex_Male = al_matrix[,,2] %*% srv_len_selex
# srv_age_selex_Female = logist(slope = 3, bins = age_bins, midpoint = 2)
# srv_age_selex_Male = logist(slope = 0.5, bins = age_bins, midpoint = 5)
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
Fmort = array(seq(0.1, 0.001, length.out = n_years), dim = c(n_years, n_fish_fleets, n_sims))
# Fishery Age Selectivity
FishAge_Selex = array(c(fish_age_selex_Female, fish_age_selex_Male),dim = c(length(age_bins), n_sexes, n_fish_fleets))
# Fishery Length Selectivity
# FishLen_Selex = array(c(fish_len_selex, fish_len_selex),dim = c(length(len_mids), n_sexes, n_fish_fleets))
# Fishery Age Comps
Fish_AgeComps = array(0, dim = c(n_years, length(age_bins), n_sexes, n_fish_fleets, n_sims)) 
Fish_Neff_Age = array(100, dim = c(n_years, n_fish_fleets)) # Age Effective Sample Size
# Fishery Length Comps
Fish_LenComps = array(0, dim = c(n_years, length(len_mids), n_sexes, n_fish_fleets, n_sims)) 
Fish_Neff_Len= array(100, dim = c(n_years, n_fish_fleets)) # Length Effective Sample Size
# Fishery Length-at-age
Fish_LAA = list()
# Fishery Length-Weight
Fish_LW = list()
# Fishery Index
Fish_Index = array(0, dim = c(n_years, n_fish_fleets, n_sims))
q_Fish = c(0.03, 0.05)
cv_Fish_Index = c(0.25, 0.25)

# Survey Selex
# SrvLen_Selex = array(c(srv_len_selex, srv_len_selex),dim = c(length(len_mids), n_sexes, n_srv_fleets))
SrvAge_Selex = array(c(srv_age_selex_Female, srv_age_selex_Male),dim = c(length(age_bins), n_sexes, n_srv_fleets))
# Survey Age Comps
Srv_AgeComps = array(0, dim = c(n_years, length(age_bins), n_sexes, n_srv_fleets, n_sims)) 
Srv_Neff_Age = array(100, dim = c(n_years, n_srv_fleets)) # Age Effective Sample Size
# Survey Length Comps
Srv_LenComps = array(0, dim = c(n_years, length(len_mids), n_sexes, n_srv_fleets, n_sims)) 
Srv_Neff_Len= array(100, dim = c(n_years, n_lens)) # Length Effective Sample Size
# Survey Length-at-age
Srv_LAA = list()
# Survey Length-Weight
Srv_LW = list()
# Survey Index
Srv_Index = array(0, dim = c(n_years, n_srv_fleets, n_sims))
q_Srv = c(0.03, 0.05)
cv_Srv_Index = c(0.25, 0.25)

# Counters
fish_counter_age = 1
fish_counter_len = 1
srv_counter_age = 1
srv_counter_len = 1

RecDevs = array(0, dim = c(n_years-1, 2))
InitDevs = array(0, dim = c(length(age_bins), 2))

# Start simulation
for(sim in 1:2) {
  
  # Initialize Population --------------------------------------------------
  # Generate recruitment deviates and deviations from equilibrium
  RecDevs[,sim] = exp(rnorm(n_years-1, mean = -sigma_rec^2/2, sd = sigma_rec))
  InitDevs[,sim] = exp(rnorm(length(age_bins), mean = -sigma_rec^2/2, sd = sigma_rec))
  
  for(s in 1:n_sexes) {
    # Get numbers at age - not plus group
    NAA[1,1:length(age_bins),s,sim] = r0 * exp(-M[s] * 0:(length(age_bins)-1)) * InitDevs[,sim] * sexRatio[s]
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
                                                 h = h,  r0 = r0, M = M[1]) * RecDevs[y - 1,sim] * sexRatio[s]
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
      } # end third sex loop
      
      ### Simulate Fishery Compositions ---------------------------------------
      # Getting numbers at age as probability across sexes (this is more reflective
      # of a real world scenario; doing it as probabilities within sexes assumes equal
      # probability of selecting a particular sex, which is likely not true)
      # Get Probability of sampling age classes
      Prob_FishAge = as.vector(CAA[y-1,,,f,sim]/sum(CAA[y-1,,,f,sim]))
      Fish_AgeComps[y-1,,,f,sim] = matrix(rmultinom(1, size = Fish_Neff_Age[y,f] * n_sexes, 
                                                    prob = Prob_FishAge), nrow = n_ages, ncol = n_sexes)
      
      # Fishery Length Compositions
      # Get Probability of sampling length bins
      Prob_FishLen = as.vector(CAL[y-1,,,f,sim]/sum(CAL[y-1,,,f,sim]))
      Fish_LenComps[y-1,,,f,sim] = matrix(rmultinom(1, size = Fish_Neff_Len[y,f] * n_sexes, 
                                                    Prob_FishLen), nrow = n_lens, ncol = n_sexes)
      
      ### Simulate Fishery Index --------------------------------------------------
      # Calculate sd for deviations
      Fish_Index_sd = sqrt(log(cv_Fish_Index[f]^2+1))
      # Sample Fishery Index with lognormal error
      Fish_Index[y-1,f,sim] = q_Fish[f] * sum(FishAge_Selex[,,f] * NAA[y-1,,,sim] * wt_at_age[y,,,sim]) * 
                              exp(rnorm(1, -Fish_Index_sd^2/2, Fish_Index_sd))

      # Get Total Catch
      Total_Catch[y-1,f,sim] = sum(Total_Catch_Sex[y-1,,f,sim])
    } # end fish fleet loop
    
# Observation Model (Survey) ----------------------------------------------
    for(sf in 1:n_srv_fleets) {
      for(s in 1:n_sexes) {
        
        ### Simulate Growth from Survey Age and Length Compositions ---------------------------
        # Assuming random sampling here
        # Get length-at-age samples
        # for (a in 1:n_ages) {
        #   if (Srv_AgeComps[y-1,a,s,sf,sim] != 0) {
        #     n_age_samples = Srv_AgeComps[y-1,a,s,sf,sim]
        #     sampled_srv_lens = rnorm(n = n_age_samples, mean = vonB[a,s], sd = vonB_cv)
        #     Srv_LAA[[srv_counter_age]] = data.frame(lens = sampled_srv_lens, ages = a, 
        #                                          sex = s, srv_fleet = sf, sim = sim)
        #     srv_counter_age = srv_counter_age + 1
        #   } # end if statement
        # } # end third age loop
        # 
        # # Get length-weight samples
        # for(l in 1:n_lens) {
        #   if (Srv_LenComps[y-1,l,s,f,sim] != 0) {
        #     n_len_samples = Srv_LenComps[y-1,l,s,sf,sim]
        #     sampled_srv_wts = rnorm(n = n_len_samples, mean = lw[l,s], sd = lw_cv)
        #     Srv_LW[[srv_counter_len]] = data.frame(wts = sampled_srv_wts, lens = l, 
        #                                         sex = s, srv_fleet = sf, sim = sim)
        #     srv_counter_len = srv_counter_len + 1
        #   } # end if statement
        # } # end second length loop
         
      } # end fourth sex loop
      
      ### Simulate Survey Compositions ---------------------------------------
      # Again, we are sampling from a multinomial as proportions across sexes, given
      # that the probabilities of selecting a given age and sex are influenced by
      # the underlying sex-ratio of the population (i.e., if there are more females than males,
      # we should techicnally be sampling more females just by random chance).
      
      # Survey Age Compositions
      Prob_SrvAge = (NAA[y-1,,,sim] * SrvAge_Selex[,,sf]) / sum(NAA[y-1,,,sim] * SrvAge_Selex[,,sf]) # Get probability of sampling ages
      Srv_AgeComps[y-1,,,sf,sim] = matrix(rmultinom(1, size = Srv_Neff_Age[y,sf] * n_sexes,  # sampling 
                                                   as.vector(Prob_SrvAge)), nrow = n_ages, ncol = n_sexes)

      # Survey Length Compositions
      # Get probability of sampling lengths
      Prob_SrvLen = vector() # re-initialize vector
      for(s in 1:n_sexes) {
        Prob_SrvLen_s = (t(al_matrix[,,s]) %*% Prob_SrvAge[,s]) 
        Prob_SrvLen = rbind(Prob_SrvLen, Prob_SrvLen_s)
      } # end fifth sex loop
      
      Prob_Len = Prob_SrvLen / sum(Prob_SrvLen) # compute prob of len sampling
      Srv_LenComps[y-1,,,sf,sim] = matrix(rmultinom(1, size = Srv_Neff_Len[y,sf] * n_sexes, # sampling
                                                    Prob_SrvLen), nrow = n_lens, ncol = n_sexes)
      
      ### Simulate Survey Index --------------------------------------------------
      # Calculate sd for deviations
      Srv_Index_sd = sqrt(log(cv_Srv_Index[sf]^2+1))
      # Sample Fishery Index with lognormal error
      Srv_Index[y-1,sf,sim] = q_Srv[sf] * sum(SrvAge_Selex[,,sf] * NAA[y-1,,,sim]) * 
                              exp(rnorm(1, -Srv_Index_sd^2/2, Srv_Index_sd))
      
    } # end survey fleet loop
  } # end year loop
  print(sim)
} # end sim loop

plot(SSB[,1], type = "l")
plot(CAL[n_years-1,,1,1,1], type = 'l')
plot(Total_Catch[-n_years,2,1])
plot(Fish_Index[-n_years,1,1], type = "l")
plot(Srv_Index[-n_years,1,1], type = "l")
plot(Fish_AgeComps[3,,1,1,1], type = "l")
plot(Srv_LenComps[20,,1,1,1], type = "l")


# Get Length weight samples
# Srv_LAA = data.table::rbindlist(Srv_LAA)
# Srv_LW = data.table::rbindlist(Srv_LW)
# 
# # Plot data
# ggplot(Srv_LAA %>% filter(sim == 2), aes(x = ages, y = lens, color = sex)) +
#   geom_point() 
# ggplot(Srv_LW %>% filter(sim == 1), aes(x = lens, y = wts, color = sex)) +
#   geom_point() 


# TMB Testing -------------------------------------------------------------

data = list(
  
  # Controls
  years = 1:(n_years-1),
  ages = 1:n_ages,
  lens = 1:n_lens,
  n_sexes = n_sexes,
  n_fish_fleets = n_fish_fleets,
  n_srv_fleets = n_srv_fleets,

  # Fishery
  obs_catch = as.matrix(Total_Catch[-n_years,,1], dim = c(50, 2)),
  catch_cv = c(1e-3, 1e-3),
  obs_fish_index = array(Fish_Index[-n_years,,1], dim = c(50, 2)),
  fish_index_cv = cv_Fish_Index,
  obs_fish_age_comps = array(Fish_AgeComps[-n_years,,,,1], dim = c(50, 30, 2, 2)),
  fish_age_comps_inputN = array(Fish_Neff_Age[-n_years,], dim = c(50, 2)),
  obs_fish_len_comps = array(Fish_LenComps[-n_years,,,,1], dim = c(50, 65, 2, 2)),
  fish_len_comps_inputN = array(Fish_Neff_Len[-n_years,], dim = c(50, 2)),
  
  # Survey
  obs_srv_index = array(Srv_Index[-n_years,,1], dim = c(50, 2)),
  srv_index_cv = cv_Srv_Index,
  obs_srv_age_comps = array(Srv_AgeComps[-n_years,,,,1], dim = c(50, 30, 2, 2)),
  srv_age_comps_inputN = array(Srv_Neff_Age[-n_years,], dim = c(50, 2)),
  obs_srv_len_comps = array(Srv_LenComps[-n_years,,,,1], dim = c(50, 65, 2, 2)),
  srv_len_comps_inputN = array(Srv_Neff_Len[-n_years,], dim = c(50, 2)),
  
  # Biologicals
  sexRatio = sexRatio,
  WAA = matrix(wt_at_age[1,,,1], ncol = 2, nrow = 30),
  MatAA = as.vector(mat_at_age[1,,1,1]),
  age_len_transition = al_matrix,
  age_len_transition_unsexed = al_matrix_unsexed,
  
  # Data Indicators
  use_catch = matrix(1, ncol = 2, nrow = 50),
  use_fish_index = matrix(1, ncol = 2, nrow = 50),
  use_srv_index = matrix(1, ncol = 2, nrow = 50),
  use_fish_age_comps = matrix(1, ncol = 2, nrow = 50),
  use_fish_len_comps = matrix(1, ncol = 2, nrow = 50),
  use_srv_age_comps = matrix(1, ncol = 2, nrow = 50),
  use_srv_len_comps = matrix(1, ncol = 2, nrow = 50),
  p_ow_sex_fish_age = 1,
  p_ow_sex_fish_len = 1,
  p_ow_sex_srv_age = 1,
  p_ow_sex_srv_len = 1,
  agg_sex_fish_age = 0,
  agg_sex_srv_age = 0
)

# Parameters
parameters = list(
  ln_M = log(M) * 0.5,
  ln_InitDevs = log(InitDevs[,1]),
  ln_RecDevs = log(RecDevs[-50,1]),
  RecPars = c(log(r0) * 0.5, h),
  ln_sigmaRec = log(sigma_rec),
  ln_q_fish = log(q_Fish),
  ln_q_srv = log(q_Srv),
  ln_Fy = matrix(log(Fmort[-n_years,,1]), ncol = 2),
  ln_fish_selpars = array(log(3), dim = c(2, 2, 2)),
  ln_srv_selpars = array(log(3), dim = c(2, 2, 2))
)

# Mapping
map = list(
  # ln_InitDevs = factor(rep(NA, length(InitDevs[,1]))),
  # ln_RecDevs = factor(rep(NA, length(RecDevs[-50,1]))),
  ln_sigmaRec = factor(NA),
  RecPars = factor(c(1, NA))
  # ln_q_fish = factor(c(NA, NA)),
  # ln_q_srv = factor(c(NA, NA)),
  # ln_Fy = factor(rep(NA, 100)),
  # ln_M = factor(c(NA, NA))
  # ln_srv_selpars = factor(rep(NA, 8))
)

library(TMB)
# setwd("src")
TMB::compile("Sex_Str_EM.cpp")
dyn.unload(dynlib('Sex_Str_EM'))
dyn.load(dynlib('Sex_Str_EM'))

# Make AD Function here
model_fxn = TMB::MakeADFun(data, parameters, map, random = NULL,
                            DLL= "Sex_Str_EM", silent = FALSE,  
                            checkParameterOrder = TRUE, tracepar = TRUE)

Opt = TMBhelper::fit_tmb( obj = model_fxn,
                          newtonsteps = 3,
                          bias.correct = FALSE,
                          getsd = TRUE,
                          savedir = paste0(getwd(),"/") )

Report = model_fxn$report(model_fxn$env$last.par.best)
sum(Report$pred_srv_len_comps[1,,,1])

ParHat = model_fxn$env$parList()

exp(Opt$SD$par.fixed[names(Opt$SD$par.fixed) == "ln_q_fish"])
exp(Opt$SD$par.fixed[names(Opt$SD$par.fixed) == "ln_q_srv"])
exp(Opt$SD$par.fixed[names(Opt$SD$par.fixed) == "ln_M"])
plot(exp(Opt$SD$par.fixed[names(Opt$SD$par.fixed) == "ln_RecDevs"]), type = "l")
lines(RecDevs[,1], col = "red")

plot(Report$NAA[3,,1], type = "l", col = "black")
# lines(Report$NAA[2,,1], type = "l", col = "black")
lines(NAA[3,,1,1])

plot(SSB[-n_years,1])
lines(Report$SSB, type = "l", col = "blue")

par(mfrow = c(2,2))
plot(Report$NAA[1,,1], type = "l", col = "red")
lines(Report$NAA[1,,2], type = "l", col = "blue")
plot(Report$pred_fish_age_comps[1,,1,1], type = "l", col = "red")
lines(Report$pred_fish_age_comps[1,,2,1], type = "l", col = "blue")
plot(Report$pred_srv_age_comps[1,,1,1], type = "l", col = "red")
lines(Report$pred_srv_age_comps[1,,2,1], type = "l", col = "blue")
plot(Report$pred_fish_len_comps[1,,1,1], type = "l", col = "red")
lines(Report$pred_fish_len_comps[1,,2,1], type = "l", col = "blue")
plot(Report$pred_srv_len_comps[1,,1,1], type = "l", col = "red")
lines(Report$pred_srv_len_comps[1,,2,1], type = "l", col = "blue")
dev.off()

Report$pred_fish_age_comps[1,,1,1]
sum(Report$pred_fish_age_comps[1,,,1])

# Catch
plot(Report$pred_catch[,2])
lines(data$obs_catch[,2])

plot(Report$Fish_Slx[1,,1,1], type = "l", col = "blue")
lines(FishAge_Selex[,1,1])
plot(Report$Fish_Slx[1,,2,2], type = "l", col = "blue")
lines(FishAge_Selex[,2,2])

plot(Report$Srv_Slx[1,,1,2], type = "l", col = "blue")
lines(SrvAge_Selex[,1,2])

plot(Report$pred_fish_index[,1])
lines(Fish_Index[-n_years,1,1])
plot(Report$pred_srv_index[,1])
lines(Srv_Index[-n_years,1,1])

# Likelihoods
sum(Report$srv_age_comp_nLL)
sum(Report$srv_index_nLL)
sum(Report$rec_nLL)
sum(Report$fish_index_nLL)
sum(Report$fish_age_comp_nLL)
sum(Report$fish_len_comp_nLL)
sum(Report$srv_age_comp_nLL)
sum(Report$srv_len_comp_nLL)
sum(Report$catch_nLL)
Report$jnLL


