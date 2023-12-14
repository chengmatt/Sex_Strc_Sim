# Purpose: To read in input parameters for the OM from an excel spreadsheet. The function calls another
# function to which creates OM objects to store values within it.
# Creator: Matthew LH. Cheng
# Date: 10/30/22

#' @param spreadsheet_path Path to input spreadsheet

read_params <- function(spreadsheet_path) {
  
  require(readxl)
  require(tidyverse)
  require(here)
  
  source(here("R", "functions", "create_OM_objects.R")) # Create OM objects 

# Controls ----------------------------------------------------------------
  ctl <- read_xlsx(spreadsheet_path, sheet = "Controls")
  
  # Read in and save objects in our environment
  n_years <<- ctl$Value[ctl$Par == "n_years"] # Number of years
  n_sims <<- ctl$Value[ctl$Par == "n_sims"] # Number of simulations
  n_sexes <<- ctl$Value[ctl$Par == "n_sexes"] # Numbers of sexes
  n_fish_fleets <<- ctl$Value[ctl$Par == "n_fish_fleets"] # Numbers of fishery fleets
  n_srv_fleets <<- ctl$Value[ctl$Par == "n_srv_fleets"] # Numbers of survey fleets

# Age and length Bins ----------------------------------------------------------------
  bins <- read_xlsx(spreadsheet_path, sheet = "Bins")
  # Read in ages
  age_bins <<- bins$ages[!is.na(bins$ages)]
  n_ages <<- length(age_bins)
  # Read in lens
  len_bins <<- bins$lens[!is.na(bins$lens)]
  n_lens <<- length(len_bins)
  len_mids <<- len_bins[1:(length(len_bins) - 1)] + diff(len_bins) / 2 # Get midpoint of lengths
  
# Maturity at age ---------------------------------------------------------
  maturity <- read_xlsx(spreadsheet_path, sheet = "Maturity_At_Age")
  mat_at_age <<- as.matrix(maturity, ncol = age_bins, nrow = n_sexes)
  
# Growth Parameters -----------------------------------------------------------
  grwth <- read_xlsx(spreadsheet_path, sheet = "Growth_Param") # females then males
  k <<- as.numeric(as.vector(grwth[1,1:2])) # brody coefficient
  L_inf <<- as.numeric(as.vector(grwth[2,1:2])) # Linf parameter
  t0 <<- as.numeric(as.vector(grwth[3,1:2])) # t0 parameter
  beta_wl <<- as.numeric(as.vector(grwth[4,1:2])) # beta parameter
  alpha_wl <<- as.numeric(as.vector(grwth[5,1:2])) # alpha parameter
  vonB_sd1 <<- as.numeric(as.vector(grwth[6,1:2])) # vonbert sd1 min len
  vonB_sd2 <<- as.numeric(as.vector(grwth[7,1:2])) # vonbert sd max len
  wl_sd <<- as.numeric(as.vector(grwth[8,1:2])) # weight length sd
  
# Recruitment + Mortality-------------------------------------------------------------
  recruitment_pars <- read_xlsx(spreadsheet_path, sheet = "Recruitment_Mortality")
  h <<- as.numeric(recruitment_pars$Value[recruitment_pars$Par == "h"]) # Steepness (Recruitment at 20% of SSB0)
  r0 <<- as.numeric(recruitment_pars$Value[recruitment_pars$Par == "r0"] ) # Virgin Recruitment
  sigma_rec <<- as.numeric(recruitment_pars$Value[recruitment_pars$Par == "sigma_rec"]) # Recruitment variability
  M_Female <<- as.numeric(recruitment_pars$Value[recruitment_pars$Par == "M_Female"]) # Female natural mortality
  M_Male <<- as.numeric(recruitment_pars$Value[recruitment_pars$Par == "M_Male"]) # Male natural mortality
  M <<- c(M_Female, M_Male) # combined male and female natural mortality


# Selectivity -------------------------------------------------------------
  selex_pars <- read_xlsx(spreadsheet_path, sheet = "Selex") # length-based selectivity
  fish_len_slope <<- as.numeric(selex_pars$Value[selex_pars$Par == "fish_len_slope"])
  fish_len_midpoint <<- as.numeric(selex_pars$Value[selex_pars$Par == "fish_len_midpoint"])
  srv_len_slope <<- as.numeric(selex_pars$Value[selex_pars$Par == "srv_len_slope"])
  srv_len_midpoint <<- as.numeric(selex_pars$Value[selex_pars$Par == "srv_len_midpoint"])
  fish_age_slope_f <<- as.numeric(selex_pars$Value[selex_pars$Par == "fish_age_slope_f"])
  fish_age_slope_m <<- as.numeric(selex_pars$Value[selex_pars$Par == "fish_age_slope_m"])
  fish_age_midpoint_f <<- as.numeric(selex_pars$Value[selex_pars$Par == "fish_age_midpoint_f"])
  fish_age_midpoint_m <<- as.numeric(selex_pars$Value[selex_pars$Par == "fish_age_midpoint_m"])
  srv_age_slope_f <<- as.numeric(selex_pars$Value[selex_pars$Par == "srv_age_slope_f"])
  srv_age_slope_m <<- as.numeric(selex_pars$Value[selex_pars$Par == "srv_age_slope_m"])
  srv_age_midpoint_f <<- as.numeric(selex_pars$Value[selex_pars$Par == "srv_age_midpoint_f"])
  srv_age_midpoint_m <<- as.numeric(selex_pars$Value[selex_pars$Par == "srv_age_midpoint_m"])
  
  print("### Input parameters have been read in and OM objects have been created ###")
} # end function
