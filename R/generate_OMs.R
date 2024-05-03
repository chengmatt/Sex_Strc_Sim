# Purpose: To generate OMs for sex-structure simulations (experiment 1)
# Creator: Matthew LH. Cheng
# Date: 10/18/23


# Set up ------------------------------------------------------------------
library(here)
library(tidyverse)
library(readxl)

# Load in all functions from the functions folder
fxn_path <- here("R", "functions")
files <- list.files(fxn_path)
for(i in 1:length(files)) source(here(fxn_path, files[i]))
dir.create(here("output", "Experiment 1"))
dir.create(here("output", "Experiment 2"))
dir.create(here("output", "Experiment 3"))

# Experiment 1 ------------------------------------------------------------

# Read in OMs for experiment 1 (Varying Neff and looking at impact of simulating proportions within vs. across)
oms_exp1 <- read_xlsx(here("input", "generate_OMs.xlsx"), sheet = "OM_Exp1")

for(i in 1:nrow(oms_exp1)) { # make sure to change sims to 750
  
  # Create file directory to save model outputs
  om_path = here("output", "Experiment 1", oms_exp1$OM_Name[i])
  dir.create(om_path)
  
  # Get specified effective sample size (constant for survey and fishery)
  Neff = oms_exp1$Neff[i] 
  comp_across_sex = oms_exp1$comp_across_sex[i] # composition parameterization
  
  # simulate data
  oms = simulate_data(spreadsheet_path = here("input", "Sablefish_Inputs.xlsx"),
                      Fish_Neff_Age = Neff, # Neff * 2
                      Fish_Neff_Len = Neff,
                      Srv_Neff_Age = Neff,
                      Srv_Neff_Len = Neff,
                      F_pattern = "Contrast",
                      comp_across_sex = comp_across_sex,
                      selex_type = "length", 
                      q_Fish = 0.025,
                      cv_Fish_Index = 0.25,
                      q_Srv = 0.025,
                      cv_Srv_Index = 0.25, 
                      sexRatio = c(0.5, 0.5), # keeping sex ratio at 50:50 here
                      growth_control = "chg_males_rel_females",
                      natmort_control = "chg_males_rel_females",
                      growth_control_fct = oms_exp1$growth_control_fct[i], # holding constant at 17.5% difference 
                      natmort_control_fct = oms_exp1$natmort_control_fct[i], # holding constant at 17.5% difference 
                      force_grwth_same_yng = TRUE)  # force minimum age to be similar
  
  # Save as RData file - ifelse for sensitivity tests
  save(oms, file = here(om_path, paste(oms_exp1$OM_Name[i], ".RData", sep = "")))
  plot_OMs(oms, path = om_path)
  
  # plot(oms$SrvAge_Selex[,1,1], ylim = c(0,1))
  # lines(oms$SrvAge_Selex[,2,1])
  # lines(oms$FishAge_Selex[,1,1])
  # lines(oms$FishAge_Selex[,2,1])
  
} # end i loop


# Experiment 2 ------------------------------------------------------------

# Read in OMs for experiment 2 (ignoring sex-structure)
oms_exp2 <- read_xlsx(here("input", "generate_OMs.xlsx"), sheet = "OM_Exp2")

for(i in 1:nrow(oms_exp2)) {
  
  # Create file directory to save model outputs
  om_path = here("output", "Experiment 2", oms_exp2$OM_Name[i])
  dir.create(om_path)
  
  # simulate data
  oms = simulate_data(spreadsheet_path = here("input", "Sablefish_Inputs.xlsx"),
                      Fish_Neff_Age = 75, # 75*2
                      Fish_Neff_Len = 75,
                      Srv_Neff_Age = 75,
                      Srv_Neff_Len = 75,
                      F_pattern = "Contrast",
                      comp_across_sex = "within",
                      selex_type = "length",
                      q_Fish = 0.025,
                      cv_Fish_Index = 0.25,
                      q_Srv = 0.025,
                      cv_Srv_Index = 0.25, 
                      sexRatio = c(0.5, 0.5), # keeping sex ratio at 50:50 here
                      growth_control = "chg_males_rel_females",
                      natmort_control = "chg_males_rel_females",
                      growth_control_fct = oms_exp2$growth_control_fct[i], 
                      natmort_control_fct = oms_exp2$natmort_control_fct[i],  
                      force_grwth_same_yng = TRUE) # force minimum age to be similar
  
  # Save as RData file - ifelse for sensitivity tests
  save(oms, file = here(om_path, paste(oms_exp2$OM_Name[i], ".RData", sep = "")))
  plot_OMs(oms, path = om_path)
  
} # end i loop


# Experiment 3 ------------------------------------------------------------

# Read in OMs for experiment 3 (sex ratio misspecification)
oms_exp3 <- read_xlsx(here("input", "generate_OMs.xlsx"), sheet = "OM_Exp3")

for(i in 1:nrow(oms_exp3)) {
  
  # Create file directory to save model outputs
  om_path = here("output", "Experiment 3", oms_exp3$OM_Name[i])
  dir.create(om_path)
  
  # get sex-ratio here
  sr = as.numeric(strsplit(oms_exp3$Sex_Ratios[i], ",")[[1]]) # females then males
  
  # simulate data
  oms = simulate_data(spreadsheet_path = here("input", "Sablefish_Inputs.xlsx"),
                      Fish_Neff_Age = 75, # z * 2
                      Fish_Neff_Len = 75,
                      Srv_Neff_Age = 75,
                      Srv_Neff_Len = 75,
                      F_pattern = "Contrast",
                      comp_across_sex = "within",
                      selex_type = "length",
                      q_Fish = 0.025,
                      cv_Fish_Index = 0.25,
                      q_Srv = 0.025,
                      cv_Srv_Index = 0.25, 
                      sexRatio = sr, 
                      growth_control = "chg_males_rel_females",
                      natmort_control = "chg_males_rel_females",
                      growth_control_fct = oms_exp3$growth_control_fct[i], # holding constant at 17.5% difference 
                      natmort_control_fct = oms_exp3$natmort_control_fct[i], # holding constant at 17.5% difference 
                      force_grwth_same_yng = TRUE)  # force minimum age to be similar
  
  # Save as RData file - ifelse for sensitivity tests
  save(oms, file = here(om_path, paste(oms_exp3$OM_Name[i], ".RData", sep = "")))
  plot_OMs(oms, path = om_path)
  
} # end i loop
