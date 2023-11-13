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

# Experiment 1 ------------------------------------------------------------

# Read in OMs for experiment 1
oms_exp1 <- read_xlsx(here("input", "generate_OMs.xlsx"), sheet = "OM_Exp1")

for(i in 1:nrow(oms_exp1)) {
  
  # Create file directory to save model outputs
  om_path = here("output", "Experiment 1", oms_exp1$OM_Name[i])
  dir.create(om_path)
  
  # simulate data
  oms = simulate_data(spreadsheet_path = here("input", "Sablefish_Inputs.xlsx"),
                Fish_Neff_Age = 100,
                Fish_Neff_Len = 100,
                Srv_Neff_Age = 100,
                Srv_Neff_Len = 100,
                F_pattern = "Contrast",
                comp_across_sex = "across",
                selex_type = "length",
                q_Fish = 0.025,
                cv_Fish_Index = 0.25,
                q_Srv = 0.05,
                cv_Srv_Index = 0.25, 
                sexRatio = c(0.5, 0.5), # keeping sex ratio at 50:50 here
                growth_control = "chg_males_rel_females",
                natmort_control = "chg_males_rel_females",
                growth_control_fct = oms_exp1$growth_control_fct[i], 
                natmort_control_fct = oms_exp1$natmort_control_fct[i],  
                force_grwth_same_yng = FALSE)
  
  # Save as RData file - ifelse for sensitivity tests
  save(oms, file = here(om_path, paste(oms_exp1$OM_Name[i], ".RData", sep = "")))
  plot_OMs(oms, path = om_path)

} # end i loop

# Experiment 2 ------------------------------------------------------------

# Read in OMs for experiment 2
oms_exp2 <- read_xlsx(here("input", "generate_OMs.xlsx"), sheet = "OM_Exp2")

for(i in 1:nrow(oms_exp2)) {
  
  # Create file directory to save model outputs
  om_path = here("output", "Experiment 2", oms_exp2$OM_Name[i])
  dir.create(om_path)
  
  # Define sex ratios here
  sr = as.numeric(strsplit(oms_exp2$Sex_Ratios[i], ",")[[1]]) # females then males
  
  # simulate data
  oms = simulate_data(spreadsheet_path = here("input", "Sablefish_Inputs.xlsx"),
                      Fish_Neff_Age = 100,
                      Fish_Neff_Len = 100,
                      Srv_Neff_Age = 100,
                      Srv_Neff_Len = 100,
                      F_pattern = "Contrast",
                      comp_across_sex = "across",
                      selex_type = "length",
                      q_Fish = 0.025,
                      cv_Fish_Index = 0.3,
                      q_Srv = 0.05,
                      cv_Srv_Index = 0.3, 
                      sexRatio = sr, # keeping sex ratio at 50:50 here
                      growth_control = "chg_males_rel_females",
                      natmort_control = "chg_males_rel_females",
                      growth_control_fct = oms_exp2$growth_control_fct[i], # holding constant at 15% difference 
                      natmort_control_fct = oms_exp2$natmort_control_fct[i], # holding constant at 15% difference 
                      force_grwth_same_yng = FALSE) 
  
  # Save as RData file - ifelse for sensitivity tests
  save(oms, file = here(om_path, paste(oms_exp2$OM_Name[i], ".RData", sep = "")))
  plot_OMs(oms, path = om_path)
  
} # end i loop
