# Purpose: To collate EM outputs
# Creator: Matthew LH. Cheng
# Date: UAF-CFOS


# Set up ------------------------------------------------------------------

library(here)
library(tidyverse)

# Experiment 1 ------------------------------------------------------------

exp1_path = here("output", "Experiment 1")
exp1_oms = list.files(exp1_path)

# storage containers
selex_all_list = list()
growth_all_list = list()
param_all_list = list()
ts_all_list = list()
conv_all_list = list()

for(i in 1:length(exp1_oms)) {
  
  # Go into a given OM folder
  om_folder = here(exp1_path, exp1_oms[i])
  em_folders = list.files(om_folder) # list out em folders
  em_folders = em_folders[!str_detect(em_folders, "RData|pdf")] # remove these
  
  # storage containers for ems - reset for every om
  selex_em_list = list()
  growth_em_list = list()
  param_em_list = list()
  ts_em_list = list()
  conv_em_list = list()
  
  for(n_em in 1:length(em_folders)) {
    em_path = here(om_folder, em_folders[n_em]) # list out ems
    selex_df = data.table::fread(here(em_path, "Selectivity.csv")) # read in selectivity
    growth_df =  data.table::fread(here(em_path, "Growth.csv")) # read in growth
    param_df =  data.table::fread(here(em_path, "Parameters.csv")) # read in parameters
    ts_df =  data.table::fread(here(em_path, "Time_Series.csv")) # read in time series
    conv_df =  data.table::fread(here(em_path, "Convergence.csv")) # read in time series
    
    # input into our list
    selex_em_list[[n_em]] = selex_df
    growth_em_list[[n_em]] = growth_df
    param_em_list[[n_em]] = param_df
    ts_em_list[[n_em]] = ts_df
    conv_em_list[[n_em]] = conv_df
    
  } # end n_em loop
  
  # Output into our "all om" list
  selex_em_df = data.table::rbindlist(selex_em_list)
  growth_em_df = data.table::rbindlist(growth_em_list)
  param_em_df = data.table::rbindlist(param_em_list)
  ts_em_df = data.table::rbindlist(ts_em_list)
  conv_em_df = data.table::rbindlist(conv_em_list)
  
  # Input dataframes into our list
  selex_all_list[[i]] = selex_em_df
  growth_all_list[[i]] = growth_em_df
  param_all_list[[i]] = param_em_df
  ts_all_list[[i]] = ts_em_df
  conv_all_list[[i]] = conv_em_df

} # end i loop

# Now output these into our environment as csvs
data.table::fwrite(data.table::rbindlist(selex_all_list), file = here("output", "Experiment_1_Selex.csv"))
data.table::fwrite(data.table::rbindlist(growth_all_list), file = here("output", "Experiment_1_Growth.csv"))
data.table::fwrite(data.table::rbindlist(param_all_list), file = here("output", "Experiment_1_Param.csv"))
data.table::fwrite(data.table::rbindlist(ts_all_list), file = here("output", "Experiment_1_TimeSeries.csv"))
data.table::fwrite(data.table::rbindlist(conv_all_list), file = here("output", "Experiment_1_Convergence.csv"))


# Experiment 2 ------------------------------------------------------------

exp2_path = here("output", "Experiment 2")
exp2_oms = list.files(exp2_path)

# storage containers
selex_all_list = list()
growth_all_list = list()
param_all_list = list()
ts_all_list = list()
sr_all_list = list()
conv_all_list = list()

for(i in 1:length(exp2_oms)) {
  
  # Go into a given OM folder
  om_folder = here(exp2_path, exp2_oms[i])
  em_folders = list.files(om_folder) # list out em folders
  em_folders = em_folders[!str_detect(em_folders, "RData|pdf")] # remove these
  
  # storage containers for ems - reset for every om
  selex_em_list = list()
  growth_em_list = list()
  param_em_list = list()
  ts_em_list = list()
  sr_em_list = list()
  conv_em_list = list()
  
  for(n_em in 1:length(em_folders)) {
    em_path = here(om_folder, em_folders[n_em]) # list out ems
    selex_df = data.table::fread(here(em_path, "Selectivity.csv")) # read in selectivity
    growth_df =  data.table::fread(here(em_path, "Growth.csv")) # read in growth
    param_df =  data.table::fread(here(em_path, "Parameters.csv")) # read in parameters
    ts_df =  data.table::fread(here(em_path, "Time_Series.csv")) # read in time series
    sr_df = data.table::fread(here(em_path, "NAA_SexRatios.csv")) # read in sex ratio stuff
    conv_df = data.table::fread(here(em_path, "Convergence.csv")) # read in sex ratio stuff
    
    # input into our list
    selex_em_list[[n_em]] = selex_df
    growth_em_list[[n_em]] = growth_df
    param_em_list[[n_em]] = param_df
    ts_em_list[[n_em]] = ts_df
    sr_em_list[[n_em]] = sr_df
    conv_em_list[[n_em]] = conv_df
    
  } # end n_em loop
  
  # Output into our "all om" list
  selex_em_df = data.table::rbindlist(selex_em_list)
  growth_em_df = data.table::rbindlist(growth_em_list)
  param_em_df = data.table::rbindlist(param_em_list)
  ts_em_df = data.table::rbindlist(ts_em_list)
  sr_em_df = data.table::rbindlist(sr_em_list)
  conv_em_df = data.table::rbindlist(conv_em_list)
  
  # Input dataframes into our list
  selex_all_list[[i]] = selex_em_df
  growth_all_list[[i]] = growth_em_df
  param_all_list[[i]] = param_em_df
  ts_all_list[[i]] = ts_em_df
  sr_all_list[[i]] = sr_em_df
  conv_all_list[[i]] = conv_em_df
  
} # end i loop

# Now output these into our environment as csvs
data.table::fwrite(data.table::rbindlist(selex_all_list), file = here("output", "Experiment_2_Selex.csv"))
data.table::fwrite(data.table::rbindlist(growth_all_list), file = here("output", "Experiment_2_Growth.csv"))
data.table::fwrite(data.table::rbindlist(param_all_list), file = here("output", "Experiment_2_Param.csv"))
data.table::fwrite(data.table::rbindlist(ts_all_list), file = here("output", "Experiment_2_TimeSeries.csv"))
data.table::fwrite(data.table::rbindlist(sr_all_list), file = here("output", "Experiment_2_SexRatio.csv"))
data.table::fwrite(data.table::rbindlist(conv_all_list), file = here("output", "Experiment_2_Convergence.csv"))
