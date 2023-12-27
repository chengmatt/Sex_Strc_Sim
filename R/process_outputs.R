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
naa_all_list = list()

for(i in 1:length(exp1_oms)) {
  
  # Go into a given OM folder
  om_folder = here(exp1_path, exp1_oms[i])
  em_folders = list.files(om_folder) # list out em folders
  em_folders = em_folders[!str_detect(em_folders, "RData|pdf")] # remove these
  load(here(om_folder, paste(exp1_oms[i], ".RData", sep = ""))) # load in OMs
  
  # Read in OM NAA dataframe
  om_NAA = reshape2::melt(oms$NAA) 
  names(om_NAA) = c("Year", "Age", "Sex", "sim", "Truth") # rename df
  om_NAA = om_NAA %>% filter(Year != max(Year)) # remove last year
  om_NAA$OM = exp1_oms[i] # denote OM
  
  # storage containers for ems - reset for every om
  selex_em_list = list()
  growth_em_list = list()
  param_em_list = list()
  ts_em_list = list()
  conv_em_list = list()
  naa_em_list = list()
  
  for(n_em in 1:length(em_folders)) {
    em_path = here(om_folder, em_folders[n_em]) # list out ems
    selex_df = data.table::fread(here(em_path, "Selectivity.csv")) # read in selectivity
    growth_df =  data.table::fread(here(em_path, "Growth.csv")) # read in growth
    param_df =  data.table::fread(here(em_path, "Parameters.csv")) # read in parameters
    ts_df =  data.table::fread(here(em_path, "Time_Series.csv")) # read in time series
    conv_df =  data.table::fread(here(em_path, "Convergence.csv")) # read in time series
    load(here(em_path, paste(em_folders[n_em], ".RData", sep = ""))) # load in EMs
    em_NAA_store = data.frame() # storage dataframe
    
    for(k in 1:length(model_list)) { # loop through model list
      em_NAA = model_list[[k]]$rep$NAA # read in numbers at age
      em_NAA_df = reshape2::melt(em_NAA) # munge into dataframe
      names(em_NAA_df) = c("Year", "Age", "Sex", "Pred") # rename columns
      em_NAA_df$sim = k  # denote simulation number
      
      # if this is an age-structured assessment only
      if(em_folders[n_em] == "Age") {
        em_NAA_df$Pred = em_NAA_df$Pred * 0.5 # Duplicating dataframe and multiply by 0.5 because that's how sexes are denoted
        em_NAA_df_m = em_NAA_df # denote a new dataframe for males
        em_NAA_df_m$Sex = 2 # denote males
        em_NAA_df = rbind(em_NAA_df, em_NAA_df_m) # rbind to males 
      } # if age structured assessment only 
      
      em_NAA_df$EM = em_folders[n_em] # denote EM
      
      # Now left_join to OM dataframe
      em_NAA_df = em_NAA_df %>% left_join(om_NAA, by = c("Year", "Age", "Sex", "sim"))
      em_NAA_df = em_NAA_df %>% left_join(conv_df %>% select(OM, EM, sim, convergence), by = c("OM", "EM", "sim"))
      em_NAA_store = rbind(em_NAA_df, em_NAA_store) # rbind to store
      if(sum(summary(model_list[[k]]$sd_rep)[,2] >= 100, na.rm = TRUE)) sdNA = TRUE else sdNA = FALSE
      conv_df$sdNA[k] = sdNA
    } # end i
    
    # input into our list
    selex_em_list[[n_em]] = selex_df
    growth_em_list[[n_em]] = growth_df
    param_em_list[[n_em]] = param_df
    ts_em_list[[n_em]] = ts_df
    conv_em_list[[n_em]] = conv_df
    naa_em_list[[n_em]] = em_NAA_store
    
  } # end n_em loop
  
  # Output into our "all om" list
  selex_em_df = data.table::rbindlist(selex_em_list)
  growth_em_df = data.table::rbindlist(growth_em_list)
  param_em_df = data.table::rbindlist(param_em_list)
  ts_em_df = data.table::rbindlist(ts_em_list)
  conv_em_df = data.table::rbindlist(conv_em_list)
  naa_em_df = data.table::rbindlist(naa_em_list)
  
  # Input dataframes into our list
  selex_all_list[[i]] = selex_em_df
  growth_all_list[[i]] = growth_em_df
  param_all_list[[i]] = param_em_df
  ts_all_list[[i]] = ts_em_df
  conv_all_list[[i]] = conv_em_df
  naa_all_list[[i]] = naa_em_df

  print(i)
} # end i loop

# Now output these into our environment as csvs
data.table::fwrite(data.table::rbindlist(selex_all_list), file = here("output", "Experiment_1_Selex.csv"))
data.table::fwrite(data.table::rbindlist(growth_all_list), file = here("output", "Experiment_1_Growth.csv"))
data.table::fwrite(data.table::rbindlist(param_all_list), file = here("output", "Experiment_1_Param.csv"))
data.table::fwrite(data.table::rbindlist(ts_all_list), file = here("output", "Experiment_1_TimeSeries.csv"))
data.table::fwrite(data.table::rbindlist(conv_all_list), file = here("output", "Experiment_1_Convergence.csv"))
data.table::fwrite(data.table::rbindlist(naa_all_list), file = here("output", "Experiment_1_NAA.csv"))


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
  load(here(om_folder, paste(exp2_oms[i], ".RData", sep = ""))) # load in OMs
  
  # storage containers for ems - reset for every om
  selex_em_list = list()
  growth_em_list = list()
  param_em_list = list()
  ts_em_list = list()
  sr_em_list = list()
  conv_em_list = list()
  
  for(n_em in 1:length(em_folders)) {
    em_path = here(om_folder, em_folders[n_em]) # list out ems
    load(here(em_path, paste(em_folders[n_em], ".RData", sep = ""))) # load in EMs
    selex_df = data.table::fread(here(em_path, "Selectivity.csv")) # read in selectivity
    growth_df =  data.table::fread(here(em_path, "Growth.csv")) # read in growth
    param_df =  data.table::fread(here(em_path, "Parameters.csv")) # read in parameters
    ts_df =  data.table::fread(here(em_path, "Time_Series.csv")) # read in time series
    sr_df = data.table::fread(here(em_path, "NAA_SexRatios.csv")) # read in sex ratio stuff
    conv_df = data.table::fread(here(em_path, "Convergence.csv")) # read in sex ratio stuff
    
    for(k in 1:length(model_list)) { # loop through model list
      if(sum(summary(model_list[[k]]$sd_rep)[,2] >= 100, na.rm = TRUE)) sdNA = TRUE else sdNA = FALSE
      conv_df$sdNA[k] = sdNA
    } # get se bounds into convergence
      
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


# Exp 3 -------------------------------------------------------------------

exp3_path = here("output", "Experiment 3")
exp3_oms = list.files(exp3_path)

# storage containers
selex_all_list = list()
growth_all_list = list()
param_all_list = list()
ts_all_list = list()
sr_all_list = list()
conv_all_list = list()

for(i in 1:length(exp3_oms)) {
  
  # Go into a given OM folder
  om_folder = here(exp3_path, exp3_oms[i])
  em_folders = list.files(om_folder) # list out em folders
  em_folders = em_folders[!str_detect(em_folders, "RData|pdf")] # remove these
  load(here(om_folder, paste(exp3_oms[i], ".RData", sep = ""))) # load in OMs
  
  # storage containers for ems - reset for every om
  selex_em_list = list()
  growth_em_list = list()
  param_em_list = list()
  ts_em_list = list()
  sr_em_list = list()
  conv_em_list = list()
  
  for(n_em in 1:length(em_folders)) {
    em_path = here(om_folder, em_folders[n_em]) # list out ems
    load(here(em_path, paste(em_folders[n_em], ".RData", sep = ""))) # load in EMs
    selex_df = data.table::fread(here(em_path, "Selectivity.csv")) # read in selectivity
    growth_df =  data.table::fread(here(em_path, "Growth.csv")) # read in growth
    param_df =  data.table::fread(here(em_path, "Parameters.csv")) # read in parameters
    ts_df =  data.table::fread(here(em_path, "Time_Series.csv")) # read in time series
    sr_df = data.table::fread(here(em_path, "NAA_SexRatios.csv")) # read in sex ratio stuff
    conv_df = data.table::fread(here(em_path, "Convergence.csv")) # read in sex ratio stuff
    
    for(k in 1:length(model_list)) { # loop through model list
      if(sum(summary(model_list[[k]]$sd_rep)[,2] >= 100, na.rm = TRUE)) sdNA = TRUE else sdNA = FALSE
      conv_df$sdNA[k] = sdNA
    } # get se bounds into convergence
    
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
data.table::fwrite(data.table::rbindlist(selex_all_list), file = here("output", "Experiment_3_Selex.csv"))
data.table::fwrite(data.table::rbindlist(growth_all_list), file = here("output", "Experiment_3_Growth.csv"))
data.table::fwrite(data.table::rbindlist(param_all_list), file = here("output", "Experiment_3_Param.csv"))
data.table::fwrite(data.table::rbindlist(ts_all_list), file = here("output", "Experiment_3_TimeSeries.csv"))
data.table::fwrite(data.table::rbindlist(sr_all_list), file = here("output", "Experiment_3_SexRatio.csv"))
data.table::fwrite(data.table::rbindlist(conv_all_list), file = here("output", "Experiment_3_Convergence.csv"))

