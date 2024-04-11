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
coverage_all_list = list()

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
  coverage_em_list = list()
  
  for(n_em in 1:length(em_folders)) {
    em_path = here(om_folder, em_folders[n_em]) # list out ems
    load(here(em_path, paste(em_folders[n_em], ".RData", sep = ""))) # load in EMs
    selex_df = data.table::fread(here(em_path, "Selectivity.csv")) # read in selectivity
    growth_df =  data.table::fread(here(em_path, "Growth.csv")) # read in growth
    param_df =  data.table::fread(here(em_path, "Parameters.csv")) # read in parameters
    ts_df =  data.table::fread(here(em_path, "Time_Series.csv")) # read in time series
    sr_df = data.table::fread(here(em_path, "NAA_SexRatios.csv")) # read in sex ratio stuff
    conv_df = data.table::fread(here(em_path, "Convergence.csv")) # read in sex ratio stuff
    coverage_df =  data.table::fread(here(em_path, "Coverage.csv")) # read in time series
    
    for(k in 1:length(model_list)) { # loop through model list
      if(sum(summary(model_list[[k]]$sd_rep)[,2] >= 100, na.rm = TRUE)) sdNA = TRUE else sdNA = FALSE
      conv_df$sdNA[k] = sdNA
    } # get se bounds into convergence
    
    # input into our list
    selex_em_list[[n_em]] = selex_df
    growth_em_list[[n_em]] = growth_df
    param_em_list[[n_em]] = param_df
    ts_em_list[[n_em]] = ts_df
    conv_em_list[[n_em]] = conv_df
    coverage_em_list[[n_em]] = coverage_df
    
  } # end n_em loop
  
  # Output into our "all om" list
  selex_em_df = data.table::rbindlist(selex_em_list)
  growth_em_df = data.table::rbindlist(growth_em_list)
  param_em_df = data.table::rbindlist(param_em_list)
  ts_em_df = data.table::rbindlist(ts_em_list)
  conv_em_df = data.table::rbindlist(conv_em_list)
  naa_em_df = data.table::rbindlist(naa_em_list)
  coverage_em_df = data.table::rbindlist(coverage_em_list)
  
  # Input dataframes into our list
  selex_all_list[[i]] = selex_em_df
  growth_all_list[[i]] = growth_em_df
  param_all_list[[i]] = param_em_df
  ts_all_list[[i]] = ts_em_df
  conv_all_list[[i]] = conv_em_df
  naa_all_list[[i]] = naa_em_df
  coverage_all_list[[i]] = coverage_em_df
  
  print(i)
} # end i loop

# Now output these into our environment as csvs
data.table::fwrite(data.table::rbindlist(selex_all_list), file = here("output", "Experiment_1_Selex.csv"))
data.table::fwrite(data.table::rbindlist(growth_all_list), file = here("output", "Experiment_1_Growth.csv"))
data.table::fwrite(data.table::rbindlist(param_all_list), file = here("output", "Experiment_1_Param.csv"))
data.table::fwrite(data.table::rbindlist(ts_all_list), file = here("output", "Experiment_1_TimeSeries.csv"))
data.table::fwrite(data.table::rbindlist(conv_all_list), file = here("output", "Experiment_1_Convergence.csv"))
data.table::fwrite(data.table::rbindlist(coverage_all_list), file = here("output", "Experiment_1_Coverage.csv"))


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
coverage_all_list = list()

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
  coverage_em_list = list()
  
  for(n_em in 1:length(em_folders)) {
    em_path = here(om_folder, em_folders[n_em]) # list out ems
    load(here(em_path, paste(em_folders[n_em], ".RData", sep = ""))) # load in EMs
    selex_df = data.table::fread(here(em_path, "Selectivity.csv")) # read in selectivity
    growth_df =  data.table::fread(here(em_path, "Growth.csv")) # read in growth
    param_df =  data.table::fread(here(em_path, "Parameters.csv")) # read in parameters
    ts_df =  data.table::fread(here(em_path, "Time_Series.csv")) # read in time series
    sr_df = data.table::fread(here(em_path, "NAA_SexRatios.csv")) # read in sex ratio stuff
    conv_df = data.table::fread(here(em_path, "Convergence.csv")) # read in sex ratio stuff
    coverage_df =  data.table::fread(here(em_path, "Coverage.csv")) # read in time series
    
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
    coverage_em_list[[n_em]] = coverage_df
    
  } # end n_em loop
  
  # Output into our "all om" list
  selex_em_df = data.table::rbindlist(selex_em_list)
  growth_em_df = data.table::rbindlist(growth_em_list)
  param_em_df = data.table::rbindlist(param_em_list)
  ts_em_df = data.table::rbindlist(ts_em_list)
  sr_em_df = data.table::rbindlist(sr_em_list)
  conv_em_df = data.table::rbindlist(conv_em_list)
  coverage_em_df = data.table::rbindlist(coverage_em_list)
  
  # Input dataframes into our list
  selex_all_list[[i]] = selex_em_df
  growth_all_list[[i]] = growth_em_df
  param_all_list[[i]] = param_em_df
  ts_all_list[[i]] = ts_em_df
  sr_all_list[[i]] = sr_em_df
  conv_all_list[[i]] = conv_em_df
  coverage_all_list[[i]] = coverage_em_df
  
} # end i loop

# Now output these into our environment as csvs
data.table::fwrite(data.table::rbindlist(selex_all_list), file = here("output", "Experiment_2_Selex.csv"))
data.table::fwrite(data.table::rbindlist(growth_all_list), file = here("output", "Experiment_2_Growth.csv"))
data.table::fwrite(data.table::rbindlist(param_all_list), file = here("output", "Experiment_2_Param.csv"))
data.table::fwrite(data.table::rbindlist(ts_all_list), file = here("output", "Experiment_2_TimeSeries.csv"))
data.table::fwrite(data.table::rbindlist(sr_all_list), file = here("output", "Experiment_2_SexRatio.csv"))
data.table::fwrite(data.table::rbindlist(conv_all_list), file = here("output", "Experiment_2_Convergence.csv"))
data.table::fwrite(data.table::rbindlist(coverage_all_list), file = here("output", "Experiment_2_Coverage.csv"))


# Experiment 2 (Numbers at Age and SSB CV) ------------------------------------

### Numbers at Age ----------------------------------------------------------
om_list <- list()
for(i in 1:length(exp2_oms)) {
  
  # Go into a given OM folder
  om_folder = here(exp2_path, exp2_oms[i])
  em_folders = list.files(om_folder) # list out em folders
  em_folders = em_folders[!str_detect(em_folders, "RData|pdf")] # remove these
  load(here(om_folder, paste(exp2_oms[i], ".RData", sep = ""))) # load in OMs
  em_list <- list()
  
  for(n_em in 1:length(em_folders)) {
    
    em_path <- here(om_folder, em_folders[n_em]) # list out ems
    load(here(em_path, paste(em_folders[n_em], ".RData", sep = ""))) # load in EMs
    
    # if age-structured only models, then append dataframe for females and calculate it assuming 0.5
    if(em_folders[n_em] %in% c("Age (AgeSel)", "Age (LenSel)")) {
      # Get NAA from simulations
      naa_em <- lapply(model_list, function(x) reshape2::melt(x$rep$NAA)) # EM NAA
      naa_em_df <- data.table::rbindlist(naa_em, idcol = "Var4") # Turn estimated NAA into dataframes 
      naa_em_df$value <- naa_em_df$value * 0.5 # apportion half of these to females
      naa_em_df <- rbind(naa_em_df, naa_em_df %>% mutate(Var3 = 2)) # now rbind! 
    } else{ # if these are not age-structured models
      # Get NAA from simulations
      naa_em <- lapply(model_list, function(x) reshape2::melt(x$rep$NAA)) # EM NAA
      naa_em_df <- data.table::rbindlist(naa_em, idcol = "Var4") # Turn estimated NAA into dataframes 
    } # end else
    
    # Residual munging
    naa_om <- reshape2::melt(oms$NAA[-dim(oms$NAA)[1],,,]) # OM NAA (remove last year)
    naa_em_df <- naa_em_df %>% mutate(EM = em_folders[n_em], OM = exp2_oms[i]) %>% rename(Pred = value)
    naa_em_df <- naa_em_df %>% left_join(naa_om %>% rename(Truth = value), by = c('Var4', "Var3", "Var2", "Var1"))
    naa_em_df <- naa_em_df %>% rename(Year = Var1, Age = Var2, Sex = Var3, Sim = Var4) # rename variables
    em_list[[n_em]] <- naa_em_df
  } # end n_em
  
  em_df <- data.table::rbindlist(em_list) # output as dataframe here
  om_list[[i]] <- em_df # input into om list
  print(i)
  
} # end i loop

# Turn OM list into dataframe
om_df <- data.table::rbindlist(om_list)
data.table::fwrite(om_df, file = here("output", "Experiment_2_NAA.csv"))


### CV in SSB ---------------------------------------------------------------
om_list <- list()
for(i in 1:length(exp2_oms)) {
  
  # Go into a given OM folder
  om_folder = here(exp2_path, exp2_oms[i])
  em_folders = list.files(om_folder) # list out em folders
  em_folders = em_folders[!str_detect(em_folders, "RData|pdf")] # remove these
  load(here(om_folder, paste(exp2_oms[i], ".RData", sep = ""))) # load in OMs
  em_list <- list()
  
  for(n_em in 1:length(em_folders)) {
    
    em_path <- here(om_folder, em_folders[n_em]) # list out ems
    load(here(em_path, paste(em_folders[n_em], ".RData", sep = ""))) # load in EMs

    # Get ssb estimates
    ssb_est <- lapply(model_list, function(x) {
      est <- x$sd_rep$value[names(x$sd_rep$value) == "SSB"]
      data.frame(Year = 1:length(est), Estimate = est)}) # EM SSB
    
    # get SE
    ssb_se <- lapply(model_list, function(x) {
      se <- x$sd_rep$sd[names(x$sd_rep$value) == "SSB"] 
      data.frame(Year = 1:length(se), SE = se)
    }) # EM SSB SE
    
    # Output to dataframe
    ssb_df <- data.table::rbindlist(ssb_est, idcol = "Sim")
    se_df <- data.table::rbindlist(ssb_se, idcol = "Sim")
    
    # Residual munging and naming
    ssb_em_df <- ssb_df %>% left_join(se_df, by = c("Year", "Sim")) %>% 
      mutate(EM = em_folders[n_em], OM = exp2_oms[i], CV = SE/Estimate)

    em_list[[n_em]] <- ssb_em_df
  } # end n_em
  
  em_df <- data.table::rbindlist(em_list) # output as dataframe here
  om_list[[i]] <- em_df # input into om list
  print(i)
  
} # end i loop

# Turn OM list into dataframe
om_df <- data.table::rbindlist(om_list)
data.table::fwrite(om_df, file = here("output", "Experiment_2_SSBCV.csv"))

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
coverage_all_list = list()

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
  coverage_em_list = list()
  
  for(n_em in 1:length(em_folders)) {
    em_path = here(om_folder, em_folders[n_em]) # list out ems
    load(here(em_path, paste(em_folders[n_em], ".RData", sep = ""))) # load in EMs
    selex_df = data.table::fread(here(em_path, "Selectivity.csv")) # read in selectivity
    growth_df =  data.table::fread(here(em_path, "Growth.csv")) # read in growth
    param_df =  data.table::fread(here(em_path, "Parameters.csv")) # read in parameters
    ts_df =  data.table::fread(here(em_path, "Time_Series.csv")) # read in time series
    sr_df = data.table::fread(here(em_path, "NAA_SexRatios.csv")) # read in sex ratio stuff
    conv_df = data.table::fread(here(em_path, "Convergence.csv")) # read in sex ratio stuff
    coverage_df = data.table::fread(here(em_path, "Coverage.csv")) # read in coverage stuff
    
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
    coverage_em_list[[n_em]] = coverage_df
    
  } # end n_em loop
  
  # Output into our "all om" list
  selex_em_df = data.table::rbindlist(selex_em_list)
  growth_em_df = data.table::rbindlist(growth_em_list)
  param_em_df = data.table::rbindlist(param_em_list)
  ts_em_df = data.table::rbindlist(ts_em_list)
  sr_em_df = data.table::rbindlist(sr_em_list)
  conv_em_df = data.table::rbindlist(conv_em_list)
  coverage_em_df = data.table::rbindlist(coverage_em_list)
  
  # Input dataframes into our list
  selex_all_list[[i]] = selex_em_df
  growth_all_list[[i]] = growth_em_df
  param_all_list[[i]] = param_em_df
  ts_all_list[[i]] = ts_em_df
  sr_all_list[[i]] = sr_em_df
  conv_all_list[[i]] = conv_em_df
  coverage_all_list[[i]] = coverage_em_df
  
} # end i loop

# Now output these into our environment as csvs
data.table::fwrite(data.table::rbindlist(selex_all_list), file = here("output", "Experiment_3_Selex.csv"))
data.table::fwrite(data.table::rbindlist(growth_all_list), file = here("output", "Experiment_3_Growth.csv"))
data.table::fwrite(data.table::rbindlist(param_all_list), file = here("output", "Experiment_3_Param.csv"))
data.table::fwrite(data.table::rbindlist(ts_all_list), file = here("output", "Experiment_3_TimeSeries.csv"))
data.table::fwrite(data.table::rbindlist(sr_all_list), file = here("output", "Experiment_3_SexRatio.csv"))
data.table::fwrite(data.table::rbindlist(conv_all_list), file = here("output", "Experiment_3_Convergence.csv"))
data.table::fwrite(data.table::rbindlist(coverage_all_list), file = here("output", "Experiment_3_Coverage.csv"))

