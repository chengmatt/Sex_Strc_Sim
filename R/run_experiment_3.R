# Purpose: To run EMs for experiment 2 as a full factorial (sex ratio differences)
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date: 10/19/23

# Set up ------------------------------------------------------------------
library(here)
library(tidyverse)
library(TMB)
library(readxl)
library(doSNOW)
library(parallel)

# Compile
# setwd("src")
TMB::compile("Sex_Str_EM.cpp")
dyn.unload(dynlib("Sex_Str_EM"))
dyn.load(dynlib('Sex_Str_EM'))

rm(list=ls()) # remove objects prior to running
ncores <- detectCores() 
# Register cluster here
cl <- makeCluster(ncores - 2)
registerDoSNOW(cl)

# Load in all functions into the environment
fxn_path <- here("R", "functions")
# Load in all functions from the functions folder
files <- list.files(fxn_path)
for(f in 1:length(files)) source(here(fxn_path, files[f]))

# Set up paths to OMs
om_path = here("output", "Experiment 3")
om_names = list.files(om_path)

# Read in OMs for experiment 1 and set up
oms_exp3 <- read_xlsx(here("input", "generate_OMs.xlsx"), sheet = "OM_Exp3") 

# Specify CV
catch_cv <- 0.025
catch_sd <- sqrt(log(catch_cv^2 + 1))

# Run Experiment 3 --------------------------------------------------------

for(n_om in 1:nrow(oms_exp3)) {
  
  # Specify the OM scenario folder
  om_scenario = here(om_path, oms_exp3$OM_Name[n_om])
  om_name = oms_exp3$OM_Name[n_om] # om name
  load(here(om_scenario, paste(om_name,".RData",sep = "")))
  list2env(oms,globalenv()) # output into global environment
  
  # Add observation error to catch common to a run
  oms$Total_Catch <- oms$Total_Catch * exp(rnorm(prod(dim(oms$Total_Catch)), -catch_sd^2/2, catch_sd ))
  oms$Total_Catch_Sex <- oms$Total_Catch_Sex * exp(rnorm(prod(dim(oms$Total_Catch_Sex)), -catch_sd^2/2, catch_sd ))
  
  # read in EM experiments
  ems_exp3 <- read_xlsx(here("input", "run_EMs.xlsx"), sheet = "EM_Exp3", na = "NA")

  # If these are variants we are testing to understand model behavior
  if(str_detect(om_name, "No") == TRUE) ems_exp3 = ems_exp3 %>% filter(str_detect(EM_Name, "Fix"))
  
  for(n_em in 1:nrow(ems_exp3)) {
    
    # Set up specifications here -------------------------------------------
    n_sexes_em = ems_exp3$n_sexes[n_em] # number of sexes to model
    sex_specific_em = ems_exp3$sex_specific[n_em] # whether this is a sex-specific assessment
    share_M_sex_em = ems_exp3$share_M_sex[n_em] # wether to share M or not
    use_fish_sexRatio_em = ems_exp3$use_fish_sexRatio[n_em] # whether to use sex ratio nLL
    use_srv_sexRatio_em = ems_exp3$use_srv_sexRatio[n_em] # whether to use sex ratio nLL
    fish_age_prop_em = ems_exp3$fish_age_prop[n_em] # whether to do proportions across or within
    srv_age_prop_em = ems_exp3$srv_age_prop[n_em] # whether to do proportions across or within
    fish_len_prop_em = ems_exp3$fish_len_prop[n_em] # whether to do proportions across or within
    srv_len_prop_em = ems_exp3$srv_len_prop[n_em] # whether to do proportions across or within
    est_sexRatio_par_em = ems_exp3$est_sexRatio_par[n_em] # whether to estimate sex Ratios
    sexRatio_al_or_y_em = ems_exp3$sexRatio_al_or_y[n_em] # if we want to fit sex ratio as within year only or both
    em_name = ems_exp3$EM_Name[n_em] # em name

    # Run Simulations here ----------------------------------------------------
    sim_models <- foreach(sim = 1:n_sims, .packages = c("TMB", "here", "tidyverse")) %dopar% {
      
      TMB::compile("Sex_Str_EM.cpp")
      dyn.load(dynlib('Sex_Str_EM'))
      
      # estimate biological weight at age
      biologicals = get_biologicals(n_sexes, n_ages, age_bins, len_bins, Srv_LAA, Srv_LW, sim = sim)
      
      # If single sex model
      if(n_sexes_em == 1) {
        waa_em = biologicals$waa_nosex
        al_matrix_em = biologicals$al_matrix_sexagg
        sexRatio_em = c(1, 1)
      } # end if
      
      # If multi-sex model
      if(n_sexes_em == 2) {
        waa_em = biologicals$waa_sex
        al_matrix_em = biologicals$al_matrix_sexsp
        sexRatio_em = c(0.5, 0.5)
      } # end if  
      
      # Prepare EM inputs into assessment
      em_inputs = prepare_EM_inputs(sim = sim,
                                    # EM_Parameterization
                                    n_sexes = n_sexes_em,
                                    sex_specific = sex_specific_em, 
                                    share_M_sex = share_M_sex_em,
                                    use_fish_sexRatio = use_fish_sexRatio_em,
                                    use_srv_sexRatio = use_srv_sexRatio_em,
                                    sexRatio_al_or_y = sexRatio_al_or_y_em,
                                    
                                    # Proportion treatment
                                    fish_age_prop = fish_age_prop_em,
                                    srv_age_prop = srv_age_prop_em,
                                    fish_len_prop = fish_len_prop_em,
                                    srv_len_prop = srv_len_prop_em,
                                    
                                    # Biologicals
                                    WAA = waa_em,
                                    age_len_transition = al_matrix_em,
                                    
                                    # Fixed controls
                                    selex_type = "length",
                                    est_sexRatio_par = est_sexRatio_par_em,
                                    sexRatio = sexRatio_em,
                                    fit_sexsp_catch = FALSE, 
                                    # Aggregating comps
                                    agg_fish_age = FALSE,
                                    agg_srv_age = FALSE, 
                                    agg_fish_len = FALSE,
                                    agg_srv_len = FALSE,
                                    use_fish_index = FALSE,
                                    # Parameter fixing
                                    fix_pars = c("h", "ln_sigmaRec", "ln_q_fish"))
      
      # run model here
      model = run_model(data = em_inputs$data, 
                        parameters = em_inputs$parameters, 
                        map = em_inputs$map, silent = TRUE, n.newton = 3)
      
      # extract quantities
      quants_df = get_quantities(biologicals = biologicals,
                                 model = model, sim = sim, om_name = om_name,
                                 em_name = em_name, n_sexes_em = n_sexes_em)
      
      # Output this into a list when we're done
      all_obj_list = list(model, quants_df$ts_df, quants_df$NAA_sr_female_df,
                          quants_df$grwth_df, quants_df$selex_all_df,
                          quants_df$conv_df, quants_df$par_df, quants_df$coverage_df) 
      
    } # end foreach loop
    
    # process outputs here
    model_list = list()
    params = data.frame()
    time_series = data.frame()
    NAA_sexratio = data.frame()
    growth_df = data.frame()
    selex_df = data.frame()
    convergence_df = data.frame()
    coverage = data.frame()
    
    # Save files and output
    for(s in 1:n_sims) {
      model_list[[s]] = sim_models[[s]][[1]]
      time_series = rbind(time_series, sim_models[[s]][[2]])
      NAA_sexratio = rbind(NAA_sexratio, sim_models[[s]][[3]])
      growth_df = rbind(growth_df, sim_models[[s]][[4]])
      selex_df = rbind(selex_df, sim_models[[s]][[5]])
      convergence_df = rbind(convergence_df, sim_models[[s]][[6]])
      params = rbind(params, sim_models[[s]][[7]])
      coverage = rbind(coverage, sim_models[[s]][[8]])
    } # end s loop
    
    # Now, save our results - create directory to store results first
    em_path_res <- here(om_path, om_name, em_name)
    dir.create(em_path_res)
    save(model_list, file = here(em_path_res, paste(em_name, ".RData", sep = ""))) # save models
    
    # save dataframes
    write.csv(time_series, here(em_path_res, 'Time_Series.csv'), row.names = FALSE)
    write.csv(NAA_sexratio, here(em_path_res, 'NAA_SexRatios.csv'), row.names = FALSE)
    write.csv(growth_df, here(em_path_res, 'Growth.csv'), row.names = FALSE)
    write.csv(selex_df %>% drop_na(), here(em_path_res, 'Selectivity.csv'), row.names = FALSE)
    write.csv(convergence_df, here(em_path_res, 'Convergence.csv'), row.names = FALSE)
    write.csv(params, here(em_path_res, 'Parameters.csv'), row.names = FALSE)
    write.csv(coverage, here(em_path_res, 'Coverage.csv'), row.names = FALSE)
    
    plot_EMs(time_series = time_series, NAA_sexratio = NAA_sexratio,
             growth_df = growth_df, selex_df = selex_df %>% drop_na(), params = params,
             path = em_path_res)
    
    # Progress
    cat(crayon::yellow("EM", n_em, "out of", nrow(ems_exp3)))
    
  } # end em loop
  
  # Progress
  cat(crayon::yellow("OM", n_om, "out of", nrow(oms_exp3)))
} # end om loop
