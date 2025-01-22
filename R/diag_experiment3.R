# Purpose: Run model diagnostics for experiment 3 (likelihood profiles and OSA residuals)
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date: 12/12/24

# Set up quick theme
theme_tj = function() {
  theme_bw() +
    theme(legend.position = "none",
          axis.text = element_text(size = 13, color = "black"),
          axis.title = element_text(size = 15, color =  "black"),
          strip.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 13),
          plot.title = element_text(size = 16))
}

# Set up ------------------------------------------------------------------ 
library(compResidual)
library(here)
library(tidyverse)
library(TMB)
library(readxl)
library(doSNOW)
library(parallel)

ncores <- detectCores() 
# Register cluster here
cl <- makeCluster(ncores - 2)
registerDoSNOW(cl)

# setwd("src")
TMB::compile("Sex_Str_EM.cpp")
dyn.load(dynlib('Sex_Str_EM'))

# Load in all functions into the environment
fxn_path <- here("R", "functions")
# Load in all functions from the functions folder
files <- list.files(fxn_path)
for(f in 1:length(files)) source(here(fxn_path, files[f]))

load(here("output", "Experiment 3", "SR (F60, M40)", "SR (F60, M40).RData")) # Load in OMs
load(here("output", "Experiment 3", "SR (F60, M40)", "FixSR", "FixSR.RData")) # Load in EMs

# Composition Residuals --------------------------------------------------------
om_srv_ages <- oms$Fish_AgeComps[-36,,1,1,1] # load in observed ages
pred_srv_ages <- model_list[[1]]$rep$pred_fish_age_comps[,,2,1] # load in predicted ages
resmult_df <- compResidual::resMulti(om_srv_ages, pred_srv_ages) # get residuals

# plot!
pdf(here("figs", "Experiment 3", "OSA_Residuals.pdf"), width = 13)
plot(resmult_df)
dev.off()

# Likelihood Profiles -----------------------------------------------------
# Read in OMs for experiment 1 and set up
oms_exp3 <- read_xlsx(here("input", "generate_OMs.xlsx"), sheet = "OM_Exp3") 
ems_exp3 <- read_xlsx(here("input", "run_EMs.xlsx"), sheet = "EM_Exp3", na = "NA")
conv_df <- read.csv(here("output", "Experiment 3", "SR (F60, M40)", "FixSR", "Convergence.csv")) %>%  
  filter(convergence == "Converged")# read in converged runs only
list2env(oms,globalenv()) # output OMs into global environment

n_em <- 2 # get EM we want to profile
set.seed(666)
sim_num <- sample(conv_df$sim, 10) # sample x simulations to re-run and profile across

### Profile R0 --------------------------------------------------------------
# progress bar
pb <- txtProgressBar(max = length(sim_num), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

r0_nLL_df = data.frame() # set up likelihood profile dataframe
r0_vary = seq(1.3, 6, 0.15) # vary R0
r0_em_list = list()

for(n_em in 1:nrow(ems_exp3)) {
  
  profiles <- foreach(sim = 1:length(sim_num), .packages = c("TMB", "here", "tidyverse"),
          .options.snow = opts) %dopar% {
    # Set up inputs
    use_fish_sexRatio_em = ems_exp3$use_fish_sexRatio[n_em] # whether to use sex ratio nLL
    use_srv_sexRatio_em = ems_exp3$use_srv_sexRatio[n_em] # whether to use sex ratio nLL
    fish_age_prop_em = ems_exp3$fish_age_prop[n_em] # whether to do proportions across or within
    srv_age_prop_em = ems_exp3$srv_age_prop[n_em] # whether to do proportions across or within
    fish_len_prop_em = ems_exp3$fish_len_prop[n_em] # whether to do proportions across or within
    srv_len_prop_em = ems_exp3$srv_len_prop[n_em] # whether to do proportions across or within
    est_sexRatio_par_em = ems_exp3$est_sexRatio_par[n_em] # whether to estimate sex Ratios
    share_M_sex_em = ems_exp3$share_M_sex[n_em] # wether to share M or not
    sexRatio_al_or_y_em = ems_exp3$sexRatio_al_or_y[n_em] # if we want to fit sex ratio as within year only or both
    sexRatio_fixed = c(0.5, 0.5) # fixed sex ratio values to mis-specify initial sex-ratios
    em_name = ems_exp3$EM_Name[n_em] # em name
    
    # estimate biologicals
    biologicals = get_biologicals(n_sexes, n_ages, age_bins, len_bins, Srv_LAA, Srv_LW, sim = sim_num[sim])
    
    # Prepare EM inputs into assessment
    em_inputs = prepare_EM_inputs(sim = sim_num[sim],
                                  # EM_Parameterization
                                  use_fish_sexRatio = use_fish_sexRatio_em,
                                  use_srv_sexRatio = use_srv_sexRatio_em,
                                  est_sexRatio_par = est_sexRatio_par_em,
                                  sexRatio = sexRatio_fixed,
                                  share_M_sex = share_M_sex_em, 
                                  sexRatio_al_or_y = sexRatio_al_or_y_em,
                                  
                                  # Proportion treatment
                                  fish_age_prop = fish_age_prop_em,
                                  srv_age_prop = srv_age_prop_em,
                                  fish_len_prop = fish_len_prop_em,
                                  srv_len_prop = srv_len_prop_em,
                                  
                                  # Fixed controls
                                  n_sexes = 2,
                                  sex_specific = TRUE, 
                                  fit_sexsp_catch = FALSE,
                                  selex_type = "length",
                                  # Aggregating comps
                                  agg_fish_age = FALSE,
                                  agg_srv_age = FALSE, 
                                  agg_fish_len = FALSE,
                                  agg_srv_len = FALSE,
                                  catch_cv = c(0.01), 
                                  use_fish_index = FALSE,
                                  # Biologicals
                                  WAA = biologicals$waa_sex,
                                  age_len_transition = biologicals$al_matrix_sexsp,
                                  # Parameter fixing
                                  fix_pars = c("h", "ln_sigmaRec", "ln_q_fish"))
    
    em_inputs$map$RecPars = factor(c(NA, NA)) # profile R0 
    like_df_all <- data.frame() # empty dataframe for storage
    
    for(i in 1:length(r0_vary)) {
      em_inputs$parameters$RecPars[1] = r0_vary[i] # input varying r0 here
      TMB::compile("Sex_Str_EM.cpp")
      dyn.load(dynlib('Sex_Str_EM'))
      
      # run model here
      model = run_model(data = em_inputs$data, 
                        parameters = em_inputs$parameters, 
                        map = em_inputs$map, silent = TRUE, n.newton = 1)
      # put data into dataframe
      like_df = data.frame(nLL = c(model$rep$jnLL, 
                                   model$rep$rec_nLL,
                                   sum(model$rep$fish_age_comp_nLL[,1,1]),
                                   sum(model$rep$fish_age_comp_nLL[,2,1]),
                                   sum(model$rep$srv_age_comp_nLL[,1,1]),
                                   sum(model$rep$srv_age_comp_nLL[,2,1]),
                                   sum(model$rep$catch_nLL),
                                   sum(model$rep$srv_index_nLL),
                                   sum(model$rep$srv_len_comp_nLL[,1,1]),
                                   sum(model$rep$srv_len_comp_nLL[,2,1]),
                                   sum(model$rep$fish_len_comp_nLL[,1,1]),
                                   sum(model$rep$fish_len_comp_nLL[,1,1]),
                                   sum( model$rep$sexRatio_fishage_nLL[,1,1]),
                                   sum( model$rep$sexRatio_srvage_nLL[,1,1]),
                                   sum( model$rep$sexRatio_fishlen_nLL[,1,1]),
                                   sum( model$rep$sexRatio_srvlen_nLL[,1,1])),
                           comp = c("Joint", "Rec", "FishAgeF", "FishAgeM",
                                    "SrvAgeF", "SrvAgeM", "Catch", "SrvIndx", 
                                    "SrvLenF", "SrvLenM", "FishLenF", "FishLenM",
                                    "SRFishAge", "SRSrvAge", "SRFishLen", "SRSrvLen"),
                           par = r0_vary[i], sim = sim_num[sim], type = "R0", em = em_name)
      
      like_df_all <- rbind(like_df_all, like_df)
    } # end i for likelihood profile loop
    like_df_all
  } # end sim for simulation loop
  r0_em_list[[n_em]] <- profiles # put profiles in here
} # end em loop

r0_profile_df <- bind_rows(r0_em_list) # output list of profiled values

### Profile Female M --------------------------------------------------------------
m_nLL_df = data.frame() # set up likelihood profile dataframe
m_vary = log(seq(0.01, 0.35, 0.01)) # vary M
m_em_list <- list()

for(n_em in 1:nrow(ems_exp3)) {
  
  profiles <- foreach(sim = 1:length(sim_num), .packages = c("TMB", "here", "tidyverse"),
                      .options.snow = opts) %dopar% {    
    # Set up inputs
    use_fish_sexRatio_em = ems_exp3$use_fish_sexRatio[n_em] # whether to use sex ratio nLL
    use_srv_sexRatio_em = ems_exp3$use_srv_sexRatio[n_em] # whether to use sex ratio nLL
    fish_age_prop_em = ems_exp3$fish_age_prop[n_em] # whether to do proportions across or within
    srv_age_prop_em = ems_exp3$srv_age_prop[n_em] # whether to do proportions across or within
    fish_len_prop_em = ems_exp3$fish_len_prop[n_em] # whether to do proportions across or within
    srv_len_prop_em = ems_exp3$srv_len_prop[n_em] # whether to do proportions across or within
    est_sexRatio_par_em = ems_exp3$est_sexRatio_par[n_em] # whether to estimate sex Ratios
    share_M_sex_em = ems_exp3$share_M_sex[n_em] # wether to share M or not
    sexRatio_al_or_y_em = ems_exp3$sexRatio_al_or_y[n_em] # if we want to fit sex ratio as within year only or both
    sexRatio_fixed = c(0.5, 0.5) # fixed sex ratio values to mis-specify initial sex-ratios
    em_name = ems_exp3$EM_Name[n_em] # em name
    
    # estimate biologicals
    biologicals = get_biologicals(n_sexes, n_ages, age_bins, len_bins, Srv_LAA, Srv_LW, sim = sim_num[sim])
    
    # Prepare EM inputs into assessment
    em_inputs = prepare_EM_inputs(sim = sim_num[sim],
                                  # EM_Parameterization
                                  use_fish_sexRatio = use_fish_sexRatio_em,
                                  use_srv_sexRatio = use_srv_sexRatio_em,
                                  est_sexRatio_par = est_sexRatio_par_em,
                                  sexRatio = sexRatio_fixed,
                                  share_M_sex = share_M_sex_em, 
                                  sexRatio_al_or_y = sexRatio_al_or_y_em,
                                  
                                  # Proportion treatment
                                  fish_age_prop = fish_age_prop_em,
                                  srv_age_prop = srv_age_prop_em,
                                  fish_len_prop = fish_len_prop_em,
                                  srv_len_prop = srv_len_prop_em,
                                  
                                  # Fixed controls
                                  n_sexes = 2,
                                  sex_specific = TRUE, 
                                  fit_sexsp_catch = FALSE,
                                  selex_type = "length",
                                  # Aggregating comps
                                  agg_fish_age = FALSE,
                                  agg_srv_age = FALSE, 
                                  agg_fish_len = FALSE,
                                  agg_srv_len = FALSE,
                                  catch_cv = c(0.01), 
                                  use_fish_index = FALSE,
                                  # Biologicals
                                  WAA = biologicals$waa_sex,
                                  age_len_transition = biologicals$al_matrix_sexsp,
                                  # Parameter fixing
                                  fix_pars = c("h", "ln_sigmaRec", "ln_q_fish"))
    
    em_inputs$map$ln_M = factor(c(NA, 1)) # profile M 
    like_df_all <- data.frame() # empty dataframe for storage
    
    for(i in 1:length(m_vary)) {
      em_inputs$parameters$ln_M[1] = m_vary[i] # input varying M here
      TMB::compile("Sex_Str_EM.cpp")
      dyn.load(dynlib('Sex_Str_EM'))
      # run model here
      model = run_model(data = em_inputs$data, 
                        parameters = em_inputs$parameters, 
                        map = em_inputs$map, silent = TRUE, n.newton = 1)
      
      # put data into dataframe
      like_df = data.frame(nLL = c(model$rep$jnLL, 
                                   model$rep$rec_nLL,
                                   sum(model$rep$fish_age_comp_nLL[,1,1]),
                                   sum(model$rep$fish_age_comp_nLL[,2,1]),
                                   sum(model$rep$srv_age_comp_nLL[,1,1]),
                                   sum(model$rep$srv_age_comp_nLL[,2,1]),
                                   sum(model$rep$catch_nLL),
                                   sum(model$rep$srv_index_nLL),
                                   sum(model$rep$srv_len_comp_nLL[,1,1]),
                                   sum(model$rep$srv_len_comp_nLL[,2,1]),
                                   sum(model$rep$fish_len_comp_nLL[,1,1]),
                                   sum(model$rep$fish_len_comp_nLL[,1,1]),
                                   sum( model$rep$sexRatio_fishage_nLL[,1,1]),
                                   sum( model$rep$sexRatio_srvage_nLL[,1,1]),
                                   sum( model$rep$sexRatio_fishlen_nLL[,1,1]),
                                   sum( model$rep$sexRatio_srvlen_nLL[,1,1])),
                           comp = c("Joint", "Rec", "FishAgeF", "FishAgeM",
                                    "SrvAgeF", "SrvAgeM", "Catch", "SrvIndx", 
                                    "SrvLenF", "SrvLenM", "FishLenF", "FishLenM",
                                    "SRFishAge", "SRSrvAge", "SRFishLen", "SRSrvLen"),
                           par = m_vary[i], sim = sim_num[sim], type = "M_F", em = em_name)
      
      like_df_all <- rbind(like_df_all, like_df)
      print(i)
    } # end i for likelihood profile loop
    like_df_all
  } # end sim for simulation loop
  m_em_list[[n_em]] <- profiles
} # end em loop

mf_profile_df <- bind_rows(m_em_list) # output list of profiled values

# bind all df here
all_nLL_df = rbind(mf_profile_df, r0_profile_df)
write.csv(all_nLL_df, here("output", "Experiment3_Likelihood_Prof.csv"))

# Plot ll profiles --------------------------------------------------------

# read in data
all_nLL_df <- read.csv(here("output", "Experiment3_Likelihood_Prof.csv"))

# Get summary statistics for likelihood profiles
r0_nLL_plot_df <- all_nLL_df %>%
  filter(type == "R0") %>% 
  filter(par <= 5.25, par >= 1.5) %>% 
  group_by(comp, sim, type, em) %>%
  mutate(nLL = nLL - min(nLL)) %>% 
  ungroup() %>% 
  group_by(comp, par, type, em) %>%
  summarize(
    median_nLL = median(nLL),
    lwr = quantile(nLL, 0.025),
    upr = quantile(nLL, 0.975)
  )

# get minimum of each data component
r0_min_nLL_df <- r0_nLL_plot_df %>%
  group_by(comp, type, em) %>%
  slice_min(order_by = median_nLL) %>%
  select(comp, par, type, em) %>%
  rename(parmin = par)

# join these two together
r0_nLL_plot_df <- r0_nLL_plot_df %>% left_join(r0_min_nLL_df, by = c("comp", "type", "em"))

# Get summary statistics for likelihood profiles
m_nLL_plot_df <- all_nLL_df %>%
  filter(type == "M_F") %>% 
  group_by(comp, sim, type, em) %>%
  mutate(nLL = nLL - min(nLL)) %>% 
  ungroup() %>% 
  group_by(comp, par, type, em) %>%
  summarize(
    median_nLL = median(nLL),
    lwr = quantile(nLL, 0.025),
    upr = quantile(nLL, 0.975)
  )

# get minimum of each data component
m_min_nLL_df <- m_nLL_plot_df %>%
  group_by(comp, type, em) %>%
  slice_min(order_by = median_nLL) %>%
  select(comp, par, type, em) %>%
  rename(parmin = par)

# join these two together
m_nLL_plot_df <- m_nLL_plot_df %>% left_join(m_min_nLL_df, by = c("comp", "type", "em"))
nLL_plot_df <- rbind(m_nLL_plot_df, r0_nLL_plot_df)

pdf(here("figs", "Experiment 3", "Likelihood_Profile.pdf"), width = 12)
ggplot(nLL_plot_df %>% filter(!str_detect(comp, "SR")), 
       aes(x = exp(par), y = median_nLL, color = comp)) +
  geom_line(size = 1.5, alpha = 0.85) +
  geom_vline(aes(xintercept = exp(parmin), color = comp), lty = 2, lwd = 1.2, alpha = 0.75) +
  coord_cartesian(ylim = c(0,100)) +
  labs(x = "Parameter Value", y = "Median nLL", color = "Data Component") +
  theme_tj() +
  facet_grid(em~type, scales = "free") +
  theme(legend.position = "top")
dev.off()


ggplot(nLL_plot_df %>% filter(!str_detect(comp, "SR")), 
       aes(x = exp(par), y = median_nLL, color = comp)) +
  geom_line(size = 1.5, alpha = 0.85) +
  geom_vline(aes(xintercept = exp(parmin), color = comp), lty = 2, lwd = 1.2, alpha = 0.75) +
  coord_cartesian(ylim = c(0,100)) +
  labs(x = "Parameter Value", y = "Median nLL", color = "Data Component") +
  theme_tj() +
  facet_grid(em~type, scales = "free") +
  theme(legend.position = "top")

