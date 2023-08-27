# Purpose: To simulate datasets for use in sex-structured simulations
# Author: Matthew LH. Cheng (UAF - CFOS)
# Date: 8/2/23

library(here)
library(tidyverse)

# Load in all functions from the functions folder
fxn_path <- here("R", "functions")
files <- list.files(fxn_path)
for(i in 1:length(files)) source(here(fxn_path, files[i]))

simulate_data(spreadsheet_path = here("input", "Sablefish_Inputs.xlsx"),
              Fish_Neff_Age = 100,
              Fish_Neff_Len = 100,
              Srv_Neff_Age = 100,
              Srv_Neff_Len = 100,
              comp_across_sex = "across",
              q_Fish = 0.025,
              cv_Fish_Index = 0.25,
              q_Srv = 0.05,
              cv_Srv_Index = 0.25)

# Get Length weight samples
Srv_LAA = data.table::rbindlist(Srv_LAA)
Srv_LW = data.table::rbindlist(Srv_LW)

ggplot(Srv_LAA %>% filter(sim == 1), aes(x = ages, y = lens, color = factor(sex))) +
  geom_point() 
ggplot(Srv_LW %>% filter(sim == 1), aes(x = lens, y = wts, color = factor(sex))) +
  geom_point() 

plot(get_WAA(LAA_obs_age = Srv_LAA$ages[Srv_LAA$sim == 1], 
             LAA_obs_len = Srv_LAA$lens[Srv_LAA$sim == 1],
             WL_obs_len = Srv_LW$lens[Srv_LW$sim == 1],
             WL_obs_wt = Srv_LW$wts[Srv_LW$sim == 1],
             ages = age_bins)[[1]])
lines(waa[,1])

# TMB Testing -------------------------------------------------------------

library(TMB)
# setwd("src")
TMB::compile("Sex_Str_EM.cpp")
dyn.unload(dynlib('Sex_Str_EM'))
dyn.load(dynlib('Sex_Str_EM'))

em_inputs = prepare_EM_inputs(sim = 1,
                              sexRatio = c(0.5, 0.5),
                              catch_cv = c(1e-3),
                              WAA = matrix(waa[,], ncol = 2, nrow = 30),
                              age_len_transition = al_matrix,
                              n_sexes = 2,
                              fish_age_prop = "across",
                              srv_age_prop = "across",
                              fish_len_prop = "across",
                              srv_len_prop = "across",
                              agg_fish_age = FALSE,
                              agg_srv_age = FALSE,
                              share_M_sex = FALSE,
                              sex_specific = TRUE,
                              fix_pars = c("h", "ln_sigmaRec"))

model_fxn = TMB::MakeADFun(em_inputs$data, em_inputs$parameters, em_inputs$map, random = NULL,
                           DLL= "Sex_Str_EM", silent = FALSE,  
                           checkParameterOrder = TRUE, tracepar = TRUE)

Opt = TMBhelper::fit_tmb( obj = model_fxn,
                          newtonsteps = 3,
                          bias.correct = FALSE,
                          getsd = TRUE,
                          savedir = NULL)

Report = model_fxn$report(model_fxn$env$last.par.best)
sum(Report$pred_srv_len_comps[1,,,1])

plot(Report$SSB, type = "l", col = "blue")
lines(SSB[-n_years,1])
exp(Opt$SD$par.fixed[names(Opt$SD$par.fixed) == "ln_M"])

ParHat = model_fxn$env$parList()
exp(Opt$SD$par.fixed[names(Opt$SD$par.fixed) == "ln_q_fish"])
exp(Opt$SD$par.fixed[names(Opt$SD$par.fixed) == "ln_q_srv"])
exp(Opt$SD$par.fixed[names(Opt$SD$par.fixed) == "ln_M"])
plot(exp(Opt$SD$par.fixed[names(Opt$SD$par.fixed) == "ln_RecDevs"]), type = "l")
lines(RecDevs[,1], col = "red")

plot(Report$NAA[50,,1], type = "l", col = "black")
# lines(Report$NAA[2,,1], type = "l", col = "black")
lines(NAA[50,,1,1])

par(mfrow = c(2,2))
plot(Report$NAA[1,,1], type = "l", col = "red")
plot(Report$NAL[50,,1], type = "l")
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
plot(Report$pred_catch[,1])
lines(em_inputs$data$obs_catch[,1])

plot(Report$Fish_Slx[1,,1,1], type = "l", col = "blue")
lines(FishAge_Selex[,1,1])

plot(Report$Srv_Slx[1,,1,1], type = "l", col = "blue")
lines(SrvAge_Selex[,1,1])

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


