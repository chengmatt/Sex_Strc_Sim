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
                F_pattern = "Increase",
                comp_across_sex = "across",
                selex_type = "length",
                q_Fish = 0.025,
                cv_Fish_Index = 0.25,
                q_Srv = 0.05,
                cv_Srv_Index = 0.25)
  
  plot(FishAge_Selex[,1,1])
  lines(FishAge_Selex[,2,1])
  plot(SrvAge_Selex[,1,1])
  lines(SrvAge_Selex[,2,1])
  plot(NAA[1,,1,5])
  
  # Get Length weight samples
  Srv_LAA = data.table::rbindlist(Srv_LAA)
  Srv_LW = data.table::rbindlist(Srv_LW)
  
  ggplot(Srv_LAA %>% filter(sim == 1), aes(x = ages, y = lens, color = factor(sex))) +
    geom_point() 
  ggplot(Srv_LW %>% filter(sim == 1), aes(x = lens, y = wts, color = factor(sex))) +
    geom_point() 
  
  # TMB Testing -------------------------------------------------------------
  
  library(TMB)
  # setwd("src")
  TMB::compile("Sex_Str_EM.cpp")
  dyn.unload(dynlib('Sex_Str_EM'))
  dyn.load(dynlib('Sex_Str_EM'))
  
  totalrec_all = data.frame()
  totalbiom_all = data.frame()
  ssb_all = data.frame()
  m_all = data.frame()
  rec_all = data.frame()
  fmsy_all = data.frame()
  
  for(sim in 1:n_sims) {
    
    # get biological information
    biologicals = get_biologicals(n_sexes, n_ages, age_bins, len_mids, Srv_LAA, Srv_LW, sim = sim)
    
    em_inputs = prepare_EM_inputs(sim = sim,
                                  sexRatio = c(0.5, 0.5),
                                  catch_cv = c(1e-3),
                                  WAA = waa,
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
                                  selex_type = "length",
                                  fix_pars = c("h", "ln_sigmaRec"))
    
    model_fxn = TMB::MakeADFun(em_inputs$data, em_inputs$parameters, em_inputs$map, random = NULL,
                               DLL= "Sex_Str_EM", silent = TRUE,  
                               checkParameterOrder = TRUE, tracepar = TRUE)
    
    # Optimize model here w/ nlminb
    mle_optim <- stats::nlminb(model_fxn$par, model_fxn$fn, model_fxn$gr,
                               control = list(iter.max = 1e5, eval.max = 1e5))
    add_newton(n.newton = 3, ad_model = model_fxn, mle_optim = mle_optim)
    model_fxn$rep <- model_fxn$report(model_fxn$env$last.par.best) # Need to pass both fixed and random effects!!!
    sd_rep <- TMB::sdreport(model_fxn)
    Report = model_fxn$rep # get report
    
    if(sum(is.na(sd_rep$sd)) == 0) {
      SSBres = data.frame(Pred = Report$SSB, True = SSB[-n_years,sim], years = 1:length(Report$SSB))
      ssb_all = rbind(SSBres, ssb_all)
      TotalBiomres = data.frame(Pred = Report$Total_Biom, True = Total_Biom[-n_years,sim], years = 1:length(Report$Total_Biom))
      totalbiom_all = rbind(TotalBiomres, totalbiom_all)
      TotalRecres = data.frame(Pred = Report$Total_Rec, True = rowSums(NAA[-n_years,1,,sim]), years = 1:length(Report$Total_Biom))
      totalrec_all = rbind(totalrec_all, TotalRecres)
      Mest = data.frame(Pred = exp(sd_rep$par.fixed[names(sd_rep$par.fixed) == "ln_M"]), True = M, Sex = rep(c("F", "M")))
      recest = data.frame(Pred = exp(sd_rep$par.fixed[names(sd_rep$par.fixed) == "RecPars"]), True = r0)
      rec_all = rbind(rec_all, recest)
      m_all = rbind(Mest, m_all)
      fmsy_est = data.frame(Pred = exp(sd_rep$par.fixed[names(sd_rep$par.fixed) == "ln_Fmsy"]), True = fmsy[sim])
      fmsy_all <- rbind(fmsy_est, fmsy_all)
    }
    
    print(sim)
  } # end sim
  
plot(Report$Fish_Slx[1,,1,1])
lines(FishAge_Selex[,1,1])
plot(Report$Fish_Slx[1,,2,1])
lines(FishAge_Selex[,2,1])

plot(Report$Srv_Slx[1,,1,1])
lines(SrvAge_Selex[,1,1])
plot(Report$Srv_Slx[1,,2,1])
lines(SrvAge_Selex[,2,1])

ggplot(m_all, aes(x = Sex, y = (Pred - True) / True)) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  coord_cartesian(ylim = c(-0.3, 0.3))

ggplot(rec_all, aes(x = "r0", y = (Pred - True) / True)) +
  geom_boxplot() +
  geom_hline(
    yintercept = 0) +
  coord_cartesian(ylim = c(-0.5, 0.5))

ggplot(fmsy_all, aes(x = "fmsy", y = (Pred - True) / True)) +
  geom_boxplot() +
  geom_hline(
    yintercept = 0) +
  coord_cartesian(ylim = c(-0.5, 0.5))

ssb_all %>% 
  group_by(years) %>% 
  mutate(RE = (Pred - True) / True) %>% 
  summarize(median = median(RE),
            lwr95 = quantile(RE, 0.025),
            upr95 = quantile(RE, 0.975)) %>% 
  ggplot(aes(x = years, y = median, ymin = lwr95, ymax = upr95)) +
  geom_line() +
  geom_ribbon(alpha = 0.3) +
  geom_hline(yintercept = 0) +
  coord_cartesian(ylim = c(-0.5, 0.5))

totalbiom_all %>% 
  group_by(years) %>% 
  mutate(RE = (Pred - True) / True) %>% 
  summarize(median = median(RE),
         lwr95 = quantile(RE, 0.025),
         upr95 = quantile(RE, 0.975)) %>% 
  ggplot(aes(x = years, y = median, ymin = lwr95, ymax = upr95)) +
  geom_line() +
  geom_ribbon(alpha = 0.3) +
  geom_hline(yintercept = 0) +
  coord_cartesian(ylim = c(-0.5, 0.5))

totalrec_all %>% 
  group_by(years) %>% 
  mutate(RE = (Pred - True) / True) %>% 
  summarize(median = median(RE),
            lwr95 = quantile(RE, 0.025),
            upr95 = quantile(RE, 0.975)) %>% 
  ggplot(aes(x = years, y = median, ymin = lwr95, ymax = upr95)) +
  geom_line() +
  geom_ribbon(alpha = 0.3) +
  geom_hline(yintercept = 0) +
  coord_cartesian(ylim = c(-0.5, 0.5))

