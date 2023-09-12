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
                Fish_Neff_Age = 30,
                Fish_Neff_Len = 30,
                Srv_Neff_Age = 30,
                Srv_Neff_Len = 30,
                F_pattern = "Contrast",
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
  bmsy_all = data.frame()
  all_profiles = data.frame()
  
  for(sim in 1:n_sims) {
    
    # get biological information
    biologicals = get_biologicals(n_sexes, n_ages, age_bins, len_mids, Srv_LAA, Srv_LW, sim = sim)
    
    em_inputs = prepare_EM_inputs(sim = sim,
                                  sexRatio = c(0.5, 0.5),
                                  catch_cv = c(1e-3),
                                  WAA = waa,
                                  age_len_transition = al_matrix,
                                  n_sexes = 2,
                                  fish_age_prop = "within",
                                  srv_age_prop = "within",
                                  fish_len_prop = "within",
                                  srv_len_prop = "within",
                                  agg_fish_age = FALSE,
                                  agg_srv_age = FALSE,
                                  share_M_sex = FALSE,
                                  sex_specific = TRUE,
                                  selex_type = "length",
                                  fix_pars = c("h", "ln_sigmaRec"))
    
    # run model here
    models = run_model(data = em_inputs$data, parameters = em_inputs$parameters, map = em_inputs$map)
    
    # if(sim %in% c(seq(1, n_sims, 10))) {
    #   profiles = like_prof(em_inputs = em_inputs, sim = sim, assessment_name = "a", 
    #                        share_M = FALSE, sex_specific = TRUE)
    #   all_profiles = rbind(profiles, all_profiles)
    # } # only do profiles for every tenth simulation
    
    # plot(NAA[30,,1,sim-1])
    # lines(models$rep$NAA[30,,1])
    
    if(sum(is.na(models$sd_rep$sd)) == 0) {
      SSBres = data.frame(Pred = models$rep$SSB, True = SSB[-n_years,sim], years = 1:length(models$rep$SSB), sim = sim)
      ssb_all = rbind(SSBres, ssb_all)
      TotalBiomres = data.frame(Pred = models$rep$Total_Biom, True = Total_Biom[-n_years,sim], years = 1:length(models$rep$Total_Biom), sim = sim)
      totalbiom_all = rbind(TotalBiomres, totalbiom_all)
      bmsy_df = data.frame(Pred = models$rep$BMSY, True = bmsy[sim], sim = sim)
      bmsy_all = rbind(bmsy_all, bmsy_df)
      TotalRecres = data.frame(Pred = models$rep$Total_Rec, True = rowSums(NAA[-n_years,1,,sim]), years = 1:length(models$rep$Total_Biom), sim = sim)
      totalrec_all = rbind(totalrec_all, TotalRecres)
      Mest = data.frame(Pred = exp(models$sd_rep$par.fixed[names(models$sd_rep$par.fixed) == "ln_M"]), True = M, Sex = rep(c("F", "M")), sim = sim)
      m_all = rbind(Mest, m_all)
      recest = data.frame(Pred = exp(models$sd_rep$par.fixed[names(models$sd_rep$par.fixed) == "RecPars"]), True = r0, sim = sim)
      rec_all = rbind(rec_all, recest)
      fmsy_est = data.frame(Pred = exp(models$sd_rep$par.fixed[names(models$sd_rep$par.fixed) == "ln_Fmsy"]), True = fmsy[sim], sim = sim)
      fmsy_all <- rbind(fmsy_est, fmsy_all)
    }
    
    print(sim)
  } # end sim

m_all_df = m_all %>% mutate(RE = (Pred - True) / True, type = paste("M", Sex, sep = "_")) %>% select(-Sex)
rec_all_df = rec_all %>% mutate(RE = (Pred - True) / True, type = "R0")
fmsy_all_df = fmsy_all %>% mutate(RE = (Pred - True) / True, type = "Fmsy")
bmsy_all_df = bmsy_all %>% mutate(RE = (Pred - True) / True, type = "bmsy")
par_all = rbind(bmsy_all_df, m_all_df, rec_all_df, fmsy_all_df) %>% 
  select(RE, type, sim) %>% 
  pivot_wider(names_from = "type", values_from = "RE") %>% select(-sim)

# ssb df
ssb_df = ssb_all %>% 
  group_by(years) %>% 
  mutate(RE = (Pred - True) / True) %>% 
  summarize(median = median(RE),
            lwr95 = quantile(RE, 0.025),
            upr95 = quantile(RE, 0.975),
            lwr50 = quantile(RE, 0.25),
            upr50 = quantile(RE, 0.75))

# total biomass
totalbiom_all_df = totalbiom_all %>% 
  group_by(years) %>% 
  mutate(RE = (Pred - True) / True) %>% 
  summarize(median = median(RE),
         lwr95 = quantile(RE, 0.025),
         upr95 = quantile(RE, 0.975),
         lwr50 = quantile(RE, 0.25),
         upr50 = quantile(RE, 0.75))

totalrec_all_df = totalrec_all %>% 
  group_by(years) %>% 
  mutate(RE = (Pred - True) / True) %>% 
  summarize(median = median(RE),
            lwr95 = quantile(RE, 0.025),
            upr95 = quantile(RE, 0.975),
            lwr50 = quantile(RE, 0.25),
            upr50 = quantile(RE, 0.75))

par(mar=c(4,6,1,2))

nf <- layout(
  matrix(c(1,2,3,rep(4, 3)), ncol=3, byrow=TRUE), 
  widths=c(1,1), 
  heights=c(2,2)
)
plot_re_ts(ssb_df, ylab = "Relative Error in SSB")
plot_re_ts(totalbiom_all_df, ylab = "Relative Error in Total Biomass")
plot_re_ts(totalrec_all_df, ylab = "Relative Error in Total Recruitment")
plot_re_par(par_all)

# all_profiles %>%
#   pivot_longer(!c(prof_val, om_val, par_name, assessment_name, sim),
#                names_to = "type", values_to = "nLL") %>%
#   drop_na() %>% 
#   group_by(type, par_name, sim) %>%
#   mutate(nLL = nLL - min(nLL, na.rm = TRUE)) %>%
#   ungroup() %>% 
#   group_by(type, par_name, prof_val) %>% 
#   mutate(median = median(nLL)) %>% 
#   filter(case_when(
#     str_detect(par_name,"ln_M") ~ prof_val < log(0.185), TRUE ~ TRUE )) %>%  
#   filter(par_name == "RecPars_1") %>%
#   ggplot(aes(x = exp(prof_val), y = nLL, group = sim)) +
#   geom_line(size = 1.25, alpha = 95, color = "grey") +
#   geom_line(aes(y = median)) +
#   geom_vline(aes(xintercept = exp(om_val)), lty = 2) +
#   facet_wrap(~type, scale = "free_x") +
#   theme(legend.position = "top")
