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

  ggplot(Srv_LAA %>% filter(sim == 1), aes(x = ages, y = lens, color = factor(sex))) +
    geom_point() 
  ggplot(Srv_LW %>% filter(sim == 1), aes(x = lens, y = wts, color = factor(sex))) +
    geom_point() 
  
  # TMB Testing -------------------------------------------------------------
  
  library(TMB)
  setwd("src")
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
  Req_all = data.frame()
  fmort_all = data.frame()
  selex_all = data.frame()
  all_profiles = data.frame()
  
  # plot(Fish_AgeComps[1,,1,1,1], type = "l", col = "red")
  # lines(Fish_AgeComps[1,,2,1,1], type = "l", col = "blue")
  # plot(Srv_AgeComps[1,,1,1,1], type = "l", col = "red")
  # lines(Srv_AgeComps[1,,2,1,1], type = "l", col = "blue")
  # plot(Fish_LenComps[1,,1,1,1], type = "l", col = "red")
  # lines(Fish_LenComps[1,,2,1,1], type = "l", col = "blue")
  # plot(Srv_LenComps[1,,1,1,1]/sum(Srv_LenComps[1,,1,1,1]), type = "l", col = "red")
  # lines(Srv_LenComps[1,,2,1,1]/sum(Srv_LenComps[1,,2,1,1]), type = "l", col = "blue")
  
  for(sim in 1:n_sims) {
    
    # get biological information
    biologicals = get_biologicals(n_sexes, n_ages, age_bins, len_mids, Srv_LAA, Srv_LW, sim = sim)
    
    em_inputs = prepare_EM_inputs(sim = sim,
                                  sexRatio = c(0.5, 0.5),
                                  catch_cv = c(1e-2),
                                  WAA = waa,
                                  age_len_transition = al_matrix,
                                  n_sexes = 2,
                                  fish_age_prop = "across",
                                  srv_age_prop = "across",
                                  fish_len_prop = "across",
                                  srv_len_prop = "across",
                                  agg_fish_age = FALSE,
                                  agg_srv_age = FALSE, 
                                  agg_fish_len = FALSE,
                                  agg_srv_len = FALSE,
                                  share_M_sex = FALSE,
                                  sex_specific = TRUE, 
                                  fit_sexsp_catch = TRUE,
                                  selex_type = "age",
                                  fix_pars = c("h", "ln_sigmaRec"))
    
    # run model here
    models = run_model(data = em_inputs$data, 
                       parameters = em_inputs$parameters, 
                       map = em_inputs$map, silent = TRUE, n.newton = 5)

    plot(models$rep$pred_catch_sexsp[,1,1], col = "red", type = "l")
    lines(Total_Catch_Sex[-31,1,1,sim], col = "red")
    plot(models$rep$pred_catch_sexsp[,2,1], col = 'blue', type = "l")
    lines(Total_Catch_Sex[-31,2,1,sim], col = "blue")
    
    plot(rowSums(Total_Catch_Sex[-31,,1,sim]))
    lines(models$rep$pred_catch_agg)
    models$rep$catch_nLL

    # plot(models$rep$pred_srv_len_agg)
    # lines(models$rep$obs_srv_len_agg / sum(models$rep$obs_srv_len_agg))
    # 
    # plot(models$rep$pred_fish_len_agg)
    # lines(models$rep$obs_fish_len_agg / sum(models$rep$obs_fish_len_agg))
    
    
    # options(mc.cores = 4)
    # models_stan = tmbstan::tmbstasn(models, init = "last.par.best")
    # rstan::traceplot(models_stan) # trace plot
    # vals = rstan::extract(models_stan)
    # rstan::stan_rhat(models_stan) # rhats
    # shinystan::launch_shinystan(models_stan) # shiny stan
    # pairs(models_stan, pars=c("ln_M", "ln_q_srv", "RecPars"))
    
    
    # plot(models$rep$pred_srv_age_as/sum(models$rep$pred_srv_age_as), type = 'l')
    # lines(models$rep$obs_srv_age_as/sum(models$rep$obs_srv_age_as), type = 'l', col = "red")
    # lines(as.vector(CAA[50,,,1,sim]/sum(CAA[50,,,1,sim])), type = 'l')
    # plot(models$rep$Fish_Slx[1,,1,1])
    # lines(models$rep$Fish_Slx[1,,2,1])
    
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
      bmsy_df = data.frame(Pred = models$rep$BMSY, True = bmsy, sim = sim)
      bmsy_all = rbind(bmsy_all, bmsy_df)
      TotalRecres = data.frame(Pred = models$rep$Total_Rec, True = rowSums(NAA[-n_years,1,,sim]), years = 1:length(models$rep$Total_Biom), sim = sim)
      totalrec_all = rbind(totalrec_all, TotalRecres)
      Mest = data.frame(Pred = exp(models$sd_rep$par.fixed[names(models$sd_rep$par.fixed) == "ln_M"]), True = M, Sex = rep(c("F", "M")), sim = sim)
      m_all = rbind(Mest, m_all)
      recest = data.frame(Pred = exp(models$sd_rep$par.fixed[names(models$sd_rep$par.fixed) == "RecPars"]), True = r0, sim = sim)
      rec_all = rbind(rec_all, recest)
      fmsy_est = data.frame(Pred = exp(models$sd_rep$par.fixed[names(models$sd_rep$par.fixed) == "ln_Fmsy"]), True = fmsy, sim = sim)
      fmsy_all <- rbind(fmsy_est, fmsy_all)
      Reqres = data.frame(Pred = models$rep$Req, True = Req, sim = sim)
      Req_all = rbind(Reqres, Req_all)
      fmort = data.frame(Pred = exp(models$sd_rep$par.fixed[names(models$sd_rep$par.fixed) == "ln_Fy"]), 
                         True = Fmort[-n_years,1,sim], years = 1:length(models$rep$Total_Biom), sim = sim)
      fmort_all = rbind(fmort_all, fmort)
      
      # Save selex estimates
      selex_f = data.frame(Pred = models$rep$Fish_Slx[1,,1,1], True = FishAge_Selex[,1,1], sim = sim, sex = "F", age = age_bins)
      selex_m = data.frame(Pred = models$rep$Fish_Slx[1,,2,1], True = FishAge_Selex[,2,1], sim = sim, sex = "M", age = age_bins)
      selex_all = rbind(selex_f, selex_m, selex_all)
      
    }
    
    print(sim)
    
  } # end sim

# plot(models$rep$Fish_Slx[20,,1,1], type = "l", col = "red")
# lines(models$rep$Fish_Slx[20,,2,1], type = "l", col = "blue")
# lines(FishAge_Selex[,1,1])
# lines(FishAge_Selex[,2,1])

m_all_df = m_all %>% mutate(RE = (Pred - True) / True, type = paste("M", Sex, sep = "_")) %>% select(-Sex)
rec_all_df = rec_all %>% mutate(RE = (Pred - True) / True, type = "R0")
fmsy_all_df = fmsy_all %>% mutate(RE = (Pred - True) / True, type = "Fmsy")
bmsy_all_df = bmsy_all %>% mutate(RE = (Pred - True) / True, type = "Bmsy")
Req_all_df = Req_all %>% mutate(RE = (Pred - True) / True, type = "Req")

par_all = rbind(Req_all_df, bmsy_all_df, m_all_df, rec_all_df, fmsy_all_df) %>% 
  select(RE, type, sim) %>%  pivot_wider(names_from = "type", values_from = "RE") %>%
  select(-sim)

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

fmort_all_df = fmort_all %>% 
  group_by(years) %>% 
  mutate(RE = (Pred - True) / True) %>% 
  summarize(median = median(RE),
            lwr95 = quantile(RE, 0.025),
            upr95 = quantile(RE, 0.975),
            lwr50 = quantile(RE, 0.25),
            upr50 = quantile(RE, 0.75))

fmort_all_df = fmort_all %>% 
  group_by(years) %>% 
  mutate(RE = (Pred - True) / True) %>% 
  summarize(median = median(RE),
            lwr95 = quantile(RE, 0.025),
            upr95 = quantile(RE, 0.975),
            lwr50 = quantile(RE, 0.25),
            upr50 = quantile(RE, 0.75))

ggplot(selex_all, aes(x = age, y = Pred, group = sim)) +
  geom_line(color = "grey") +
  geom_line(aes(y = True, color = sex), size = 1.5) +
  facet_wrap(~sex) +
  labs(x = "Age", y = "Selex") + 
  ggthemes::scale_color_excel_new() +
  ggthemes::theme_excel_new()

selex_all%>% 
  group_by(sex, age) %>% 
  mutate(RE = (Pred - True) / True) %>% 
  summarize(median = median(RE),
            lwr95 = quantile(RE, 0.025),
            upr95 = quantile(RE, 0.975),
            lwr50 = quantile(RE, 0.25),
            upr50 = quantile(RE, 0.75)) %>% 
  ggplot(aes(x = age, y = median, ymin = lwr95, ymax = upr95)) +
  geom_line() +
  geom_ribbon(alpha = 0.5) +
  facet_wrap(~sex) +
  geom_hline(yintercept = 0) +
  coord_cartesian(ylim = c(-0.75, 0.75)) 
  # ggthemes::scale_color_excel_new() +
  # ggthemes::theme_excel_new()
  
par(mar=c(4,6,1,2))

nf <- layout(
  matrix(c(1,2,3,4,rep(5, 4)), ncol=4, byrow=TRUE), 
  widths=c(1,1), 
  heights=c(2,2)
)

plot_re_ts(ssb_df, ylab = "Relative Error in SSB")
plot_re_ts(totalbiom_all_df, ylab = "Relative Error in Total Biomass")
plot_re_ts(totalrec_all_df, ylab = "Relative Error in Total Recruitment")
plot_re_ts(fmort_all_df, ylab = "Relative Error in Total Fishing Mortality")
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
