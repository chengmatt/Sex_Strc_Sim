  # Purpose: To simulate datasets for use in sex-structured simulations
  # Author: Matthew LH. Cheng (UAF - CFOS)
  # Date: 8/2/23
  
  library(here)
  library(tidyverse)
  library(compResidual)
  
  # Load in all functions from the functions folder
  fxn_path <- here("R", "functions")
  files <- list.files(fxn_path)
  for(i in 1:length(files)) source(here(fxn_path, files[i]))
  
  oms = simulate_data(spreadsheet_path = here("input", "Sablefish_Inputs_test.xlsx"),
                Fish_Neff_Age = 200,
                Fish_Neff_Len = 200,
                Srv_Neff_Age = 200,
                Srv_Neff_Len = 200,
                F_pattern = "Contrast",
                comp_across_sex = "across",
                selex_type = "length",
                q_Fish = 0.025,
                cv_Fish_Index = 0.25,
                q_Srv = 0.05,
                sexRatio = c(0.5, 0.5),
                cv_Srv_Index = 0.1, 
                growth_control = "chg_males_rel_females",
                natmort_control = "chg_males_rel_females",
                growth_control_fct = 0.9, 
                natmort_control_fct = 0.9,
                force_grwth_same_yng = FALSE
                )
  
  plot(vonB[,1])
  lines(vonB[,2])
  
  # get biological information
  biologicals = get_biologicals(n_sexes, n_ages, age_bins, len_mids, Srv_LAA, Srv_LW, sim = 1)
  plot(biologicals$al_matrix_sexsp[6,,1])
  lines(biologicals$al_matrix_sexagg[1,,1])
  
  plot(oms$FishAge_Selex[,1,], ylim = c(0,1))
  lines(oms$FishAge_Selex[,2,])
  
  plot(oms$SrvAge_Selex[,1,], ylim = c(0,1))
  lines(oms$SrvAge_Selex[,2,])
  
  # HCR_proj_catch[sim]
  get_proj_catch(fmsy_val = fmsy, bmsy_val = bmsy, sex_ratio = 1, 
                 n_ages = n_ages, n_sexes = 1, term_NAA = rowSums(NAA[n_years-1,,,sim]), 
                 term_SSB = SSB[n_years-1, sim], term_F_Slx = FishAge_Selex, r0 = r0,
                 term_F = Fmort[n_years - 1,,sim], M_s = M, WAA = waa, MatAA = mat_at_age)
  
  get_proj_catch(fmsy_val = fmsy, bmsy_val = bmsy, sex_ratio = sexRatio, 
                 n_ages = n_ages, n_sexes = n_sexes, term_NAA = NAA[n_years-1,,,sim], 
                 term_SSB = SSB[n_years-1, sim], term_F_Slx = FishAge_Selex, r0 = r0,
                 term_F = Fmort[n_years - 1,,sim], M_s = M, WAA = waa, MatAA = mat_at_age)
  
  ggplot(Srv_LAA %>% filter(sim == 8), aes(x = ages, y = lens, color = factor(sex))) +
    geom_point(alpha = 0.05) +
    # facet_wrap(~sex) +
    geom_smooth()
  
  # ggplot(Srv_LW %>% filter(sim == 1), aes(x = lens, y = wts, color = factor(sex))) +
  #   geom_point()
  
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
  Req_all = data.frame()
  fmort_all = data.frame()
  selex_all = data.frame()
  all_profiles = data.frame()
  sexratio = vector()
  sbpr_msy = vector()
  est_hcr_catch = vector()
  conv = vector()
  list2env(oms,globalenv()) # output into global environment
  
  b = reshape2::melt(cov2cor(models$sd_rep$cov.fixed)) %>% 
    filter(value != 1)
  
  for(sim in 1:n_sims) {
    
    # get biological information
    biologicals = get_biologicals(n_sexes, n_ages, age_bins, len_mids, Srv_LAA, Srv_LW, sim = sim)
    
    plot(waa[,1])
    # lines(biologicals$waa_nosex)
    plot(biologicals$waa_sex[,1])
    lines(biologicals$waa_sex[,2])
    
    em_inputs = prepare_EM_inputs(sim = sim,
                                  # sex-parameterizations
                                  n_sexes = 1,
                                  sex_specific = FALSE, 
                                  share_M_sex = TRUE,
                                  sexRatio = c(1, 1),
                                  est_sexRatio_par = FALSE,
                                  use_fish_sexRatio = FALSE,
                                  use_srv_sexRatio = FALSE,
                                  fit_sexsp_catch = FALSE,
                                  sexRatio_al_or_y_em = NA,
                                  
                                  # selectivity
                                  selex_type = "length",
                                  
                                  # Biologicals
                                  WAA = biologicals$waa_nosex,
                                  age_len_transition = biologicals$al_matrix_sexagg,
                                  # WAA = biologicals$waa_nosex,
                                  # age_len_transition = biologicals$al_matrix_sexagg,
                                  
                                  # Fishery proportion treatment
                                  fish_age_prop = "within",
                                  srv_age_prop = "within",
                                  fish_len_prop = "within",
                                  srv_len_prop = "within",
                                  
                                  # Aggregating comps
                                  agg_fish_age = FALSE,
                                  agg_srv_age = FALSE, 
                                  agg_fish_len = FALSE,
                                  agg_srv_len = FALSE,
                                  catch_cv = c(1e-2),
                                  use_fish_index = FALSE,

                                  # Parameter fixing
                                  fix_pars = c("h", "ln_sigmaRec", "ln_q_fish"))
    
    # em_inputs$parameters$ln_Fy = matrix(log(Fmort[-n_years,1,,sim]))
    
    # run model here
    models = run_model(data = em_inputs$data, 
                       parameters = em_inputs$parameters, 
                       map = em_inputs$map, silent = T, n.newton = 3)
    
    par(mfrow = c(1,2))
    plot(oms$SrvAge_Selex[,1,1])
    lines(oms$SrvAge_Selex[,2,1])
    lines(models$rep$Srv_Slx[1,,1,1])

    plot(oms$FishAge_Selex[,1,1])
    lines(oms$FishAge_Selex[,2,1])
    lines(models$rep$Fish_Slx[1,,1,1])

    if(models$sd_rep$pdHess == TRUE & max(abs(models$sd_rep$gradient.fixed)) < 0.1) conv[sim] = "Converged"
    else conv[sim] = "Not Converg"
    # if(sim %in% c(seq(1, n_sims, 10))) {
    #   profiles = like_prof(em_inputs = em_inputs, sim = sim, assessment_name = "a",
    #                        share_M = FALSE, sex_specific = TRUE)
    #   all_profiles = rbind(profiles, all_profiles)
    # } # only do profiles for every tenth simulation
    # 
    # if(sum(is.na(models$sd_rep$sd)) == 0) {
    #   SSBres = data.frame(Pred = models$rep$SSB, True = SSB[-n_years,sim], years = 1:length(models$rep$SSB), sim = sim)
    #   ssb_all = rbind(SSBres, ssb_all)
    #   TotalBiomres = data.frame(Pred = models$rep$Total_Biom, True = Total_Biom[-n_years,sim], years = 1:length(models$rep$Total_Biom), sim = sim)
    #   totalbiom_all = rbind(TotalBiomres, totalbiom_all)
    #   TotalRecres = data.frame(Pred = models$rep$Total_Rec, True = rowSums(NAA[-n_years,1,,sim]), years = 1:length(models$rep$Total_Biom), sim = sim)
    #   totalrec_all = rbind(totalrec_all, TotalRecres)
    #   Mest = data.frame(Pred = exp(models$sd_rep$par.fixed[names(models$sd_rep$par.fixed) == "ln_M"]), True = M, Sex = rep(c("F", "M")), sim = sim)
    #   m_all = rbind(Mest, m_all)
    #   recest = data.frame(Pred = exp(models$sd_rep$par.fixed[names(models$sd_rep$par.fixed) == "RecPars"]), True = r0, sim = sim)
    #   rec_all = rbind(rec_all, recest)
    #   # fmort = data.frame(Pred = exp(models$sd_rep$par.fixed[names(models$sd_rep$par.fixed) == "ln_Fy"]), 
    #                      # Truł = Fmort[-n_years,1,sim], years = 1:length(models$rep$Total_Biom), sim = sim)
    #   # fmort_all = rbind(fmort_all, fmort)
    #   
    #   # # Save selex estimates
    #   selex_f = data.frame(Pred = models$rep$Fish_Slx[1,,1,1], True = FishAge_Selex[,1,1], sim = sim, sex = "F", age = age_bins)
    #   # selex_m = data.frame(Pred = models$rep$Fish_Slx[1,,2,1], True = FishAge_Selex[,2,1], sim = sim, sex = "M", age = age_bins)
    #   # selex_all = rbind(selex_f, selex_m, selex_all)
    #   
    #   fmsy_val = get_Fmsy(ln_Fmsy = log(0.1), M = Mest$Pred[1], selex = selex_f$Pred,
    #            waa = as.matrix(waa[,1]), mat_at_age = mat_at_age[,1], ages = age_bins)[[1]]
    #   fmsy_est = data.frame(Pred = fmsy_val, True = fmsy, sim = sim)
    #   fmsy_all <- rbind(fmsy_est, fmsy_all)l
    #   SBPR_MSY_est = get_SBPR(M = Mest$Pred[1], selex = selex_f$Pred, Trial_F = fmsy_val,
    #                       waa = as.matrix(waa[,1]), mat_at_age = mat_at_age[,1], ages = age_bins)$SBPR_sum
    #   Req_est = get_Req(SBPR_Fmsy = SBPR_MSY_est, waa = as.matrix(waa[,1]), mat_at_age = mat_at_age[,1], ages = age_bins)
    #   bmsy_est = SBPR_MSY_est * Req_est 
    #   Reqres = data.frame(Pred = Req_est, True = Req, sim = sim)
    #   Req_all = rbind(Reqres, Req_all)
    #   bmsy_df = data.frame(Pred = bmsy_est, True = bmsy, sim = sim)
    #   bmsy_all = rbind(bmsy_all, bmsy_df)
    #   
    #   est_hcr_catch[sim] =  get_proj_catch(fmsy_val = fmsy_val,
    #                                 bmsy_val = bmsy_est,
    #                                 sex_ratio = models$rep$init_sexRatios,
    #                                 n_ages = n_ages, n_sexes = em_inputs$data$n_sexes,
    #                                 term_NAA = models$rep$NAA[n_years-1,,],
    #                                 term_SSB = models$rep$SSB[n_years - 1],
    #                                 term_F_Slx = models$rep$Fish_Slx[n_years-1,,,],
    #                                 term_F = exp(models$sd_rep$par.fixed[names(models$sd_rep$par.fixed) == "ln_Fy"])[n_years - 1],
    #                                 M_s = exp(models$sd_rep$par.fixed[names(models$sd_rep$par.fixed) == "ln_M"]),
    #                                 r0 = exp(models$sd_rep$par.fixed[names(models$sd_rep$par.fixed) == "RecPars"]),
    #                                 WAA = as.matrix(waa[,1]), MatAA = mat_at_age)
    # 
    #     get_proj_catch(fmsy_val = fmsy, bmsy_val = bmsy, sex_ratio = 1, 
    #                  n_ages = n_ages, n_sexes = 1, term_NAA = rowSums(NAA[n_years-1,,,sim]), 
    #                  term_SSB = SSB[n_years-1, sim], term_F_Slx = FishAge_Selex, r0 = r0,
    #                  term_F = Fmort[n_years - 1,,sim], M_s = M, WAA = waa, MatAA = mat_at_age)
    #   
    #   median((est_hcr_catch[1:sim] - HCR_proj_catch[1:sim]) / HCR_proj_catch[1:sim], na.rm = T)
    #   
    #  ratio_par = models$sd_rep$par.fixed[names(models$sd_rep$par.fixed) == "logit_init_sexRatio"]
    #   # sexratio[sim] =  0 + (1 - 0) * (1 / (1 + exp(-ratio_par)))
    # }
    print(sim)
  } # end sim

m_all_df = m_all %>% mutate(RE = (Pred - True) / True, type = paste("M", Sex, sep = "_")) %>% select(-Sex)
rec_all_df = rec_all %>% mutate(RE = (Pred - True) / True, type = "R0")
fmsy_all_df = fmsy_all %>% mutate(RE = (Pred - True) / True, type = "Fmsy")
bmsy_all_df = bmsy_all %>% mutate(RE = (Pred - True) / True, type = "Bmsy")
Req_all_df = Req_all %>% mutate(RE = (Pred - True) / True, type = "Req")
par_all = rbind(Req_all_df, bmsy_all_df, m_all_df, rec_all_df, fmsy_all_df) %>% 
  select(RE, type, sim) %>%  pivot_wider(names_from = "type", values_from = "RE") %>%
  select(-sim)
median(rec_all_df$RE,na.rm = T)
median(m_all_df$RE,na.rm = T)
median(bmsy_all_df$RE)

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

selex_all %>% 
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
  geom_hline(yintercept = 0) 
  # coord_cartesian(ylim = c(-0.75, 0.75)) 
  # ggthemes::scale_color_excel_new() +
  # ggthemes::theme_excel_new()
  
par(mar=c(4,6,1,2))

nf <- layout(
  matrix(c(1,2,3,4,rep(5, 4)), ncol=4, byrow=TRUE), 
  widths=c(1,1), 
  heights=c(2,2)
)
# plot(NAA[1,,1,sim-1], type = "l")
# lines(models$rep$NAA[1,,1])

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

pars = read.csv("/Users/matthewcheng/Desktop/UAF Thesis/Sex_Strc_Sim/output/Experiment 1/Growth_M (0,0)/Age/Parameters.csv")
a = pars %>% 
  filter(Type == "Bmsy") %>% 
  mutate(RE = (Pred - Truth) /Truth)
median(a$RE)
