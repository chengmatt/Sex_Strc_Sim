#' Title Get estimated quantities
#'
#' @param model Model object
#' @param sim Simulation number
#' @param biologicals biologicals (growth)
#' @param om_name OM name
#' @param em_name EM name
#' @param n_sexes_em Number of modelled sexes
#'
#' @return
#' @export
#'
#' @examples
get_quantities = function(biologicals, model, sim, om_name, em_name, n_sexes_em) {
  
# Convergence -------------------------------------------------------------
gradient = max(model$sd_rep$gradient.fixed)
pdHess = model$sd_rep$pdHess
if(gradient < 0.01 & pdHess == TRUE) conv = "Converged"
else conv = "Not Converged"

# Put into dataframe to output
conv_df = data.frame(EM = em_name, OM = om_name, sim = sim, 
                     gradient = gradient, pdHess = pdHess,
                     convergence = conv)

# Get time series estimates -----------------------------------------------

# Total biomass
total_biomass_df = data.frame(Truth = Total_Biom[-n_years,sim], Pred = model$rep$Total_Biom,
                              Convergence = conv_df$convergence, sim = sim, EM = em_name, OM = om_name,
                              Type = "Total Biomass", Years = 1:length(model$rep$Total_Biom))
  
# SSB
total_ssb_df = data.frame(Truth = SSB[-n_years,sim], Pred = model$rep$SSB,
                              Convergence = conv_df$convergence, sim = sim, EM = em_name, OM = om_name,
                              Type = "Spawning Stock Biomass", Years = 1:length(model$rep$SSB))
# Total Recruitment
total_rec_df = data.frame(Truth = rowSums(NAA[-n_years,1,,sim]), Pred = model$rep$Total_Rec,
                          Convergence = conv_df$convergence, sim = sim, EM = em_name, OM = om_name,
                          Type = "Total Recruitment", Years = 1:length(model$rep$Total_Rec))

# Fishing Mortality
total_fmort_df = data.frame(Truth = Fmort[-n_years,1,sim], 
                          Pred = exp(model$sd_rep$par.fixed[names(model$sd_rep$par.fixed) == "ln_Fy"]),
                          Convergence = conv_df$convergence, sim = sim, EM = em_name, OM = om_name,
                          Type = "Total Fishing Mortality", Years = 1:length(Fmort[-n_years,1,sim]))

ts_df = rbind(total_biomass_df, total_ssb_df, total_rec_df, total_fmort_df)

# Sex Ratios
# Sum across 3rd dimension to get female sex ratio
true_NAA_female_sr = reshape2::melt( NAA[-n_years,,1,sim] / apply(NAA[-n_years,,,sim], c(1, 2), sum)) # true sex ratio
names(true_NAA_female_sr) = c("Years", "Age", "Truth")
pred_NAA_female_sr = reshape2::melt(model$rep$NAA[,,1] / apply(model$rep$NAA, c(1,2), sum)) # predicted sex ratio
names(pred_NAA_female_sr) = c("Years", "Age", "Pred")
NAA_sr_female_df = true_NAA_female_sr %>% left_join(pred_NAA_female_sr, by = c("Years", "Age")) %>% 
  mutate(Convergence = conv_df$convergence, sim = sim, EM = em_name, OM = om_name, Type = "NAA Female Sex Ratio")

# Biological Quantities + Selectivity -----------------------------------------------
if(n_sexes_em == 1) {
  # Get growth estimates
  pred_grwth = reshape2::melt(biologicals$waa_nosex)[-2]
  names(pred_grwth) = c("Age", "Pred")
  true_grwth = reshape2::melt(waa)
  names(true_grwth) = c("Age", "Sex", "True")
  grwth_df = true_grwth %>% left_join(pred_grwth, by = "Age") # left join 
  grwth_df = grwth_df %>% mutate(Convergence = conv_df$convergence, sim = sim, EM = em_name, OM = om_name,
                                 Type = "Weight at age")
  
  # Get fishery selectivity estimates
  fish_pred_selex = reshape2::melt(as.matrix(model$rep$Fish_Slx[1,,,1]))[-2]
  names(fish_pred_selex) = c("Age", "Pred")
  fish_true_selex = reshape2::melt(FishAge_Selex[,,1])
  names(fish_true_selex) = c("Age", "Sex", "True")
  fish_selex_df = fish_true_selex %>% left_join(fish_pred_selex, by = c("Age")) # left join 
  fish_selex_df = fish_selex_df %>% mutate(Convergence = conv_df$convergence, sim = sim, EM = em_name, OM = om_name,
                                 Type = "Fishery Selectivity")
  
  # get survey selectivity estimates
  srv_pred_selex = reshape2::melt(as.matrix(model$rep$Srv_Slx[1,,,1]))[-2]
  names(srv_pred_selex) = c("Age", "Pred")
  srv_true_selex = reshape2::melt(SrvAge_Selex[,,1])
  names(srv_true_selex) = c("Age", "Sex", "True")
  srv_selex_df = srv_true_selex %>% left_join(srv_pred_selex, by = c("Age")) # left join 
  srv_selex_df = srv_selex_df %>% mutate(Convergence = conv_df$convergence, sim = sim, EM = em_name, OM = om_name,
                                           Type = "Survey Selectivity")
  
  
} else{
  # Get growth estimates
  pred_grwth = reshape2::melt(biologicals$waa_sex)
  names(pred_grwth) = c("Age", "Sex", "Pred")
  true_grwth = reshape2::melt(waa)
  names(true_grwth) = c("Age", "Sex", "True")
  grwth_df = true_grwth %>% left_join(pred_grwth, by = c("Age", "Sex")) # left join 
  grwth_df = grwth_df %>% mutate(Convergence = conv_df$convergence, sim = sim, EM = em_name, OM = om_name,
                                 Type = "Weight at age")
  
  # Get fishery selectivity estimates
  fish_pred_selex = reshape2::melt(as.matrix(model$rep$Fish_Slx[1,,,1]))
  names(fish_pred_selex) = c("Age", "Sex", "Pred")
  fish_true_selex = reshape2::melt(FishAge_Selex[,,1])
  names(fish_true_selex) = c("Age", "Sex", "True")
  fish_selex_df = fish_true_selex %>% left_join(fish_pred_selex, by = c("Age", "Sex")) # left join 
  fish_selex_df = fish_selex_df %>% mutate(Convergence = conv_df$convergence, sim = sim, EM = em_name, OM = om_name,
                                           Type = "Fishery Selectivity")
  
  # get survey selectivity estimates
  srv_pred_selex = reshape2::melt(as.matrix(model$rep$Srv_Slx[1,,,1]))
  names(srv_pred_selex) = c("Age", "Sex", "Pred")
  srv_true_selex = reshape2::melt(SrvAge_Selex[,,1])
  names(srv_true_selex) = c("Age", "Sex", "True")
  srv_selex_df = srv_true_selex %>% left_join(srv_pred_selex, by = c("Age", "Sex")) # left join 
  srv_selex_df = srv_selex_df %>% mutate(Convergence = conv_df$convergence, sim = sim, EM = em_name, OM = om_name,
                                         Type = "Survey Selectivity")
} # end if else for binding selectivity and growth estimates froms single and multi-sex models

selex_all_df = rbind(srv_selex_df, fish_selex_df)

# Get reference points, HCR catch, and other parameters --------------------------------------

# Define selectivity for all sexes here
selex_all = as.matrix(model$rep$Fish_Slx[1,,,1])

# Define which waa to use
if(n_sexes_em == 1) waa_ref_use = biologicals$waa_nosex
if(n_sexes_em == 2) waa_ref_use = biologicals$waa_sex

# Get natural mortality
Mest = exp(model$sd_rep$par.fixed[names(model$sd_rep$par.fixed) == "ln_M"])
M_df = data.frame(Pred = Mest, Truth = M, Type = c("M_F", "M_M"), Convergence = conv_df$convergence, sim = sim, EM = em_name, OM = om_name)

# Get R0
# Meaning of R0 changes with and without sex
r0est = exp(model$sd_rep$par.fixed[names(model$sd_rep$par.fixed) == "RecPars"])
r0_df = data.frame(Pred = r0est, Truth = r0, Type = "R0", Convergence = conv_df$convergence, sim = sim, EM = em_name, OM = om_name)

# Get Fmsy
fmsy_est = get_Fmsy(ln_Fmsy = log(0.1), M = M_df$Pred[1], selex = selex_all[,1],
                    waa = waa_ref_use[,1], mat_at_age = mat_at_age[,1], ages = age_bins)[[1]]
fmsy_df = data.frame(Pred = fmsy_est, Truth = fmsy, Type = "Fmsy", Convergence = conv_df$convergence, sim = sim, EM = em_name, OM = om_name)

# Get bmsy
SBPR_MSY_est = get_SBPR(M = M_df$Pred[1], selex = selex_all[,1], Trial_F = fmsy_est, # first get sbpr msy
                        waa = waa_ref_use[,1], mat_at_age = mat_at_age[,1], ages = age_bins)$SBPR_sum
# Get req
Req_est = get_Req(SBPR_Fmsy = SBPR_MSY_est, waa = waa_ref_use[,1], mat_at_age = mat_at_age[,1], ages = age_bins)
bmsy_est = SBPR_MSY_est * Req_est # get bmsy here
bmsy_df = data.frame(Pred = bmsy_est, Truth = bmsy, Type = "Bmsy", Convergence = conv_df$convergence, sim = sim, EM = em_name, OM = om_name)

# Get initial sex ratio information
init_sr_est = model$rep$init_sexRatios
init_sr_df = data.frame(Pred = init_sr_est, Truth = sexRatio, Type = c("Female Sex Ratio", "Male Sex Ratio"), 
                        Convergence = conv_df$convergence, sim = sim, EM = em_name, OM = om_name)

par_df = rbind(init_sr_df, bmsy_df, fmsy_df, r0_df, M_df)

  return(
    list(ts_df = ts_df, NAA_sr_female_df = NAA_sr_female_df,
         grwth_df = grwth_df, selex_all_df = selex_all_df,
         conv_df = conv_df, par_df = par_df
         )
  )

} # end function
