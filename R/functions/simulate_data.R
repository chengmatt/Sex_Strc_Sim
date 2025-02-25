# Purpose: Function to simulate data for a sex-and age-structured population with ability to sample
# sex-specific age and length comps for multiple fisheries and surveys
# Date: 8/27/23
# Creator: Matthew LH. Cheng (UAF-CFOS)

simulate_data = function(spreadsheet_path,
                         Fish_Neff_Age = 100,
                         Fish_Neff_Len = 100,
                         Srv_Neff_Age = 100,
                         Srv_Neff_Len = 100,
                         F_pattern = "Contrast",
                         comp_across_sex = "within",
                         q_Fish = 0.025,
                         cv_Fish_Index = 0.25,
                         q_Srv = 0.05,
                         selex_type,
                         cv_Srv_Index = 0.25,
                         catch_cv = 0.025,
                         growth_control = NA,
                         natmort_control = NA,
                         growth_control_fct = 1,
                         natmort_control_fct = 1,
                         force_grwth_same_yng = FALSE,
                         sexRatio
                         ) {
  
  # Load in all functions from the functions folder
  fxn_path <- here("R", "functions")
  files <- list.files(fxn_path)
  for(i in 1:length(files)) source(here(fxn_path, files[i]))
  
  # Read in parameters
  read_params(spreadsheet_path)
  
  # Create objects for storage
  NAA = array(0, dim = c(n_years, n_ages, n_sexes, n_sims)) # Numbers at age
  NAL = array(0, dim = c(n_years, length(len_bins), n_sexes, n_sims)) # Numbers at length
  FAA = array(0, dim = c(n_years, n_ages, n_sexes, n_fish_fleets, n_sims)) # Fishery Mortality at Age
  ZAA = array(0, dim = c(n_years, n_ages, n_sexes, n_sims)) # Total Mortality at Age
  CAA = array(0, dim = c(n_years, n_ages, n_sexes, n_fish_fleets, n_sims)) # Catch at Age
  CAL = array(0, dim = c(n_years, length(len_bins), n_sexes, n_fish_fleets, n_sims)) # Catch at Length
  Total_Catch_Sex = array(0, dim = c(n_years, n_sexes, n_fish_fleets, n_sims)) # Sex-Specific Catch
  Total_Catch = array(0, dim = c(n_years, n_fish_fleets, n_sims)) # Aggregated Catch
  SSB = array(0, dim = c(n_years, n_sims)) # Spawning stock biomass
  Total_Biom = array(0, dim = c(n_years, n_sims)) # Total Biomass
  RecDevs = array(0, dim = c(n_years-1, n_sims)) # recruitment deviates
  InitDevs = array(0, dim = c(n_ages, n_sims)) # initial deviates
  
  # Fishery Containers
  Fish_AgeComps = array(0, dim = c(n_years, n_ages, n_sexes, n_fish_fleets, n_sims)) 
  Fish_Neff_Age = array(Fish_Neff_Age, dim = c(n_years, n_fish_fleets)) # Age Effective Sample Size
  Fish_LenComps = array(0, dim = c(n_years, length(len_bins), n_sexes, n_fish_fleets, n_sims)) 
  Fish_Neff_Len= array(Fish_Neff_Len, dim = c(n_years, n_fish_fleets)) # Length Effective Sample Size
  Fish_Index = array(0, dim = c(n_years, n_fish_fleets, n_sims))
  FishAge_Selex = array(0, dim = c(n_ages, n_sexes, n_fish_fleets))
  Fmort = array(0, dim = c(n_years, n_fish_fleets, n_sims)) # fishing mortality container
  HCR_proj_catch = vector()

  # Survey Containers
  Srv_AgeComps = array(0, dim = c(n_years, n_ages, n_sexes, n_srv_fleets, n_sims)) 
  Srv_Neff_Age = array(Srv_Neff_Age, dim = c(n_years, n_srv_fleets)) # Age Effective Sample Size
  Srv_LenComps = array(0, dim = c(n_years, length(len_bins), n_sexes, n_srv_fleets, n_sims)) 
  Srv_Neff_Len= array(Srv_Neff_Len, dim = c(n_years, n_srv_fleets)) # Length Effective Sample Size
  Srv_Index = array(0, dim = c(n_years, n_srv_fleets, n_sims))
  SrvAge_Selex = array(0, dim = c(n_ages, n_sexes, n_srv_fleets))
  Srv_LAA = vector("list", length = n_sexes * (n_years-1) * n_srv_fleets) # pre-allocate vector list size
  Srv_LW = vector("list", length = n_sexes * (n_years-1) * n_srv_fleets) # pre-allocate vector list size
  
  # if statement for if we want to vary linf and k relative to females
  if(growth_control == "chg_males_rel_females") { 
    k[2] = k[1] * growth_control_fct
    L_inf[2]= L_inf[1] * growth_control_fct
  }
  
  # changing natural mortality males relative to females
  if(natmort_control == "chg_males_rel_females")  M[2] = M[1] * natmort_control_fct

  # counters for sampling survey length-at-age and weight-length
  srv_counter_age = 1
  srv_counter_len = 1
  
  # Across or within sexes for sampling compositional data
  if(comp_across_sex == "across") comp_across_sex = 0
  if(comp_across_sex == "within") comp_across_sex = 1
  
# Construct Length, Weight, and Age-Length Matrix Relationships -----------
  
  # Construct Weight-Length relationship
  wl_female = alpha_wl[1] * len_bins^beta_wl[1]
  wl_male = alpha_wl[1] * len_bins^beta_wl[1]
  wl <<- matrix(c(wl_female, wl_male), ncol = n_sexes)
  
  # Construct vonB LAA 
  vonB_Female = vonB_est(age_bins = age_bins, k = k[1], L_inf = L_inf[1], t0 = t0[1], sd = 0) # females
  vonB_Male = vonB_est(age_bins = age_bins, k = k[2], L_inf = L_inf[2], t0 = t0[2], sd = 0) # males
  
  if(force_grwth_same_yng == TRUE) { # to minimze t0
    
    # set up
    t0_trial = seq(-10, -3, 0.001) # define trial values to search through
    minimize = vector()

    for(i in 1:length(t0_trial)) {
      vonB_Male = vonB_est(age_bins = age_bins, k = k[2], L_inf = L_inf[2], t0 = t0_trial[i], sd = 0) # males
      minimize[i] = sum(vonB_Female[1] - vonB_Male[1])^2 # minimize using ssq
    } # end i loop
    
    t0[2] = t0_trial[which.min(minimize)] # redefine the k for males 
    vonB_Male = vonB_est(age_bins = age_bins, k = k[2], L_inf = L_inf[2], t0 = t0[2], sd = 0) # males - now redfine this and bind
    
  } # end if forcing young individuals to similar growth 
  
  # Bind females and males
  vonB <<- matrix(c(vonB_Female, vonB_Male), ncol = n_sexes)
  
  # Get WAA relationship
  winf_Female = alpha_wl[1] * L_inf[1]^beta_wl[1]
  winf_Male = alpha_wl[2] * L_inf[2]^beta_wl[2]
  waa_Female = winf_Female * (1 - exp(-k[1] * (age_bins - t0[1])))^beta_wl[1]
  waa_Male = winf_Male * (1 - exp(-k[2] * (age_bins - t0[2])))^beta_wl[2]
  waa <<- matrix(c(waa_Female, waa_Male), ncol = n_sexes)
  
  # Get vonB sd with linear interpolation
  vonBf_sd_vec = vonB_sd1[1] + (((vonB_Female - vonB_Female[1]) / (L_inf[1] - vonB_Female[1])) * (vonB_sd2[1] - vonB_sd1[1]))
  
  # make sd constant between sex for males if no diff in growth
  if(growth_control_fct == 1) vonBm_sd_vec = vonB_sd1[1] + (((vonB_Male - vonB_Male[1]) / (L_inf[1] - vonB_Male[1])) * (vonB_sd2[1] - vonB_sd1[1]))
  if(growth_control_fct != 1) vonBm_sd_vec = vonB_sd1[2] + (((vonB_Male - vonB_Male[1]) / (L_inf[2] - vonB_Male[1])) * (vonB_sd2[2] - vonB_sd1[2]))
  
  # Bind together to make a matrix
  vonB_sd = matrix(cbind(vonBf_sd_vec, vonBm_sd_vec), ncol = n_sexes) # females then males
  
  # Construct age-length transition matrices 
  al_matrix_Female <<- get_al_trans_matrix(age_bins = age_bins, len_bins = len_bins,
                                         mean_length = vonB_Female, sd = vonBf_sd_vec) # Female
  al_matrix_Male <<- get_al_trans_matrix(age_bins = age_bins, len_bins = len_bins, 
                                       mean_length = vonB_Male, sd = vonBm_sd_vec) # Male
  al_matrix <<- array(c(al_matrix_Female, al_matrix_Male), 
                    dim = c(length(age_bins), length(len_bins), n_sexes)) # combine age-length matrices
  

# Selectivity Parameterization --------------------------------------------

  # Construct selectivity age-based and length-based
  if(selex_type == "length") {
    # fishery selex
    fish_len_selex = matrix(logist(slope = fish_len_slope,
                                   bins = rep(len_bins,n_fish_fleets), midpoint = fish_len_midpoint),
                            nrow = length(len_bins), ncol = n_fish_fleets)
    for(f in 1:n_fish_fleets) { # convert length to age
      for(s in 1:n_sexes) { 
        FishAge_Selex[,s,f] = al_matrix[,,s] %*% fish_len_selex[,f] 
      } # end s loop for sexes
    } # end f loop for fishery fleets
    
    # survey selex
    srv_len_selex = matrix(logist(slope = srv_len_slope,
                                  bins = rep(len_bins,n_fish_fleets), midpoint = srv_len_midpoint),
                           nrow = length(len_bins), ncol = n_srv_fleets)
    for(sf in 1:n_srv_fleets) { # convert length to age
      for(s in 1:n_sexes) { # standardize to 1
        SrvAge_Selex[,s,sf] = al_matrix[,,s] %*% srv_len_selex[,sf]
      } # end s loop for sexes
    } # end f loop for survey fleets
  } # end length based selectivity
  
  if(selex_type == "age") { # age-based selectivity
    fish_age_selex_Female = logist(slope = fish_age_slope_f, bins = age_bins, midpoint = fish_age_midpoint_f)
    fish_age_selex_Male = logist(slope = fish_age_slope_m, bins = age_bins, midpoint = fish_age_midpoint_m)
    srv_age_selex_Female = logist(slope = srv_age_slope_f, bins = age_bins, midpoint = srv_age_midpoint_f)
    srv_age_selex_Male = logist(slope = srv_age_slope_m, bins = age_bins, midpoint = srv_age_midpoint_m)
    FishAge_Selex = array(c(fish_age_selex_Female, fish_age_selex_Male),dim = c(length(age_bins), n_sexes, n_fish_fleets))
    SrvAge_Selex = array(c(srv_age_selex_Female, srv_age_selex_Male),dim = c(length(age_bins), n_sexes, n_srv_fleets))
  } # end if for age-based selectivity 
  
  
  # Get reference points + set up Fs ----------------------------------------
  fmsy = get_Fmsy(ln_Fmsy = log(0.1),  M = M[1],  selex = FishAge_Selex[,1,1], 
                       waa = waa[,1], mat_at_age = mat_at_age[,1], ages = age_bins)[[1]]
  # get bmsy
  SBPR_MSY = get_SBPR(M = M[1], selex = FishAge_Selex[,1,1], Trial_F = fmsy, 
                      waa = waa[,1], mat_at_age = mat_at_age[,1], ages = age_bins)$SBPR_sum
  Req = get_Req(SBPR_Fmsy = SBPR_MSY, waa = waa[,1], mat_at_age = mat_at_age[,1], ages = age_bins)
  bmsy = SBPR_MSY * Req
  
  # Set up fishing mortality
  Fmort[] = f_pattern_scenarios(F_pattern = "Contrast", n_years = n_years, fmsy = fmsy) # Specify fishing mortality scenarios

# Start Simulation --------------------------------------------------------

  for(sim in 1:n_sims) {
    
    # Generate recruitment deviates and deviations from equilibrium
    RecDevs[,sim] = exp(rnorm(n_years-1, mean = -sigma_rec^2/2, sd = sigma_rec))
    InitDevs[,sim] = exp(rnorm(length(age_bins), mean = -sigma_rec^2/2, sd = sigma_rec))

    for(s in 1:n_sexes) {
      # Get numbers at age - not plus group
      NAA[1,,s,sim] = r0 * exp(-M[s] * 0:(length(age_bins)-1)) * InitDevs[,sim] * sexRatio[s]
      # Convert NAA to numbers at length
      NAL[1,,s,sim] = t(al_matrix[,,s]) %*% NAA[1,,s,1]
    } # end first sex loop
    
    # Calculate SSB at time t = 1
    SSB[1,sim] = sum(NAA[1,,1,sim] * waa[,1] * mat_at_age[,1])
    Total_Biom[1, sim] = sum(NAA[1,,,sim] * waa[,]) # get total biomass at time t1

# Population Projection ---------------------------------------------------

    for(y in 2:n_years) {
      for(a in 1:n_ages) {
        for(s in 1:n_sexes) {
          
          # Calculate Deaths from Fishery
          FAA[y - 1,a,s,,sim] = Fmort[y-1,,sim] * FishAge_Selex[a,s,] # Fishing Mortality at Age
          # Calculate Total Mortality
          ZAA[y - 1,a,s,sim] = M[s] + sum(FAA[y - 1,a,s,,sim])
          
          # Project Population Forward ----------------------------------------------
          
          # Recruitment
          if(a == 1) {
            NAA[y,1,s,sim] = beverton_holt_recruit(ssb = SSB[y - 1, sim],  h = h, 
                                                   r0 = r0, M = M[1]) * RecDevs[y - 1,sim] * sexRatio[s]
          } # end if recruitment age
          
          if(a > 1 & a < n_ages) {
            NAA[y,a,s,sim] = NAA[y - 1,a - 1,s,sim] * exp(-ZAA[y - 1,a - 1,s,sim])
          } # end if between recruitment age and plus group
          
          # Calculate abundance for plus group (Mortality of MaxAge-1 + Mortality of MaxAge)
          if(a == n_ages) { 
            NAA[y,a,s,sim] = (NAA[y - 1,a - 1,s,sim] * exp(-ZAA[y - 1,a - 1,s,sim])) + 
              (NAA[y - 1,a,s,sim] * exp(-ZAA[y - 1,a,s,sim]))
            
            # Convert NAA to numbers at length
            NAL[y,,s,sim] = t(al_matrix[,,s]) %*% NAA[y,,s,1]
          } # end if for plus group
        } # end second sex loop
      } # end age loop
      
      # Calculate SSB here
      SSB[y,sim] = sum(NAA[y,,1,sim] * waa[,1] * mat_at_age[,1])
      Total_Biom[y,sim] = sum(NAA[y,,,sim] * waa[,]) # get total biomass here
      
      # Observation Model (Fishery) -------------------------------------------------------
      for(f in 1:n_fish_fleets) {
        for(s in 1:n_sexes) {
          # Get catch at age
          CAA[y-1,,s,f,sim] = (FAA[y-1,,s,f,sim] / ZAA[y - 1,,s,sim]) * 
            (NAA[y-1,,s,sim] * (1 - exp(-ZAA[y - 1,,s,sim])))
          # Get Catch at Length
          CAL[y-1,,s,f,sim] = t(al_matrix[,,s]) %*% CAA[y-1,,s,f,sim]
          # Get Total Catch by Sex
          Total_Catch_Sex[y-1,s,f,sim] = CAA[y-1,,s,f,sim] %*% waa[,s]
        } # end third sex loop

          ### Simulate Fishery Compositions (Within sexes) ---------------------------------------
          if(comp_across_sex == 1) {
            
            # Get sex-ratios (age-comps)
            p_af = sum(CAA[y-1,,1,f,sim]) / sum(CAA[y-1,,,f,sim])
            pa = rbinom(Fish_Neff_Age[y,f] * n_sexes, 1 ,p_af)
            Fish_AgeComps[y-1,,1,f,sim] = rmultinom(1, size = sum(pa),  prob = CAA[y-1,,1,f,sim])
            Fish_AgeComps[y-1,,2,f,sim] = rmultinom(1, size = (Fish_Neff_Age[y,f] * n_sexes) - sum(pa), prob = CAA[y-1,,2,f,sim])
            
            # Fishery Length Compositions
            p_lf = sum(CAL[y-1,,1,f,sim]) / sum(CAL[y-1,,,f,sim])
            pl = rbinom(Fish_Neff_Len[y,f] * n_sexes, 1 ,p_lf)
            Fish_LenComps[y-1,,1,f,sim] = rmultinom(1, size = sum(pl),  prob = CAL[y-1,,1,f,sim])
            Fish_LenComps[y-1,,2,f,sim] = rmultinom(1, size = (Fish_Neff_Len[y,f] * n_sexes) - sum(pl), prob = CAL[y-1,,2,f,sim])

          } # if fishery comps are simulated within sexes

        ### Simulate Fishery Compositions (Across sexes) ---------------------------------------
        if(comp_across_sex == 0) {
          Prob_FishAge = as.vector(CAA[y-1,,,f,sim]/sum(CAA[y-1,,,f,sim]))
          Fish_AgeComps[y-1,,,f,sim] = matrix(rmultinom(1, size = Fish_Neff_Age[y,f] * n_sexes, 
                                                        prob = Prob_FishAge), nrow = n_ages, ncol = n_sexes)
          
          # Fishery Length Compositions
          # Get Probability of sampling length bins
          Prob_FishLen = as.vector(CAL[y-1,,,f,sim]/sum(CAL[y-1,,,f,sim]))
          Fish_LenComps[y-1,,,f,sim] = matrix(rmultinom(1, size = Fish_Neff_Len[y,f] * n_sexes, 
                                                        Prob_FishLen), nrow = length(len_bins), ncol = n_sexes)
        } # if fishery comps are simulated across sexes
        
        ### Simulate Fishery Index --------------------------------------------------
        # Calculate sd for deviations
        Fish_Index_sd = sqrt(log(cv_Fish_Index[f]^2+1))
        # Sample Fishery Index with lognormal error
        Fish_Index[y-1,f,sim] = q_Fish[f] * sum(FishAge_Selex[,,f] * NAA[y-1,,,sim] * waa) * 
          exp(rnorm(1, -Fish_Index_sd^2/2, Fish_Index_sd))

        # Sum up sex specification catches and then add error
        Total_Catch[y-1,f,sim] = sum(Total_Catch_Sex[y-1,,f,sim])
      } # end fish fleet loop
      
      # Observation Model (Survey) ----------------------------------------------
      for(sf in 1:n_srv_fleets) {
        ### Simulate Survey Compositions (Within Sexes) ---------------------------------------
        if(comp_across_sex == 1) {
          
            # Survey Age Compositions
            p_af = sum(NAA[y-1,,1,sim] * SrvAge_Selex[,1,sf]) / sum((NAA[y-1,,,sim] * SrvAge_Selex[,,sf])) # Get probability of sampling ages
            pa = rbinom(Srv_Neff_Age[y,sf] * n_sexes, 1 ,p_af)
            Srv_AgeComps[y-1,,1,sf,sim] = rmultinom(1, size = sum(pa), NAA[y-1,,1,sim] * SrvAge_Selex[,1,sf])
            Srv_AgeComps[y-1,,2,sf,sim] = rmultinom(1, size = (Srv_Neff_Age[y,sf] * n_sexes) - sum(pa), NAA[y-1,,2,sim] * SrvAge_Selex[,2,sf])
            
            # Survey Length Compositions
            # Get probability of sampling lengths
            Prob_SrvLen_F = (t(al_matrix[,,1]) %*% (NAA[y-1,,1,sim] * SrvAge_Selex[,1,sf])) 
            Prob_SrvLen_M = (t(al_matrix[,,2]) %*% (NAA[y-1,,2,sim] * SrvAge_Selex[,2,sf])) 
            p_lf = sum(Prob_SrvLen_F) / sum(Prob_SrvLen_M, Prob_SrvLen_F)
            pl = rbinom(Srv_Neff_Len[y,sf] * n_sexes, 1 , p_lf)
            Srv_LenComps[y-1,,1,sf,sim] = rmultinom(1, size = sum(pl), Prob_SrvLen_F)
            Srv_LenComps[y-1,,2,sf,sim] = rmultinom(1, size = (Srv_Neff_Len[y,sf] * n_sexes) - sum(pl), Prob_SrvLen_M)
            
          } # if survey comps are simulated within sexes

        ### Simulate Survey Compositions (Across Sexes) ---------------------------------------
        if(comp_across_sex == 0) {
          # Survey Age Compositions
          Prob_SrvAge = (NAA[y-1,,,sim] * SrvAge_Selex[,,sf]) / sum(NAA[y-1,,,sim] * SrvAge_Selex[,,sf]) # Get probability of sampling ages
          Srv_AgeComps[y-1,,,sf,sim] = matrix(rmultinom(1, size = Srv_Neff_Age[y,sf] * n_sexes,  # sampling 
                                                        as.vector(Prob_SrvAge)), nrow = n_ages, ncol = n_sexes)
          
          # Survey Length Compositions
          # Get probability of sampling lengths
          Prob_SrvLen = vector() # re-initialize vector
          for(s in 1:n_sexes) {
            Prob_SrvLen_s = (t(al_matrix[,,s]) %*% Prob_SrvAge[,s]) 
            Prob_SrvLen = rbind(Prob_SrvLen, Prob_SrvLen_s)
          } # end fifth sex loop
          
          Prob_Len = Prob_SrvLen / sum(Prob_SrvLen) # compute prob of len sampling
          Srv_LenComps[y-1,,,sf,sim] = matrix(rmultinom(1, size = Srv_Neff_Len[y,sf] * n_sexes, # sampling
                                                        Prob_SrvLen), nrow = length(len_bins), ncol = n_sexes)
        } # if survey comps are simulated across sexes
        
        ### Simulate Survey Index --------------------------------------------------
        # Calculate sd for deviations
        Srv_Index_sd = sqrt(log(cv_Srv_Index[sf]^2+1))
        # Sample Fishery Index with lognormal error
        Srv_Index[y-1,sf,sim] = q_Srv[sf] * sum(SrvAge_Selex[,,sf] * NAA[y-1,,,sim]) * 
          exp(rnorm(1, -Srv_Index_sd^2/2, Srv_Index_sd))
        
        ### Simulate Growth from Survey Age and Length Compositions ---------------------------
        for(s in 1:n_sexes) {
          
          # Extract composition values here
          ageCompValue = Srv_AgeComps[y-1,,s,sf,sim]
          non_zero_age = which(ageCompValue != 0)
          lenCompValue = Srv_LenComps[y-1,,s,sf,sim]
          non_zero_len = which(lenCompValue != 0)
          
          # Get survey length at age
          Srv_LAA_tmp <- lapply(non_zero_age, function(a) {
            n_age_samples <- ageCompValue[a]
            sampled_srv_lens <- rnorm(n = n_age_samples, mean = vonB[a, s], sd = as.numeric(vonB_sd[a,s]))
            data.frame(lens = sampled_srv_lens, ages = age_bins[a], sex = s, srv_fleet = sf, sim = sim)
          })
          
          # Get survey length weight relationship
          Srv_LW_tmp = lapply(non_zero_len, function(l) {
            n_len_samples <- lenCompValue[l]
            sampled_srv_wts = wl[l,s] * exp(rnorm(n_len_samples, -sqrt(log(wl_sd[s]^2 +1 ))^2/2, sqrt(log(wl_sd[s]^2 +1 )) ))
            data.frame(wts = sampled_srv_wts, lens = len_bins[l], sex = s, srv_fleet = sf, sim = sim)
          })
          
          # put laa and lw into list when done
          Srv_LAA_tmp <- do.call(rbind, Srv_LAA_tmp)
          Srv_LAA[[srv_counter_age]] = Srv_LAA_tmp
          Srv_LW_tmp <- do.call(rbind, Srv_LW_tmp)
          Srv_LW[[srv_counter_len]] = Srv_LW_tmp
          
          # update counters
          srv_counter_age = srv_counter_age + 1
          srv_counter_len = srv_counter_len + 1
          
          
        } # end fifth sex loop
      } # end survey fleet loop
      
      if(y == n_years - 1) {
        # Get projected catch here
        HCR_proj_catch[sim] = get_proj_catch(fmsy_val = fmsy, bmsy_val = bmsy, sex_ratio = sexRatio, 
                                             n_ages = n_ages, n_sexes = n_sexes, term_NAA = NAA[y,,,sim], 
                                             term_SSB = SSB[y, sim], term_F_Slx = FishAge_Selex, r0 = r0,
                                             term_F = Fmort[y,,sim], M_s = M, WAA = waa, MatAA = mat_at_age)
      } # get HCR catch
      
    } # end year loop
    
    print(sim)
  } # end sim loop
  
  # Output to environment
  NAA <<- NAA
  NAL <<- NAL
  FAA <<- FAA
  ZAA <<- ZAA
  CAA <<- CAA
  CAL <<- CAL
  Total_Catch_Sex <<- Total_Catch_Sex
  Total_Catch <<- Total_Catch
  SSB <<- SSB
  Total_Biom <<- Total_Biom
  RecDevs <<- RecDevs
  InitDevs <<- InitDevs
  
  # Fishery Containers
  Fish_AgeComps <<- Fish_AgeComps
  Fish_Neff_Age <<- Fish_Neff_Age
  Fish_LenComps <<- Fish_LenComps
  Fish_Neff_Len <<- Fish_Neff_Len
  Fish_Index <<- Fish_Index
  FishAge_Selex <<- FishAge_Selex
  
  # Survey Containers
  Srv_AgeComps <<- Srv_AgeComps
  Srv_Neff_Age <<- Srv_Neff_Age
  Srv_LenComps <<- Srv_LenComps
  Srv_Neff_Len <<- Srv_Neff_Len
  Srv_Index <<- Srv_Index
  SrvAge_Selex <<- SrvAge_Selex
  
  # Get Length weight samples
  Srv_LAA = data.table::rbindlist(Srv_LAA)
  Srv_LW = data.table::rbindlist(Srv_LW)
  Srv_LAA <<- Srv_LAA
  Srv_LW <<- Srv_LW

  # Biologicals and Selex and unceratinty
  M <<- M
  vonB <<- vonB
  sexRatio <<- sexRatio
  wl <<- wl
  waa <<- waa
  al_matrix_Female <<- al_matrix_Female
  al_matrix_Male <<- al_matrix_Male
  al_matrix <<- al_matrix
  FishAge_Selex <<- FishAge_Selex
  SrvAge_Selex <<- SrvAge_Selex
  Fmort <<- Fmort
  cv_Srv_Index <<- cv_Srv_Index
  cv_Fish_Index <<- cv_Fish_Index
  q_Fish <<- q_Fish
  q_Srv <<- q_Srv
  fmsy <<- fmsy
  bmsy <<- bmsy
  Req <<- Req 
  SBPR_MSY <<- SBPR_MSY
  HCR_proj_catch <<- HCR_proj_catch
  
  # Output these as a list
  oms <- list(
    n_fish_fleets = n_fish_fleets,
    n_srv_fleets = n_srv_fleets,
    len_bins = len_bins,
    n_ages = n_ages,
    n_sexes = n_sexes,
    age_bins = age_bins,
    len_bins = len_bins,
    NAA = NAA,
    NAL = NAL,
    FAA = FAA,
    ZAA = ZAA,
    CAA = CAA,
    CAL = CAL,
    Total_Catch_Sex = Total_Catch_Sex,
    Total_Catch = Total_Catch,
    SSB = SSB,
    Total_Biom = Total_Biom,
    RecDevs = RecDevs,
    InitDevs = InitDevs,
    mat_at_age = mat_at_age,
    Fish_AgeComps = Fish_AgeComps,
    Fish_Neff_Age = Fish_Neff_Age,
    Fish_LenComps = Fish_LenComps,
    Fish_Neff_Len = Fish_Neff_Len,
    Fish_Index = Fish_Index,
    FishAge_Selex = FishAge_Selex,
    Srv_AgeComps = Srv_AgeComps,
    Srv_Neff_Age = Srv_Neff_Age,
    Srv_LenComps = Srv_LenComps,
    Srv_Neff_Len = Srv_Neff_Len,
    Srv_Index = Srv_Index,
    SrvAge_Selex = SrvAge_Selex,
    Srv_LAA = Srv_LAA,
    Srv_LW = Srv_LW,
    M = M,
    r0 = r0,
    h = h,
    sigma_rec = sigma_rec,
    vonB = vonB,
    sexRatio = sexRatio,
    wl = wl,
    waa = waa,
    al_matrix_Female = al_matrix_Female,
    al_matrix_Male = al_matrix_Male,
    al_matrix = al_matrix,
    FishAge_Selex = FishAge_Selex,
    SrvAge_Selex = SrvAge_Selex,
    srv_len_selex = srv_len_selex,
    fish_len_selex = fish_len_selex,
    Fmort = Fmort,
    cv_Srv_Index = cv_Srv_Index,
    cv_Fish_Index = cv_Fish_Index,
    q_Fish = q_Fish,
    q_Srv = q_Srv,
    fmsy = fmsy,
    bmsy = bmsy,
    Req = Req,
    beta_wl = beta_wl,
    alpha_wl = alpha_wl,
    k = k,
    t0 = t0,
    n_sims = n_sims,
    n_years = n_years,
    L_inf = L_inf,
    vonBf_sd_vec = vonBf_sd_vec,
    vonBm_sd_vec = vonBm_sd_vec,
    SBPR_MSY = SBPR_MSY,
    HCR_proj_catch = HCR_proj_catch)
  
  return(oms)
} # end function