# Purpose: To prepare inputs for sex-structured EMs
# Author: Matthew LH. Cheng (UAF-CFOS)
# Date 8/24/25

#' Title
#'
#' @param catch_cv Vector of catch CVs 
#' @param WAA Vector of weight-at-age
#' @param age_len_transition Array of age length transition (a,l,s)
#' @param fish_age_prop Treatment of proportions for fishery age (within sex and across sex)
#' @param fish_len_prop Treatment of proportions for fishery lens (within sex and across sex)
#' @param srv_age_prop Treatment of proportions for survey age (within sex and across sex)
#' @param srv_len_prop Treatment of proportions for survey lens (within sex and across sex)
#' @param agg_fish_age Whether or not to aggregate fishery ages
#' @param agg_srv_age Whether or not to aggregate survey ages
#' @param use_catch TRUE or FALSE, for using catch data
#' @param use_fish_index TRUE or FALSE, for using fishery index
#' @param use_srv_index TRUE or FALSE, for using survey index
#' @param use_fish_age_comps TRUE or FALSE, for using fishery age comps data
#' @param use_fish_len_comps TRUE or FALSE, for using fishery len comps data
#' @param use_srv_age_comps TRUE or FALSE, for using survey age comps data
#' @param use_srv_len_comps TRUE or FALSE, for using survey len comps data
#' @param fix_pars Vector of characters to fix parameters
#' @param share_M_sex TRUE OR else... for sharing Ms across sex
#' @param sex_specific Whether or not this is a sex-specific assessment (TRUE OR FALSE)
#' @param n_sexes Number of sexes to model
#' @param sim simulation index
#' @param sexRatio Sex Ratio
#'
#' @return
#' @export
#'
#' @examples
prepare_EM_inputs = function(sim,
                             catch_cv,
                             WAA,
                             sexRatio,
                             n_sexes,
                             age_len_transition,
                             fish_age_prop,
                             fish_len_prop,
                             srv_age_prop,
                             srv_len_prop,
                             agg_fish_age,
                             agg_srv_age,
                             use_catch = TRUE,
                             use_fish_index = TRUE,
                             use_srv_index = TRUE,
                             use_fish_age_comps = TRUE,
                             use_fish_len_comps = TRUE,
                             use_srv_age_comps = TRUE,
                             use_srv_len_comps = TRUE,
                             fix_pars,
                             share_M_sex,
                             sex_specific){
  
  if(sex_specific == TRUE & n_sexes == 1) stop("Sex-specific assessment, but n_sexes = 1")
  if(sex_specific == FALSE & n_sexes == 2) stop("Sex-aggregated assessment, but n_sexes = 2")
  
  # List objects to output
  data = list()
  parameters = list()
  map = list()

# Controls ----------------------------------------------------------------
  data$years = 1:(n_years-1)  
  data$ages = 1:n_ages
  data$lens = 1:n_lens
  data$n_sexes = n_sexes
  data$n_fish_fleets = n_fish_fleets
  data$n_srv_fleets = n_srv_fleets

# Fishery Data Inputs -----------------------------------------------------
  data$obs_catch = as.matrix(Total_Catch[-n_years,,sim], dim = c(n_years-1, n_fish_fleets))
  data$obs_fish_index = array(Fish_Index[-n_years,,sim], dim = c(n_years-1, n_fish_fleets))
  data$fish_age_comps_inputN = array(Fish_Neff_Age[-n_years,sim], dim = c(n_years-1, n_fish_fleets))
  data$fish_len_comps_inputN = array(Fish_Neff_Len[-n_years,sim], dim = c(n_years-1, n_fish_fleets))
  data$catch_cv = catch_cv
  data$fish_index_cv = cv_Fish_Index
  
  # Sex-specific Assessment with sex-specific comps
  if(sex_specific == TRUE & n_sexes > 1) data$obs_fish_age_comps = array(Fish_AgeComps[-n_years,,,,sim], dim = c(n_years-1, n_ages, n_sexes, n_fish_fleets))
  if(sex_specific == TRUE & n_sexes > 1) data$obs_fish_len_comps = array(Fish_LenComps[-n_years,,,,sim], dim = c(n_years-1, n_lens, n_sexes, n_fish_fleets))
  
  # Sex-aggregated assessment with sex-aggregated comps
  if(sex_specific == FALSE & n_sexes == 1) data$obs_fish_age_comps = array(apply(Fish_AgeComps[-n_years,,,,sim], c(1, 2, 4), sum), 
                                                                           dim = c(n_years-1, n_ages, n_sexes, n_fish_fleets))
  if(sex_specific == FALSE & n_sexes == 1) data$obs_fish_len_comps = array(apply(Fish_LenComps[-n_years,,,,sim], c(1, 2, 4), sum), 
                                                                           dim = c(n_years-1, n_lens, n_sexes, n_fish_fleets))

# Survey Data Inputs ------------------------------------------------------
  data$obs_srv_index = array(Srv_Index[-n_years,,sim], dim = c(n_years-1, n_srv_fleets))
  data$srv_age_comps_inputN = array(Srv_Neff_Age[-n_years,sim], dim = c(n_years-1, n_srv_fleets))
  data$srv_len_comps_inputN = array(Srv_Neff_Len[-n_years,sim], dim = c(n_years-1, n_srv_fleets))
  data$srv_index_cv = cv_Srv_Index
  
  # Sex-Specific Asessment with sex-specific comps
  if(sex_specific == TRUE & n_sexes > 1) data$obs_srv_age_comps = array(Srv_AgeComps[-n_years,,,,sim], dim = c(n_years-1, n_ages, n_sexes, n_srv_fleets))
  if(sex_specific == TRUE & n_sexes > 1) data$obs_srv_len_comps = array(Srv_LenComps[-n_years,,,,sim], dim = c(n_years-1, n_lens, n_sexes, n_fish_fleets))
  
  # Sex-aggregated assessment with sex-aggregated comps
  if(sex_specific == FALSE & n_sexes == 1) data$obs_srv_age_comps = array(apply(Srv_AgeComps[-n_years,,,,sim], c(1, 2, 4), sum), 
                                                                           dim = c(n_years-1, n_ages, n_sexes, n_srv_fleets))
  if(sex_specific == FALSE & n_sexes == 1) data$obs_srv_len_comps = array(apply(Srv_LenComps[-n_years,,,,sim], c(1, 2, 4), sum), 
                                                                           dim = c(n_years-1, n_lens, n_sexes, n_srv_fleets))

# Biological Inputs -------------------------------------------------------
  data$sexRatio = sexRatio
  data$MatAA = as.vector(mat_at_age[1,,1,1])
  data$WAA = WAA
  data$age_len_transition = age_len_transition

# Data Indicators ---------------------------------------------------------
  
  # Fishery
  # fishery catch
  if(use_catch == TRUE) data$use_catch = matrix(1, ncol = n_fish_fleets, nrow = n_years - 1)
  else data$use_catch = matrix(0, ncol = n_fish_fleets, nrow = n_years - 1)
  # fishery index
  if(use_fish_index == TRUE) data$use_fish_index = matrix(1, ncol = n_fish_fleets, nrow = n_years - 1)
  else data$use_fish_index = matrix(0, ncol = n_fish_fleets, nrow = n_years - 1)
  # fishery age comps
  if(use_fish_age_comps == TRUE) data$use_fish_age_comps = matrix(1, ncol = n_fish_fleets, nrow = n_years - 1)
  else data$use_fish_age_comps = matrix(0, ncol = n_fish_fleets, nrow = n_years - 1)
  # fishery length comps
  if(use_fish_len_comps == TRUE) data$use_fish_len_comps = matrix(1, ncol = n_fish_fleets, nrow = n_years - 1)
  else data$use_fish_len_comps = matrix(0, ncol = n_fish_fleets, nrow = n_years - 1)

  if(fish_age_prop == "within") data$p_ow_sex_fish_age = 0 # within sexes
  if(fish_len_prop == "within") data$p_ow_sex_fish_len = 0 # within sexes
  if(fish_age_prop == "across") data$p_ow_sex_fish_age = 1 # across sexes
  if(fish_len_prop == "across") data$p_ow_sex_fish_len = 1 # across sexes
  if(agg_fish_age == FALSE) data$agg_sex_fish_age = 0 # don't aggregate fish age comps
  if(agg_fish_age == TRUE) data$agg_sex_fish_age = 1 # aggregate fish age comps
  
  # Survey
  # survey index
  if(use_srv_index == TRUE) data$use_srv_index = matrix(1, ncol = n_srv_fleets, nrow = n_years - 1)
  else data$use_srv_index = matrix(0, ncol = n_srv_fleets, nrow = n_years - 1)
  # survey age comps
  if(use_srv_age_comps == TRUE)  data$use_srv_age_comps = matrix(1, ncol = n_srv_fleets, nrow = n_years - 1)
  else data$use_srv_age_comps = matrix(0, ncol = n_srv_fleets, nrow = n_years - 1)
  # survey length comps
  if(use_srv_len_comps == TRUE) data$use_srv_len_comps = matrix(1, ncol = n_srv_fleets, nrow = n_years - 1)
  else data$use_srv_len_comps = matrix(0, ncol = n_srv_fleets, nrow = n_years - 1)

  if(srv_age_prop == "within") data$p_ow_sex_srv_age = 0 # within sexes
  if(srv_len_prop == "within") data$p_ow_sex_srv_len = 0 # within sexes
  if(srv_age_prop == "across") data$p_ow_sex_srv_age = 1 # across sexes
  if(srv_len_prop == "across") data$p_ow_sex_srv_len = 1 # across sexes
  if(agg_srv_age == FALSE) data$agg_sex_srv_age = 0 # don't aggregate survey age comps
  if(agg_srv_age == TRUE) data$agg_sex_srv_age = 1 # aggregate survey age comps
  
# Parameters --------------------------------------------------------------
  parameters$ln_M = vector(length = n_sexes)
  for(s in 1:n_sexes) parameters$ln_M[s] = log(M[s])
  parameters$ln_InitDevs = log(InitDevs[,sim])
  parameters$ln_RecDevs = log(RecDevs[-50,sim])
  parameters$RecPars = c(log(r0),h)
  parameters$ln_sigmaRec = log(sigma_rec)
  parameters$ln_q_fish = log(q_Fish)
  parameters$ln_q_srv = log(q_Srv)
  parameters$ln_Fy = matrix(log(Fmort[-n_years,,sim]), ncol = n_fish_fleets)
  parameters$ln_fish_selpars = array(log(3), dim = c(n_sexes, n_fish_fleets, 2)) # only for logistic (last dim = number of selex pars)
  parameters$ln_srv_selpars = array(log(3), dim = c(n_sexes, n_srv_fleets, 2))  # only for logistic (last dim = number of selex pars)
  
# Mapping -----------------------------------------------------------------
  # fixing steepness
  if(sum(fix_pars %in% c("h")) == 1) {
    map$RecPars <- factor(c(1, NA))
    # Remove steepness from fix_pars vector so it goes through the next loop properly
    fix_pars <- fix_pars[fix_pars != 'h']
  } # end if for steepness
  
  # if we want to share a single M between sexes
  if(share_M_sex == TRUE & sex_specific == TRUE & n_sexes > 1) {
    map$ln_M = factor(1, 1)
  } # end if for sharing M between sexes
  
  # Loop through to map parameters
  for(i in 1:length(fix_pars)) {
    # Get parameter length here
    par_length <- length(unlist(parameters[names(parameters) == fix_pars[i]]))
    # Now, stick the map parameter into a list
    map_par <- list( factor(rep(NA, par_length)) )
    names(map_par) <- fix_pars[i] # name the list
    # Now, append this to our map list
    map <- c(map_par, map)
  } # end i

  return(list(data = data, parameters = parameters, map = map))
} # end function
