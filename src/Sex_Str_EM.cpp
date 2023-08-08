// General single species age-and sex-structured stock assessment
// that allows for fitting to length comps, and accommodates multiple
// fishery fleets and survey fleets.
// Sex index follows females then males
// Creator: Matthew LH. Cheng (UAF-CFOS)
// Date updated: 8/6/23

#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density; // Define namespace to use SEPARABLE, AR1, SCALE
  using namespace Eigen; // Define namespace for Eigen functions (i.e., sparse matrix)
  
  // DATA SECTION ----------------------------
  // Model Dimensions ------------------------
  DATA_VECTOR(years); // Vector of years
  DATA_VECTOR(ages); // Vector of age bins
  DATA_VECTOR(lens); // Vector of length bins
  DATA_INTEGER(n_sexes); // Number of sexes
  DATA_INTEGER(n_fish_fleets); // Number of fishery fleets
  DATA_INTEGER(n_srv_fleets); // Number of survey fleets
  
  int n_years = years.size(); // Number of years
  int n_ages = ages.size(); // Number of age bins
  int n_lens = lens.size(); // Number of length bins
  
  // Observations ----------------------------
  // Fishery Observations --------------------
  DATA_MATRIX(obs_catch); // Array of catch from each fleet; n_years x n_fish_fleets 
  DATA_VECTOR(catch_cv); // Vector of CVs for catch; n_fish_fleets
  DATA_MATRIX(obs_fish_index); // Matrix of fishery indices; n_years x n_fish_fleets
  DATA_VECTOR(fish_index_cv); // Vector of CVs for fishery index; n_fish_fleets
  DATA_ARRAY(obs_fish_age_comps); // Array of fishery observed age comps; n_years x n_ages x n_sexes x n_fish_fleets 
  DATA_ARRAY(fish_age_comps_inputN); // Array of input sample sizes for fishery age comps; n_years x n_sexes x n_fish_fleets 
  DATA_ARRAY(obs_fish_len_comps); // Array of fishery observed length comps; n_years x n_lens x n_sexes x n_fish_fleets 
  DATA_ARRAY(fish_len_comps_inputN); // Array of input sample sizes for fishery length comps; n_years x n_sexes x n_fish_fleets 
  
  // Survey Observations -------------------
  DATA_MATRIX(obs_srv_index); // Matrix of fishery indices; n_years x n_srv_fleets
  DATA_VECTOR(srv_index_cv); // Vector of CVs for survey index; n_srv_fleets
  DATA_ARRAY(obs_srv_age_comps); // Array of observed survey age comps; n_years x n_ages x n_sexes x n_srv_fleets 
  DATA_ARRAY(srv_age_comps_inputN); // Array of input sample sizes survey survey age comps; n_years x n_sexes x n_srv_fleets 
  DATA_ARRAY(obs_srv_len_comps); // Array of observed survey length comps; n_years x n_lens x n_sexes x n_srv_fleets 
  DATA_ARRAY(srv_len_comps_inputN); // Array of input sample sizes for survey length comps; n_years x n_sexes x n_srv_fleets 
  
  // Biological Observations (Fixed) ----------------
  DATA_VECTOR(sexRatio); // Vector of sex ratio; n_sexes (Females, Males)
  DATA_MATRIX(WAA); // Array for weight-at-age relationship (n_ages, n_sexes)
  DATA_VECTOR(MatAA); // Vector of maturity-at-age
  DATA_ARRAY(age_len_transition); // Array for age length transition matrix (n_ages, n_lens, n_sexes)
  DATA_MATRIX(age_len_transition_unsexed); // Matrix for age length transition matrix (n_ages, n_lens)
  
  // Data Indicators -------------------------------
  // Indicator 0 == not fitting, 1 == fit
  DATA_IMATRIX(use_catch); // Using catch data; n_years, n_fish_fleets
  DATA_IMATRIX(use_fish_index); // Using fishery index data; n_years, n_fish_fleets
  DATA_IMATRIX(use_srv_index); // Using survey index data; n_years, n_srv_fleets
  DATA_IMATRIX(use_fish_age_comps); // Using fishery age comp data; n_years, n_fish_fleets
  DATA_IMATRIX(use_fish_len_comps); // Using fishery len comp data; n_years, n_fish_fleets
  DATA_IMATRIX(use_srv_age_comps); // Using survey age comp data; n_years, n_srv_fleets
  DATA_IMATRIX(use_srv_len_comps); // Using survey len comp data; n_years, n_srv_fleets
  DATA_INTEGER(p_ow_sex_fish_age); // Fishery Proportions over sexes, or within lengths; 0 = within sex; 1 = over sex
  DATA_INTEGER(p_ow_sex_fish_len); // Fishery Proportions over sexes, or within lengths; 0 = within sex; 1 = over sex
  DATA_INTEGER(p_ow_sex_srv_age); // Survey Proportions over sexes, or within lengths; 0 = within sex; 1 = over sex
  DATA_INTEGER(p_ow_sex_srv_len); // Survey Proportions over sexes, or within lengths; 0 = within sex; 1 = over sex
  
  // PARAMETER SECTION ---------------------
  // Biological Parameters -----------------
  PARAMETER_VECTOR(ln_M); // Log Natural Mortality; n_sexes
  PARAMETER_VECTOR(ln_InitDevs); // Log Deviations from initial equilibrium age-structure; n_ages - 1
  PARAMETER_VECTOR(ln_RecDevs); // Log Deviations from recruitment model; n_years - 1
  PARAMETER_VECTOR(RecPars); // Log Recruitment parameters; 2 (log_R0 and h for Beverton-Holt)
  PARAMETER(ln_sigmaRec); // Log sigma for recruitment
  
  // Index of Abundance Parameters ---------
  PARAMETER_VECTOR(ln_q_fish); // Vector of fishery catchability coefficients; n_fish_fleets
  PARAMETER_VECTOR(ln_q_srv); // Vector of survey catchability coefficients; n_srv_fleets
  
  // Selectivity and Removal Parameters -----
  PARAMETER_MATRIX(ln_Fy); // Fishing Mortality Parameters; n_years x n_fish_fleets
  PARAMETER_ARRAY(ln_fish_selpars); // Fishery Selectivity Parameters; n_sexes x n_fish_fleets x n_fish_selpars (2; ln_a50 and ln_k)
  PARAMETER_ARRAY(ln_srv_selpars); // Fishery Selectivity Parameters; n_sexes x n_fish_fleets x n_srv_selpars (2; ln_a50 and ln_k)
  
  // PARAMETER TRANSFORMATIONS --------------------
  vector<Type> M = exp(ln_M); // Natural Mortality in normal Space
  vector<Type> InitDevs = exp(ln_InitDevs); // Deviations from Initial equilibrium age-structure in normal space
  vector<Type> RecDevs = exp(ln_RecDevs); // Deviations from Initial equilibrium age-structure in normal space
  Type sigmaRec2 = pow(exp(ln_sigmaRec), 2); // Recruitment Variability in normal Space
  vector<Type> q_fish = exp(ln_q_fish); // Fishery catchability in normal space
  vector<Type> q_srv = exp(ln_q_srv); // Survey catchability in normal space
  
  // PREDICTED QUANTITES & CONTAINERS -----------
  // Predicted Quantities ------------------------
  matrix<Type> pred_catch(n_years, n_fish_fleets); // Predicted fishery catches
  array<Type> pred_fish_age_comps(obs_fish_age_comps.dim); // Predicted Fishery Age Comps
  array<Type> pred_fish_len_comps(obs_fish_len_comps.dim); // Predicted Fishery Length Comps
  array<Type> pred_srv_age_comps(obs_srv_age_comps.dim); // Predicted survey age comps
  array<Type> pred_srv_len_comps(obs_srv_len_comps.dim); // Predicted survey length comps
  matrix<Type> pred_fish_index(n_years, n_fish_fleets); // Predicted fishery indices
  matrix<Type> pred_srv_index(n_years, n_srv_fleets); // Predicted survey indices
  
  // Storage Containers -------------------------
  array<Type> NAA(n_years, n_ages, n_sexes); // Numbers at age
  array<Type> NAL(n_years, n_lens, n_sexes); // Numbers at length
  array<Type> CAA(n_years, n_ages, n_sexes, n_fish_fleets); // Catch-at-age
  array<Type> CAL(n_years, n_lens, n_sexes, n_fish_fleets); // Catch-at-length
  array<Type> FAA(n_years, n_ages, n_sexes, n_fish_fleets); // Fishing mortality-at-age
  array<Type> Total_FAA(n_years, n_ages, n_sexes); // Total Fishing mortality-at-age
  array<Type> ZAA(n_years, n_ages, n_sexes); // Total mortality-at-age
  array<Type> SAA(n_years, n_ages, n_sexes); // Survival at age
  vector<Type> SSB(n_years); // Spawning stock biomass
  vector<Type> SBPR_N(n_ages); // Numbers Per recruit container
  vector<Type> SBPR_SSB0(n_ages); // Spawning Biomass Per recruit container
  vector<Type> Total_Rec(n_years); // Total Recruitment
  vector<Type> Total_Biom(n_years); // Total Biomass
  array<Type> Fish_Slx(n_years, n_ages, n_sexes, n_fish_fleets); // Fishery Selectivities
  array<Type> Srv_Slx(n_years, n_ages, n_sexes, n_srv_fleets); // Survey Selectivities
  array<Type> Total_Fish_Age_Numbers(n_years, n_sexes, n_fish_fleets); // Store Total Fishery Numbers for Proportions
  array<Type> Total_Fish_Len_Numbers(n_years, n_sexes, n_fish_fleets); // Store Total Survey Length Numbers for Proportions
  array<Type> Total_Srv_Age_Numbers(n_years, n_sexes, n_srv_fleets); // Store Total Survey Age Numbers for Proportions
  array<Type> Total_Srv_Len_Numbers(n_years, n_sexes, n_srv_fleets); // Store Total Survey Length Numbers for Proportions
  
  // INITIALIZE LIKELIHOODS ---------------------
  matrix<Type> catch_nLL(n_years, n_fish_fleets); // Catch likelihood
  matrix<Type> fish_index_nLL(n_years, n_fish_fleets); // Fishery index likelihood
  matrix<Type> srv_index_nLL(n_years, n_srv_fleets);  // Survey index likelihood
  array<Type> fish_age_comp_nLL(n_years, n_sexes, n_fish_fleets); // Fishery age comps likelihood
  array<Type> srv_age_comp_nLL(n_years, n_sexes, n_srv_fleets); // Survey age comps likelihood
  array<Type> fish_len_comp_nLL(n_years, n_sexes, n_fish_fleets); // Fishery length comps likelihood
  array<Type> srv_len_comp_nLL(n_years, n_sexes, n_srv_fleets); // Survey length comps likelihood
  Type rec_nLL = 0; // Recruitment likelihood penalty 
  Type sexRatio_nLL = 0; // Sex-Ratio likelihood penalty
  Type jnLL = 0; // Joint Negative log Likelihood
  
  // Set containers to zeros
  SSB.setZero();
  NAA.setZero();
  NAL.setZero();
  CAA.setZero();
  CAL.setZero();
  pred_catch.setZero();
  pred_fish_index.setZero();
  pred_srv_index.setZero();
  pred_fish_age_comps.setZero();
  pred_fish_len_comps.setZero();
  pred_srv_age_comps.setZero();
  pred_srv_len_comps.setZero();
  Total_Biom.setZero();
  Total_Rec.setZero();
  
  // Set likelihood components to zeros
  catch_nLL.setZero();
  fish_index_nLL.setZero();
  srv_index_nLL.setZero();
  fish_age_comp_nLL.setZero();
  srv_age_comp_nLL.setZero();
  fish_len_comp_nLL.setZero();
  srv_len_comp_nLL.setZero();
  
  // MODEL SECTION -------------------------
  
  // Compute Deaths and Removals -----------
  for(int y = 0; y < n_years; y++) {
    for(int a = 0; a < n_ages; a++) {
      for(int s = 0; s < n_sexes; s++) {
        for(int f = 0; f < n_fish_fleets; f++) {
          
          // Compute Fishery Selectivity 
          Type fish_slope = exp(ln_fish_selpars(s,f,0)); 
          Type fish_a50 = exp(ln_fish_selpars(s,f,1)); 
          Fish_Slx(y,a,s,f) = Type(1) / (Type(1) + exp(-fish_slope * (ages(a) - fish_a50)));
          
          // Calculate fishing mortality-at-age
          FAA(y,a,s,f) = exp(ln_Fy(y,f)) * Fish_Slx(y,a,s,f);
          Total_FAA(y,a,s) += FAA(y,a,s,f);
          
        } // end fish fleet loop
        
        // Calculate survival and total mortality-at-age
        ZAA(y,a,s) = Total_FAA(y,a,s) + M(s); // Total Mortality
        SAA(y,a,s) = exp(Type(-1.0) * ZAA(y,a,s)); // Survival
        
      } // end sex loop
    } // end age loop
  } // end year loop
  
  // Initialize the population -------------
  for(int s = 0; s < n_sexes; s++) {
    for(int a = 0; a < n_ages; a++) {
      if(a < n_ages - 1) { // not plus-group
        NAA(0, a, s) = exp(RecPars(0)) * exp(-M(s) * Type(a)) * InitDevs(a) * sexRatio(s); 
      } else{ // plus group
        NAA(0,n_ages - 1,s) = (exp(RecPars(0)) * exp(-M(s) * Type(a))) / (1 - exp(-M(s))) * sexRatio(s); 
      } 
      // Calculate starting SSB
      if(s == 0) SSB(0) += NAA(0, a, 0) * MatAA(a) * WAA(a,0); 
      // Compute Numbers-at-length in year 1
      vector<Type> NAL_tmp_vec = age_len_transition.col(s).transpose().matrix() * NAA.col(s).transpose().col(0).matrix(); 
      for(int l = 0; l < n_lens; l++) NAL(0,l,s) = NAL_tmp_vec(l); // Loop through to input
    } // end age loop
  } // end sex loop
  
  // SSB0 Calculations ---------------------
  // Loop through spawning biomass per recruit calculations
  for(int a = 0; a < n_ages; a++) {
    if(a == 0) SBPR_N(a) = Type(1);
    if(a > 0 && a < n_ages - 1) SBPR_N(a) = SBPR_N(a - 1) * exp(-M(0)); 
    if(a == n_ages - 1) SBPR_N(a) = ((SBPR_N(a - 1) * exp(-M(0))) / (1 - exp(-M(0))));
    SBPR_SSB0(a) = SBPR_N(a) * MatAA(a) * WAA(a, 0); 
  } // age loop
  
  // Get Virgin Equilibrium SSB here
  Type ssb0 = SBPR_SSB0.sum() * exp(RecPars(0)); // SSB0
  
  // Project population forward -------------
  for(int y = 1; y < n_years; y++){
    for(int s = 0; s < n_sexes; s++) {
      
      // Recruitment
      Type R0 = exp(RecPars(0)); // Virgin Recruitment
      Type h = RecPars(1); // Steepness in normal space
      // Define BH components
      Type BH_first_part = Type(4) * h * R0 * SSB(y-1);
      Type BH_sec_part = (ssb0 * (Type(1) - h)) + SSB(y-1) * ((Type(5)*h) - 1);
      Type BH_tmp_rec = BH_first_part/BH_sec_part; // Deterministic BH
      // Get recruitment with process error here
      NAA(y, 0, s) =  BH_tmp_rec * RecDevs(y-1) * sexRatio(s);
      Total_Rec(y) += NAA(y, 0, s); // Get total recruitment
      
      for(int a = 1; a < n_ages; a++) {
        // Project ages forward
        if(a < n_ages-1) {
          NAA(y,a,s) = NAA(y-1,a-1,s) * SAA(y-1,a-1,s);
        } // not plus group
        if(a == n_ages-1)  {
          NAA(y,n_ages-1,s) = (NAA(y-1,n_ages-1,s) * SAA(y-1,n_ages-1,s)) + // plus group decrement
                              (NAA(y-1,n_ages-2,s) * SAA(y-1,n_ages-2,s)); // increment into plus-group
        } // plus group
        
        if(s == 0) SSB(y) += NAA(y,a,0) * MatAA(a) * WAA(a,s); // Calculate SSB here (Indexing 0 for females)
      } // end age loop
      
      // Compute Numbers-at-length 
      vector<Type> NAL_tmp_vec = age_len_transition.col(s).transpose().matrix() * NAA.col(s).transpose().col(y).matrix(); 
      for(int l = 0; l < n_lens; l++) NAL(y,l,s) = NAL_tmp_vec(l); // Loop through to input
      
    } // end sex loop
  } // end year loop
  
  // Calculate Catch Quantities -----------
  for(int f = 0; f < n_fish_fleets; f++) {
    for(int y = 0; y < n_years; y++) {
      for(int s = 0; s < n_sexes; s++) {
        for(int a = 0; a < n_ages; a++) {
          // Baranov's Catch Equation - catch at age
          CAA(y,a,s,f) = NAA(y,a,s) * (Type(1.0) - exp(-ZAA(y,a,s))) * 
                        (FAA(y,a,s,f) / ZAA(y,a,s));
          // Get aggregated catch 
          pred_catch(y,f) += CAA(y,a,s,f) * WAA(a,s);
          
        } // age loop
        // Compute catch-at-length
        vector<Type> CAL_tmp_vec = age_len_transition.col(s).transpose().matrix() * 
                                   CAA.col(f).col(s).transpose().col(y).matrix(); 
        for(int l = 0; l < n_lens; l++) CAL(y,l,s,f) = CAL_tmp_vec(l); // Loop through to input
      } // sex loop
    } // year loop
  } // end fishery fleet
  
  // Calculate Fishery Index ----------------------------
  for(int f = 0; f < n_fish_fleets; f++) {
    for(int y = 0; y < n_years; y++) { 
      for(int s = 0; s < n_sexes; s++) {
        for(int a = 0; a < n_ages; a++) {
          // Increment to get total index, conditional on fishery selectivity
          pred_fish_index(y,f) += NAA(y,a,s) * WAA(a,s) * Fish_Slx(y,a,s,f);
        } // a loop
        // Apply catchability
        pred_fish_index(y,f) *= exp(ln_q_fish(f)); 
      } // s loop
    } // y loop
  } // f loop 
  
  // Calculate Survey Index ----------------------------
  for(int sf = 0; sf < n_srv_fleets; sf++) {
    for(int y = 0; y < n_years; y++) {
      for(int s = 0; s < n_sexes; s++) {
        
        // Extract parameters for survey selectivity
        Type srv_slope = exp(ln_srv_selpars(s,sf,0)); 
        Type srv_a50 = exp(ln_srv_selpars(s,sf,1)); 
        
        for(int a = 0; a < n_ages; a++) {
          // Calculate Survey Selectivity
          Srv_Slx(y,a,s,sf) = Type(1) / (Type(1) + exp(-srv_slope * (ages(a) - srv_a50)));
          // Increment to get total index, conditional on survey selectivity
          pred_srv_index(y,sf) += NAA(y,a,s) * Srv_Slx(y,a,s,sf);
        } // a loop
        // Apply catchability
        pred_srv_index(y,sf) *= exp(ln_q_srv(sf));
        
      } // s loop
    } // y loop
  } // sf loop
  
  
  // Calculate Compositions and Related Quantities ------
  // Fishery Compositions 
  for(int y = 0; y < n_years; y++) {
    for(int f = 0; f < n_fish_fleets; f++) {
      for(int s = 0; s < n_sexes; s++) {
        
        // Computing Age Compositions
        for(int a = 0; a < n_ages; a++) {
          // Get predicted comps here prior to normalizing w/ catch at age
          pred_fish_age_comps(y,a,s,f) = CAA(y,a,s,f);
          // Increment to get total numbers at age for a given fleet
          Total_Fish_Age_Numbers(y,s,f) += pred_fish_age_comps(y,a,s,f);
          
          // Normalize to sum to 1 (Proportions within a given sex)
          if(a == n_ages - 1 && p_ow_sex_fish_age == 0) {
            for(int a = 0; a < n_ages; a++) {
              pred_fish_age_comps(y,a,s,f) /= Total_Fish_Age_Numbers(y,s,f);
            } // a loop
          } // if a = plus group and proportion within sexes
          
          // Normalize to sum to 1 (Proportions across sexes)
          if(a == n_ages - 1 && s == n_sexes - 1 && p_ow_sex_fish_age == 1) {
            // Get total fishery age numbers across sexes
            Type Total_fishAge_across_sex_tmp = sum(Total_Fish_Age_Numbers.col(f).transpose().col(y));
            for(int s = 0; s < n_sexes; s++) {
              for(int a = 0; a < n_ages; a++) {
                pred_fish_age_comps(y,a,s,f) /= Total_fishAge_across_sex_tmp;
              } // a loop
            } // s loop
          } // if a = plus group, s = n_sexes, and proportions across sexes
        } // a loop
        
        // Computing Length Compositions
        for(int l = 0; l < n_lens; l++) {
          // Get predicted comps here for catch-at-length (prior to normalization)
          pred_fish_len_comps(y,l,s,f) = CAL(y,l,s,f);
          // Get total len numbers for normalizing
          Total_Fish_Len_Numbers(y,s,f) += pred_fish_len_comps(y,l,s,f);
          
          // Normalize to sum to 1 (Proportions within sexes)
          if(l == n_lens - 1 && p_ow_sex_fish_len == 1) {
            for(int l = 0; l < n_lens; l++) {
              pred_fish_len_comps(y,l,s,f) /= Total_Fish_Len_Numbers(y,s,f);
            } // l loop
          } // if l = n_lens - 1 and proportions within sexes
          
          // Normalize to sum to 1 (Proportions across sexes)
          if(l == n_lens - 1 && s == n_sexes - 1 && p_ow_sex_fish_len == 1) {
            // Extract total fishery lengths across exes
            Type Total_fishLen_across_sex_tmp = sum(Total_Fish_Len_Numbers.col(f).transpose().col(y));
            for(int s = 0; s < n_sexes; s++) {
              for(int l = 0; l < n_lens; l ++) {
                pred_fish_len_comps(y,l,s,f) /= Total_fishLen_across_sex_tmp;
              } // l loop
            } // s loop
          } // if l == n_lens - 1, s == n_sexes - 1, and proportions across sexes
          
        } // l loop
        
      } // s loop
    } // f loop
  } // y loop
  
  // Survey Compositions 
  for(int y = 0; y < n_years; y++) {
    for(int sf = 0; sf < n_srv_fleets; sf++) {
      for(int s = 0; s < n_sexes; s++) {
        
        // Computing Age Compositions
        for(int a = 0; a < n_ages; a++) {
          // Get predicted comps here prior to normalizing w/ catch at age
          pred_srv_age_comps(y,a,s,sf) = NAA(y,a,s) * Srv_Slx(y,a,s,sf);
          // Increment to get total numbers at age for a given fleet
          Total_Srv_Age_Numbers(y,s,sf) += pred_srv_age_comps(y,a,s,sf);
          // Normalize to sum to 1
          if(a == n_ages - 1) {
            for(int a = 0; a < n_ages; a++) {
              pred_srv_age_comps(y,a,s,sf) /= Total_Srv_Age_Numbers(y,s,sf);
            } // a loop
          } // if a = plus group
        } // a loop
        
        // Computing Length Compositions
        vector<Type> srv_lens_tmp = age_len_transition.col(s).transpose().matrix() * 
          pred_srv_age_comps.col(sf).col(s).transpose().col(y).matrix(); 
        for(int l = 0; l < n_lens; l++) pred_srv_len_comps(y,l,s,sf) = srv_lens_tmp(l); // Loop through to input
        
      } // s loop
    } // sf loop
  } // y loop
  
  // LIKELIHOOD SECTION --------------------
  // Catch Likelihood ----------------------
  vector<Type> catch_sd(n_fish_fleets); // Empty container to store sd calculation
  // Catch observed w/ minimal error
  for(int f = 0; f < n_fish_fleets; f++) {
    // Calculate sd
    catch_sd(f) = sqrt(log( (pow(catch_cv(f),2) + 1)));
    for(int y = 0; y < n_years; y ++) {
      // Get likelihood here
      catch_nLL(y, f) -= use_catch(y,f) * dnorm(log(obs_catch(y, f)), 
                         log(pred_catch(y, f)) - pow(catch_sd(f), 2)/2, 
                         catch_sd(f), true);
    } // fish fleet loop
  } // year loop
  
  // Fishery Index Likelihoods -----------
  vector<Type> fish_index_sd(n_fish_fleets); // Empty container to store sd calculation
  for(int f = 0; f < n_fish_fleets; f++) {
    // Calculate sd
    fish_index_sd(f) = sqrt(log( (pow(fish_index_cv(f),2) + 1)));
    for(int y = 0 ; y < n_years; y++) {
      // Get likelihood here
      fish_index_nLL(y,f) -= use_fish_index(y,f) * dnorm(log(obs_fish_index(y, f)), 
                             log(pred_fish_index(y, f)) - pow(fish_index_sd(f),2)/2, 
                             fish_index_sd(f), true);
    } // y loop
  } // f loop
  
  // Survey Index Likelihoods -----------
  vector<Type> srv_index_sd(n_srv_fleets); // Empty container to store sd calculation
  for(int sf = 0; sf < n_srv_fleets; sf++) {
    // Calculate sd
    srv_index_sd(sf) = sqrt(log( (pow(srv_index_cv(sf),2) + 1)));
    for(int y = 0 ; y < n_years; y++) {
      // Get likelihood here
      srv_index_nLL(y,sf) -= use_srv_index(y,sf) * dnorm(log(obs_srv_index(y, sf)),
                             log(pred_srv_index(y, sf)) - pow(srv_index_sd(sf),2)/2,
                             srv_index_sd(sf), true);
    } // y loop
  } // sf loop
  
  // Fishery Composition Likelihoods -----
  for(int f = 0; f < n_fish_fleets; f++) {
    for(int y = 0 ; y < n_years; y++) {
      
    // Proportions across sexes for fishery ages or lengths
    if(p_ow_sex_fish_age == 1 || p_ow_sex_fish_len == 1) {
      
      // Age-Compositions
      if(p_ow_sex_fish_age == 1) {
        vector<Type> obs_fish_age_as = obs_fish_age_comps.col(f).transpose().col(y).transpose().vec(); // Pull out observed vector
        vector<Type> pred_fish_age_as = pred_fish_age_comps.col(f).transpose().col(y).transpose().vec(); // Pull out observed vector
        fish_age_comp_nLL(y,0,f) -= use_fish_age_comps(y,f) * dmultinom(obs_fish_age_as, pred_fish_age_as, true); // Get likelihood here
      } // if proportions across sex for fishery ages
      
      // Length-Compositions
      if(p_ow_sex_fish_age == 1) {
        vector<Type> obs_fish_len_as = obs_fish_len_comps.col(f).transpose().col(y).transpose().vec(); // Pull out observed vector
        vector<Type> pred_fish_len_as = pred_fish_len_comps.col(f).transpose().col(y).transpose().vec(); // Pull out observed vector
        fish_age_comp_nLL(y,0,f) -= use_fish_len_comps(y,f) * dmultinom(obs_fish_len_as, pred_fish_len_as, true); // Get likelihood here
      } // if proportions across sex for fishery lengths
        
    } // if either fishery age or fishery lengths are proportions across sexes
      
      // Proportions within sexes for fishery ages or lengths
      if(p_ow_sex_fish_age == 0 || p_ow_sex_fish_len == 0) {
        
        for(int s = 0; s < n_sexes; s++) {
          // Age-Compositions
          if(p_ow_sex_fish_age == 0) {
            vector<Type> obs_fish_age_ws = obs_fish_age_comps.col(f).col(s).transpose().col(y); // Pull out observed vector
            vector<Type> pred_fish_age_ws = pred_fish_age_comps.col(f).col(s).transpose().col(y); // Pull out predicted vector
            fish_age_comp_nLL(y,s,f) -= use_fish_age_comps(y,f) * dmultinom(obs_fish_age_ws, pred_fish_age_ws, true); // Get likelihood here
          } // if proportions within sex for fishery ages
          
          // Length-Compositions
          if(p_ow_sex_fish_len == 0) {
            vector<Type> obs_fish_len_ws = obs_fish_len_comps.col(f).col(s).transpose().col(y); // Pull out observed vector
            vector<Type> pred_fish_len_ws = pred_fish_len_comps.col(f).col(s).transpose().col(y); // Pull out predicted vector
            fish_len_comp_nLL(y,s,f) -= use_fish_len_comps(y,f) * dmultinom(obs_fish_len_ws, pred_fish_len_ws, true);
          } // if proportions within sex for fishery lengths
          
        } // s loop
      } // if either fishery age or fishery lengths are proportions within sexes
      
    } // y loop
  } // f loop
  
  // Survey Composition Likelihoods ------
  for(int sf = 0; sf < n_srv_fleets; sf++) {
    for(int y = 0 ; y < n_years; y++) {
      for(int s = 0; s < n_sexes; s++) {
        
        // Pre-processing - extract out quantities
        vector<Type> obs_srv_age = obs_srv_age_comps.col(sf).col(s).transpose().col(y); // Pull out observed vector
        vector<Type> pred_srv_age = pred_srv_age_comps.col(sf).col(s).transpose().col(y); // Pull out predicted vector
        vector<Type> obs_srv_len = obs_srv_len_comps.col(sf).col(s).transpose().col(y); // Pull out observed vector
        vector<Type> pred_srv_len = pred_srv_len_comps.col(sf).col(s).transpose().col(y); // Pull out predicted vector
        
        // Get likelihood here
        srv_age_comp_nLL(y,s,sf) -= use_srv_age_comps(y,sf) * dmultinom(obs_srv_age, pred_srv_age, true);
        srv_len_comp_nLL(y,s,sf) -= use_srv_len_comps(y,sf) * dmultinom(obs_srv_len, pred_srv_len, true);
        
      } // s loop
    } // y loop
  } // sf loop
  
  // Penalty for Population Initialization
  for(int a = 0; a < ln_InitDevs.size(); a ++) {
    rec_nLL -= dnorm(ln_InitDevs(a), -sigmaRec2/Type(2), exp(ln_sigmaRec));
  } 
  
  // Penalty for recruitment deviates
  for(int a = 0; a < ln_RecDevs.size(); a ++) {
    rec_nLL -= dnorm(ln_RecDevs(a), -sigmaRec2/Type(2), exp(ln_sigmaRec));
  } 
  
  // Compute joint likelihood
  jnLL = rec_nLL + sum(catch_nLL) + sum(fish_index_nLL) + sum(fish_age_comp_nLL) +
    sum(fish_len_comp_nLL) + sum(srv_index_nLL) + sum(srv_age_comp_nLL) +
    sum(srv_len_comp_nLL);
  
  // REPORT SECTION ------------------------
  REPORT(NAA);
  REPORT(NAL);
  REPORT(CAA);
  REPORT(CAL);
  REPORT(SSB);
  REPORT(pred_catch);
  REPORT(pred_fish_index);
  REPORT(pred_fish_age_comps);
  REPORT(pred_fish_len_comps);
  REPORT(pred_srv_index);
  REPORT(pred_srv_age_comps);
  REPORT(pred_srv_len_comps);
  REPORT(Fish_Slx);
  REPORT(Srv_Slx);
  
  REPORT(jnLL);
  REPORT(rec_nLL);
  REPORT(catch_nLL);
  REPORT(fish_index_nLL);
  REPORT(fish_age_comp_nLL);
  REPORT(fish_len_comp_nLL);
  REPORT(srv_index_nLL);
  REPORT(srv_age_comp_nLL);
  REPORT(srv_len_comp_nLL);
  
  // DERIVED QUANTITIES -------------------
  ADREPORT(SSB);
  
  return(jnLL);
  
} // end objective function
