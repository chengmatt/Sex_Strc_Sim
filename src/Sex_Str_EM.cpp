// General single species age-and sex-structured stock assessment
// that allows for fitting to length comps, and accommodates multiple
// fishery fleets and survey fleets.
// Sex index follows females then males
// Creator: Matthew LH. Cheng (UAF-CFOS)
// Date updated: 8/6/23

#include <TMB.hpp>
#include "Utility_Fxns.hpp"
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density; // Define namespace to use SEPARABLE, AR1, SCALE
  using namespace Eigen; // Define namespace for Eigen functions (i.e., sparse matrix)
  
  // DATA SECTION ----------------------------
  // Model Dimensions ------------------------
  DATA_VECTOR(years); // Vector of years
  DATA_VECTOR(ages); // Vector of age bins
  DATA_VECTOR(len_mids); // Vector of length bins (using midpoints)
  DATA_INTEGER(n_sexes); // Number of sexes
  DATA_INTEGER(n_fish_fleets); // Number of fishery fleets
  DATA_INTEGER(n_srv_fleets); // Number of survey fleets
  
  int n_years = years.size(); // Number of years
  int n_ages = ages.size(); // Number of age bins
  int n_lens = len_mids.size(); // Number of length bins
  
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
  DATA_INTEGER(agg_sex_fish_age); // Aggregating sexes for fishery ages; 0 == no aggregation, 1 == aggregate
  DATA_INTEGER(agg_sex_srv_age); // Aggregating sexes for survey ages; 0 == no aggregation, 1 == aggregate
  DATA_INTEGER(selex_type); // Selectivity type == 0 (length-based selectivity), == 1 (age-based selectivity)
  
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
  
  // Reference Point Parameters --------------
  // PARAMETER(ln_Fmsy); // Fmsy parameter in log space

  // PARAMETER TRANSFORMATIONS --------------------
  vector<Type> M = exp(ln_M); // Natural Mortality in normal Space
  Type sigmaRec = exp(ln_sigmaRec); // Recruiment sd in normal space
  Type sigmaRec2 = pow(sigmaRec, 2); // Recruitment Variability in normal Space
  // Type Fmsy = exp(ln_Fmsy); // Fmsy in normal space
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
  array<Type> Srv_AA(n_years, n_ages, n_sexes, n_srv_fleets); // Survey catch-at-age
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
  vector<Type> catch_sd(n_fish_fleets); // Empty container to store sd calculation
  vector<Type> fish_index_sd(n_fish_fleets); // Empty container to store sd calculation
  vector<Type> srv_index_sd(n_srv_fleets); // Empty container to store sd calculation

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
  if(selex_type == 0) {
    vector <Type> tmp_fishlens_slx(n_lens); // temporary container for fishery length selectivity
    vector <Type> tmp_fishage_slx(n_ages); // temporary container for fishery age selectivity
    for(int y = 0; y < n_years; y++) {
      for(int f = 0; f < n_fish_fleets; f++) {
        vector <Type> tmp_fish_selpars = ln_fish_selpars.transpose().col(f);
        for(int l = 0; l < n_lens; l++) tmp_fishlens_slx(l) = Logist(len_mids(l), tmp_fish_selpars);
        for(int s = 0; s < n_sexes; s++) {
          tmp_fishage_slx = age_len_transition.col(s).matrix() * tmp_fishlens_slx.matrix();
          for(int a = 0; a < n_ages; a++) Fish_Slx(y,a,s,f) = tmp_fishage_slx(a);
        } // end s loop
      } // end f loop
    } // end y loop
  } // end if for length-based selectivity
  
  for(int y = 0; y < n_years; y++) {
    for(int a = 0; a < n_ages; a++) {
      for(int s = 0; s < n_sexes; s++) {
        for(int f = 0; f < n_fish_fleets; f++) {
          
          // Compute Fishery Selectivity (age-based)
          if(selex_type == 1) {
            vector <Type> tmp_fish_selpars = ln_fish_selpars.transpose().col(s).col(f);
            Fish_Slx(y,a,s,f) = Logist(ages(a), tmp_fish_selpars);
          } // age-based selectivity
          
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
      // Get starting NAA
      NAA(0, a, s) = exp(RecPars(0)) * exp(-M(s) * Type(a)) * 
                          exp(ln_InitDevs(a)) * sexRatio(s); 
      
      if(a == 0) Total_Rec(0) += NAA(0, 0, s); // Get total recruitment - increment to get total recruitment
      Total_Biom(0) += NAA(0, a, s) * WAA(a,s); // Get total biomass

      // Calculate starting SSB
      if(s == 0 && n_sexes > 1) SSB(0) += NAA(0, a, 0) * MatAA(a) * WAA(a,0); // sex-specific assessment
      if(s == 0 && n_sexes == 1) SSB(0) += (NAA(0, a, 0) * MatAA(a) * WAA(a,0)) * 0.5; // sex-aggregated assessment
      
      // Compute Numbers-at-length in year 1
      vector<Type> NAL_tmp_vec = Convert_AL(age_len_transition, NAA, s, 0, 0, n_lens, 0); 
      for(int l = 0; l < n_lens; l++) NAL(0,l,s) = NAL_tmp_vec(l); // Loop through to input
    } // end age loop
  } // end sex loop
  
  // SSB0 Calculations ---------------------
  SBPR_N = Get_SBPR_N(n_ages, M(0)); // M(0) to index female natural mortality
  for(int a = 0; a < n_ages; a++) SBPR_SSB0(a) = SBPR_N(a) * MatAA(a) * WAA(a,0); // Compute SBPR0
  Type ssb0 = SBPR_SSB0.sum() * exp(RecPars(0)); // Get Virgin Equilibrium SSB here

  // Project population forward -------------
  for(int y = 1; y < n_years; y++){
    for(int s = 0; s < n_sexes; s++) {
      
      // Recruitment 
      Type ln_Det_BH_Rec = Get_Det_BH_Rec(RecPars, ssb0, SSB(y-1)); // deterministic beverton-holt recruitment
      NAA(y, 0, s) =  exp(ln_Det_BH_Rec + ln_RecDevs(y-1)) * sexRatio(s); // recruitment with process error
      Total_Rec(y) += NAA(y, 0, s); // Get total recruitment - increment to get total recruitment
      
      for(int a = 1; a < n_ages; a++) {
        // Project ages forward
        if(a < n_ages-1) {
          NAA(y,a,s) = NAA(y-1,a-1,s) * SAA(y-1,a-1,s);
        } // not plus group
        if(a == n_ages-1)  {
          NAA(y,n_ages-1,s) = (NAA(y-1,n_ages-1,s) * SAA(y-1,n_ages-1,s)) + // plus group decrement
                              (NAA(y-1,n_ages-2,s) * SAA(y-1,n_ages-2,s)); // increment into plus-group
        } // plus group
      } // end age loop
      
      // Compute Numbers-at-length 
      vector<Type> NAL_tmp_vec = Convert_AL(age_len_transition, NAA, s, y, 0, n_lens, 0); // sex-specific assessment
      for(int l = 0; l < n_lens; l++) NAL(y,l,s) = NAL_tmp_vec(l); // Loop through to input
      
      // Do residual calculations (total biomass and SSB)
      for(int a = 0; a < n_ages; a++) {
        Total_Biom(y) += NAA(y, a, s) * WAA(a,s);
        if(s == 0 && n_sexes > 1) SSB(y) += NAA(y,a,0) * MatAA(a) * WAA(a,0); // Calculate SSB here (Indexing 0 for females) (sex-specific assessment)
        if(s == 0 && n_sexes == 1) SSB(y) += (NAA(y, a, 0) * MatAA(a) * WAA(a,0)) * 0.5; // sex-aggregated assessment
      } // end a loop
      
    } // end sex loop
  } // end year loop
  
  // Calculate Catch Quantities -----------
  // vector<Type> CAL_tmp_vec(n_lens);
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
        vector<Type> CAL_tmp_vec = Convert_AL(age_len_transition, CAA, s, y, f, n_lens, 1); // sex-specific assessment
        for(int l = 0; l < n_lens; l++) CAL(y,l,s,f) = CAL_tmp_vec(l); // Loop through to input
      } // sex loop
    } // year loop
  } // end fishery fleet
  
  // Calculate Fishery Index ----------------------------
  for(int f = 0; f < n_fish_fleets; f++) {
    for(int s = 0; s < n_sexes; s++) {
      for(int y = 0; y < n_years; y++) { 
        for(int a = 0; a < n_ages; a++) {
          // Increment to get total index, conditional on fishery selectivity
          pred_fish_index(y,f) += exp(ln_q_fish(f)) * NAA(y,a,s) * 
                                  WAA(a,s) * Fish_Slx(y,a,s,f);
        } // a loop
      } // y loop
    } // s loop
  } // f loop 
  
  // Calculate Survey Index ----------------------------
  if(selex_type == 0) { // length-based survey selectivity
    vector <Type> tmp_srvlens_slx(n_lens); // temporary container for fishery length selectivity
    vector <Type> tmp_srvage_slx(n_ages); // temporary container for fishery age selectivity
    for(int y = 0; y < n_years; y++) {
      for(int sf = 0; sf < n_srv_fleets; sf++) {
        vector <Type> tmp_srv_selpars = ln_srv_selpars.transpose().col(sf);
        for(int l = 0; l < n_lens; l++) tmp_srvlens_slx(l) = Logist(len_mids(l), tmp_srv_selpars);
        for(int s = 0; s < n_sexes; s++) {
          tmp_srvage_slx = age_len_transition.col(s).matrix() * tmp_srvlens_slx.matrix();
          for(int a = 0; a < n_ages; a++) Srv_Slx(y,a,s,sf) = tmp_srvage_slx(a);
        } // end s loop
      } // end sf loop
    } // end y loop
  } // end if for length-based selectivity
  
  for(int sf = 0; sf < n_srv_fleets; sf++) {
    for(int s = 0; s < n_sexes; s++) {
      for(int y = 0; y < n_years; y++) {
        for(int a = 0; a < n_ages; a++) {
          
          // Calculate Survey Selectivity (age-based selectivity)
          if(selex_type == 1) {
            vector <Type> tmp_srv_selpars = ln_srv_selpars.transpose().col(s).col(sf);
            Srv_Slx(y,a,s,sf) = Logist(ages(a), tmp_srv_selpars); 
          } // if age-based survey selectivity
          
          // Increment to get total index, conditional on survey selectivity
          Srv_AA(y,a,s,sf) = NAA(y,a,s) * Srv_Slx(y,a,s,sf);
          pred_srv_index(y,sf) += exp(ln_q_srv(sf)) * Srv_AA(y,a,s,sf);
          
        } // a loop 
      } // y loop
    } // s loop
  } // sf loop  
  
  
  // Calculate Compositions and Related Quantities ------
  // Fishery Age Compositions
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
            for(int a = 0; a < n_ages; a++) pred_fish_age_comps(y,a,s,f) /= Total_Fish_Age_Numbers(y,s,f);
          } // if a = plus group and proportion within sexes
          
          // Normalize to sum to 1 (Proportions across sexes)
          if(a == n_ages - 1 && s == n_sexes - 1 && p_ow_sex_fish_age == 1) {
            // Get total fishery age numbers across sexes
            Type Total_fishAge_as_tmp = sum(Total_Fish_Age_Numbers.col(f).transpose().col(y));
            for(int s = 0; s < n_sexes; s++) {
              for(int a = 0; a < n_ages; a++) {
                pred_fish_age_comps(y,a,s,f) /= Total_fishAge_as_tmp;
              } // a loop
            } // s loop
          } // if a = plus group, s = n_sexes, and proportions across sexes
          
        } // a loop
      } // s loop
    } // f loop
  } // y loop
  
  // Fishery Length Compositions
  for(int y = 0; y < n_years; y++) {
    for(int f = 0; f < n_fish_fleets; f++) {
      for(int s = 0; s < n_sexes; s++) {
        
        // Pull out vector for lengths
        vector<Type> fish_lens_tmp = CAL.col(f).col(s).transpose().col(y).vec();
        
        // Computing Length Compositions
        for(int l = 0; l < n_lens; l++) {
          // Get predicted comps here for catch-at-length (prior to normalization)
          pred_fish_len_comps(y,l,s,f) = fish_lens_tmp(l);
          // Get total len numbers for normalizing
          Total_Fish_Len_Numbers(y,s,f) += pred_fish_len_comps(y,l,s,f);
          
          // Normalize to sum to 1 (Proportions within sexes)
          if(l == n_lens - 1 && p_ow_sex_fish_len == 0) {
            for(int l = 0; l < n_lens; l++) pred_fish_len_comps(y,l,s,f) /= Total_Fish_Len_Numbers(y,s,f);
          } // if l = n_lens - 1 and proportions within sexes
          
          // Normalize to sum to 1 (Proportions across sexes)
          if(l == n_lens - 1 && s == n_sexes - 1 && p_ow_sex_fish_len == 1) {
            // Extract total fishery lengths across exes
            Type Total_fishLen_as_tmp = sum(Total_Fish_Len_Numbers.col(f).transpose().col(y));
            for(int s = 0; s < n_sexes; s++) {
              for(int l = 0; l < n_lens; l ++) {
                pred_fish_len_comps(y,l,s,f) /= Total_fishLen_as_tmp;
              } // l loop
            } // s loop
          } // if l == n_lens - 1, s == n_sexes - 1, and proportions across sexes
          
        } // l loop
      } // sex loop
    } // fleet loop
  } // year loop
  
  // Survey Age Compositions 
  for(int y = 0; y < n_years; y++) {
    for(int sf = 0; sf < n_srv_fleets; sf++) {
      for(int s = 0; s < n_sexes; s++) {
        
        // Computing Age Compositions
        for(int a = 0; a < n_ages; a++) {
          // Get predicted comps here prior to normalizing w/ survey catch at age
          pred_srv_age_comps(y,a,s,sf) = Srv_AA(y,a,s,sf);
          // Increment to get total numbers at age for a given fleet
          Total_Srv_Age_Numbers(y,s,sf) += pred_srv_age_comps(y,a,s,sf);
          
          // Normalize to sum to 1 (Proportions within sexes)
          if(a == n_ages - 1 && p_ow_sex_srv_age == 0) {
            for(int a = 0; a < n_ages; a++) pred_srv_age_comps(y,a,s,sf) /= Total_Srv_Age_Numbers(y,s,sf);
          } // if a = plus group, and proportions within sexes
          
          // Normalize to sum to 1 (Proportions across sexes)
          if(a == n_ages - 1 && s == n_sexes - 1 && p_ow_sex_srv_age == 1) {
            // Extract survey ages to sum across sexes
            Type Total_srvAge_across_sex_tmp = sum(Total_Srv_Age_Numbers.col(sf).transpose().col(y));
            for(int s = 0; s < n_sexes; s++) {
              for(int a = 0; a < n_ages; a ++) {
                pred_srv_age_comps(y,a,s,sf) /= Total_srvAge_across_sex_tmp;
              } // a loop
            } // s loop
          } // if a == plus group, sex == n_sexes - 1, and proportions across sexes
          
        } // a loop
      } // s loop
    } // sf loop
  } // y loop
  
  // Survey length compositions
  for(int y = 0; y < n_years; y++) {
    for(int sf = 0; sf < n_srv_fleets; sf++) {
      for(int s = 0; s < n_sexes; s++) {
        // Computing Length Compositions (Getting Expected Lengths)
        vector<Type> srv_lens_tmp = Convert_AL(age_len_transition, Srv_AA, s, y, sf, n_lens, 1);
        
        for(int l = 0; l < n_lens; l++) {
          // Input temporary vector here
          pred_srv_len_comps(y,l,s,sf) = srv_lens_tmp(l);
          Total_Srv_Len_Numbers(y,s,sf) += pred_srv_len_comps(y,l,s,sf); // Sum up total survey lengths
          
          // Proportions within sexes (Normalize to 1)
          if(l == n_lens - 1 && p_ow_sex_srv_len == 0) {
            for(int l = 0; l < n_lens; l++) pred_srv_len_comps(y,l,s,sf) /= Total_Srv_Len_Numbers(y,s,sf);
          } // if l == n_lens - 1 and proportions within sexes
          
          // Normalize to sum to 1 (Proportions across sexes)
          if(l == n_lens - 1 && s == n_sexes - 1 && p_ow_sex_srv_len == 1) {
            // Extract total fishery lengths across exes
            Type Total_srvLen_as_tmp = sum(Total_Srv_Len_Numbers.col(sf).transpose().col(y));
            for(int s = 0; s < n_sexes; s++) {
              for(int l = 0; l < n_lens; l ++) {
                pred_srv_len_comps(y,l,s,sf) /= Total_srvLen_as_tmp;
                } // l loop
              } // s loop
            } // if l == n_lens - 1, s == n_sexes - 1, and proportions across sexes
          
          } // l loop
        } // sex loop
      } // survey fleet loop
    } // year loop
  
  // LIKELIHOOD SECTION --------------------
  // Catch Likelihood ----------------------
  for(int f = 0; f < n_fish_fleets; f++) {
    catch_sd(f) = cv_to_sd(catch_cv(f)); // Calculate sd
    for(int y = 0; y < n_years; y ++) {
      // Get likelihood here
      catch_nLL(y, f) -= use_catch(y,f) * dnorm(log(obs_catch(y, f)), 
                         log(pred_catch(y, f)) - pow(catch_sd(f), 2)/2, 
                         catch_sd(f), true);
    } // fish fleet loop
  } // year loop
  
  // Fishery Index Likelihoods -----------
  for(int f = 0; f < n_fish_fleets; f++) {
    fish_index_sd(f) = cv_to_sd(fish_index_cv(f)); // Calculate sd
    for(int y = 0 ; y < n_years; y++) {
      // Get likelihood here
      fish_index_nLL(y,f) -= use_fish_index(y,f) * dnorm(log(obs_fish_index(y, f)), 
                             log(pred_fish_index(y, f)) - pow(fish_index_sd(f),2)/2, 
                             fish_index_sd(f), true);
    } // y loop
  } // f loop
  
  // Survey Index Likelihoods -----------
  for(int sf = 0; sf < n_srv_fleets; sf++) {
    srv_index_sd(sf) = cv_to_sd(srv_index_cv(sf)); // Calculate sd
    for(int y = 0 ; y < n_years; y++) {
      // Get likelihood here
      srv_index_nLL(y,sf) -= use_srv_index(y,sf) * dnorm(log(obs_srv_index(y, sf)),
                             log(pred_srv_index(y, sf)) - pow(srv_index_sd(sf),2)/2,
                             srv_index_sd(sf), true);
    } // y loop
  } // sf loop
  
  // Fishery Composition Likelihoods -----
  // Fishery Age Compositions
  for(int f = 0; f < n_fish_fleets; f++) {
    for(int y = 0 ; y < n_years; y++) {

      // Age-Compositions Proportions across sex
      if(p_ow_sex_fish_age == 1) {
        vector<Type> obs_fish_age_as = obs_fish_age_comps.col(f).transpose().col(y).transpose().vec(); // Pull out observed vector
        vector<Type> pred_fish_age_as = pred_fish_age_comps.col(f).transpose().col(y).transpose().vec(); // Pull out observed vector
        fish_age_comp_nLL(y,0,f) -= use_fish_age_comps(y,f) * dmultinom(obs_fish_age_as, pred_fish_age_as, true); // Get likelihood here
      } // if proportions across sex for fishery ages
      
        // Proportions within sexes for fishery ages or lengths
        vector<Type> obs_fish_age_agg(n_ages); obs_fish_age_agg.setZero(); // storage containers for aggregated comps
        vector<Type> pred_fish_age_agg(n_ages); pred_fish_age_agg.setZero();// storage containers for aggregated comps
        
        if(p_ow_sex_fish_age == 0) {
          for(int s = 0; s < n_sexes; s++) {
            
            // Sex-specific comps
            if(agg_sex_fish_age == 0) {
              vector<Type> obs_fish_age_ws = obs_fish_age_comps.col(f).col(s).transpose().col(y); // Pull out observed vector
              vector<Type> pred_fish_age_ws = pred_fish_age_comps.col(f).col(s).transpose().col(y); // Pull out predicted vector
              fish_age_comp_nLL(y,s,f) -= use_fish_age_comps(y,f) * dmultinom(obs_fish_age_ws, pred_fish_age_ws, true); // Get likelihood here
            } // if sex-specific comps
            
            // Sex-aggregated comps
            if(agg_sex_fish_age == 1) {
              matrix<Type> obs_fish_age_mat = obs_fish_age_comps.col(f).transpose().col(y).transpose().matrix(); // dim = n_ages x n_sexes
              matrix<Type> pred_fish_age_mat = pred_fish_age_comps.col(f).transpose().col(y).transpose().matrix(); // dim = n_ages x n_sexes
              obs_fish_age_agg = obs_fish_age_mat.rowwise().sum(); // Sum across rows
              pred_fish_age_agg = pred_fish_age_mat.rowwise().sum(); // Sum across rows
              if(s == n_sexes - 1) {
                pred_fish_age_agg /= n_sexes; // divide by the number of sexes
                fish_age_comp_nLL(y,0,f) -= use_fish_age_comps(y,f) * dmultinom(obs_fish_age_agg, pred_fish_age_agg, true); // Get likelihood here
              } // if s == n_sexes-1, compute likelihood
            } // Aggregated comps
            
          } // s loop
        } // if proportions within sex for fishery ages
        
    } // y loop
  } // f loop
  
  // Fishery length compositions
  for(int f = 0; f < n_fish_fleets; f++) {
    for(int y = 0 ; y < n_years; y++) {

      // Length-Compositions Proportions across sex
      if(p_ow_sex_fish_len == 1) {
        vector<Type> obs_fish_len_as = obs_fish_len_comps.col(f).transpose().col(y).transpose().vec(); // Pull out observed vector
        vector<Type> pred_fish_len_as = pred_fish_len_comps.col(f).transpose().col(y).transpose().vec(); // Pull out observed vector
        fish_len_comp_nLL(y,0,f) -= use_fish_len_comps(y,f) * dmultinom(obs_fish_len_as, pred_fish_len_as, true); // Get likelihood here
      } // if proportions across sex for fishery lengths
      
      // Proportions within sexes for fishery ages or lengths
      if(p_ow_sex_fish_len == 0) {
        for(int s = 0; s < n_sexes; s++) {
          vector<Type> obs_fish_len_ws = obs_fish_len_comps.col(f).col(s).transpose().col(y); // Pull out observed vector
          vector<Type> pred_fish_len_ws = pred_fish_len_comps.col(f).col(s).transpose().col(y); // Pull out predicted vector
          fish_len_comp_nLL(y,s,f) -= use_fish_len_comps(y,f) * dmultinom(obs_fish_len_ws, pred_fish_len_ws, true);
          } // sex loop
        } // if proportions within sex for fishery lengths
        
      } // y loop
    } // f loop

  // Survey Composition Likelihoods ------
  // Survey Age Compositions
  for(int sf = 0; sf < n_srv_fleets; sf++) {
    for(int y = 0 ; y < n_years; y++) {
      
      // Age-Compositions Proportions across sexes
      if(p_ow_sex_srv_age == 1) {
        vector<Type> obs_srv_age_as = obs_srv_age_comps.col(sf).transpose().col(y).transpose().vec(); // Pull out observed vector
        vector<Type> pred_srv_age_as = pred_srv_age_comps.col(sf).transpose().col(y).transpose().vec(); // Pull out observed vector
        srv_age_comp_nLL(y,0,sf) -= use_srv_age_comps(y,sf) * dmultinom(obs_srv_age_as, pred_srv_age_as, true); // Get likelihood here
      } // if proportions across sex for fishery ages
      
      vector<Type> obs_srv_age_agg(n_ages); obs_srv_age_agg.setZero(); // storage containers for aggregated comps
      vector<Type> pred_srv_age_agg(n_ages); pred_srv_age_agg.setZero();// storage containers for aggregated comps
      
      if(p_ow_sex_srv_age == 0) {
        for(int s = 0; s < n_sexes; s++) {
          // Sex-specific comps
          if(agg_sex_srv_age == 0) {
            vector<Type> obs_srv_age_ws = obs_srv_age_comps.col(sf).col(s).transpose().col(y); // Pull out observed vector
            vector<Type> pred_srv_age_ws = pred_srv_age_comps.col(sf).col(s).transpose().col(y); // Pull out predicted vector
            srv_age_comp_nLL(y,s,sf) -= use_srv_age_comps(y,sf) * dmultinom(obs_srv_age_ws, pred_srv_age_ws, true); // Get likelihood here
            } // end if for sex-specific comps
          
          if(agg_sex_srv_age == 1) {
            // Extract out quantities
            matrix<Type> obs_srv_age_mat = obs_srv_age_comps.col(sf).transpose().col(y).transpose().matrix(); // dim = n_ages x n_sexes
            matrix<Type> pred_srv_age_mat = pred_srv_age_comps.col(sf).transpose().col(y).transpose().matrix(); // dim = n_ages x n_sexes
            obs_srv_age_agg = obs_srv_age_mat.rowwise().sum(); // Sum across rows
            pred_srv_age_agg = pred_srv_age_mat.rowwise().sum(); // Sum across rows
            if(s == n_sexes - 1) {
              pred_srv_age_agg /= n_sexes; // divide by the number of sexes
              srv_age_comp_nLL(y,0,sf) -= use_srv_age_comps(y,sf) * dmultinom(obs_srv_age_agg, pred_srv_age_agg, true); // Get likelihood here
              } // if s == n_sexes-1, compute likelihood
            } // Aggregated comps
            
          } // s loop
        } // if proportions within sex for fishery ages

    } // y loop
  } // sf loop
  
  // Survey length compositions
  for(int sf = 0; sf < n_srv_fleets; sf++) {
    for(int y = 0 ; y < n_years; y++) {
      
      // Length-Compositions Proportions across sexes
      if(p_ow_sex_srv_len == 1) {
        vector<Type> obs_srv_len_as = obs_srv_len_comps.col(sf).transpose().col(y).transpose().vec(); // Pull out observed vector
        vector<Type> pred_srv_len_as = pred_srv_len_comps.col(sf).transpose().col(y).transpose().vec(); // Pull out observed vector
        srv_len_comp_nLL(y,0,sf) -= use_srv_len_comps(y,sf) * dmultinom(obs_srv_len_as, pred_srv_len_as, true); // Get likelihood here
      } // if proportions across sex for fishery lengths

      if(p_ow_sex_srv_len == 0) {
        for(int s = 0; s < n_sexes; s++) {
          vector<Type> obs_srv_len_ws = obs_srv_len_comps.col(sf).col(s).transpose().col(y); // Pull out observed vector
          vector<Type> pred_srv_len_ws = pred_srv_len_comps.col(sf).col(s).transpose().col(y); // Pull out predicted vector
          srv_len_comp_nLL(y,s,sf) -= use_srv_len_comps(y,sf) * dmultinom(obs_srv_len_ws, pred_srv_len_ws, true);
        } // s loop
      } // if proportions within sex for fishery lengths
      
    } // y loop
  } // sf loop

  // Penalty for Population Initialization
  for(int a = 0; a < ln_InitDevs.size(); a ++) {
    rec_nLL -= dnorm(ln_InitDevs(a), -sigmaRec2/Type(2), sigmaRec, true);
  } 
  // Penalty for recruitment deviates
  for(int a = 0; a < ln_RecDevs.size(); a ++) {
    rec_nLL -= dnorm(ln_RecDevs(a), -sigmaRec2/Type(2), sigmaRec, true);
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
  REPORT(ssb0);
  REPORT(Total_Rec);
  REPORT(Total_Biom);
  
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




