// Purpose: An assortment of utility functions
// Creator: Matthew LH. Cheng (UAF-CFOS)
// Date: 8/9/23

// Function to compute logistic function
template <class Type> 
Type Logist(Type age, // age integer
            vector<Type> ln_pars // vector of parameters n = 2
              ) {
  // Transform parameters
  Type slope = exp(ln_pars(0));
  Type midpoint = exp(ln_pars(1));
  Type Selex = Type(1) / (Type(1) + exp(-slope * (age - midpoint)));
  return Selex;
} // end fxn

// Function to compute logistic function using the a50 and a95 parameterization
template <class Type> 
Type Logist_19(Type age, // age integer
            vector<Type> ln_pars // vector of parameters n = 2
) {
  // Transform parameters
  Type a50 = exp(ln_pars(0));
  Type a95 = exp(ln_pars(1));
  Type Selex = Type(1.0) / (Type(1.0) + pow(Type(19.0), (a50 - age) / a95));
  return Selex;
} // end fxn

// Get SBPR in numbers
template <class Type> 
vector<Type> Get_SBPR_N(int n_ages, // integer of max age
                        Type M // natural mortality
                        ) {
  // SBPR container
  vector<Type> SBPR_N(n_ages); 
  // Loop through to get SBPR in numbers
  for(int a = 0; a < n_ages; a++) {
    if(a == 0) SBPR_N(a) = Type(1);
    if(a > 0 && a < n_ages - 1) SBPR_N(a) = SBPR_N(a - 1) * exp(-M); 
    if(a == n_ages - 1) SBPR_N(a) = ((SBPR_N(a - 1) * exp(-M)) / (1 - exp(-M)));
  } // age loop
  return SBPR_N;
} // end function

// Get Beverton-Holt recruitment
template <class Type>
Type Get_Det_BH_Rec(vector<Type> RecPars, // vector of recruitment parameters (ln_R0, and h)
                    Type ssb0,// spawning stock biomass viring
                    Type SSB // ssb from previous year
                    ) {
  // Recruitment
  Type R0 = exp(RecPars(0)); // Virgin Recruitment in normal space
  Type h = RecPars(1); // Steepness in normal space
  // Define BH components
  Type BH_first_part = Type(4) * h * R0 * SSB;
  Type BH_sec_part = (ssb0 * (Type(1) - h)) + SSB * ((Type(5)*h) - 1);
  Type BH_tmp_rec = BH_first_part/BH_sec_part; // Deterministic BH
  return BH_tmp_rec;
} // end function

// Go from CV to SD
template <class Type>
Type cv_to_sd(Type cv) {
  Type sd = 0;
  sd = sqrt(log( ((cv * cv) + 1)));
  return sd;
} // end function

// Go from age to length using age-transition matrix for NAA/CAA/SurveyAA
template<class Type>
vector<Type> Convert_AL(array<Type> age_len_trans_mat, // age length transition matrix, dim = age, length, sex
                              array<Type> numbers, // numbers at age, CAA or SurveyAA, dim = years, age, sex, fleet
                              int sex, // sex
                              int year, // year
                              int fleet, // optional input if we are doing CAA/survey
                              int n_lens, // nmber of length bins
                              int type // indicator for numbers at age or CAA/survey (0 = NAA, 1 = CAA/Survey)
                              ) {
  vector<Type> lens(n_lens); // initialize vector
  if(type == 0) { // for just numbers at age
    // matrix multiplication to convert age to len
    lens = age_len_trans_mat.col(sex).transpose().matrix() * 
      numbers.col(sex).transpose().col(year).matrix(); 
  } // if this is numbers at age
  if(type == 1) { // for CAA or survey
    lens = age_len_trans_mat.col(sex).transpose().matrix() *
      numbers.col(fleet).col(sex).transpose().col(year).matrix();
  } // if this is CAA or for a survey
  return lens;
} // end function

// Get Spawning Biomass Per Recruit
template<class Type>
Type Get_SBPR(Type F, // trial F value
             vector<Type> selex, // vector of age-specific selectivities
             Type M, // Female natural mortality
             vector<Type> waa, // vector of weight at ages
             vector<Type> MatAA, // vector of maturity at ages
             vector<Type> ages // vector of ages
               ) {
  
  Type Na = 1; // Define initial population
  vector<Type> Za = (selex * F) + M; // Get total mortality
  vector<Type> Srva = exp(-Za); // Get survival
  Type SBPR = Na * exp(-Za(0)) * waa(0) * MatAA(0); // initialize SPR with first age class
  
  // Loop through to get the rest of the quantities
  for(int a = 1; a < (ages.size() * 4); a++) {
    if(a >= ages.size() ) { // if a is at the max iteration
      SBPR += Na * exp(-Za(ages.size() - 1)) * waa(ages.size() - 1) * MatAA(ages.size() - 1);
      Na *= Srva(ages.size() - 1);
    } else {
      SBPR += Na * exp(-Za(a)) * waa(a) * MatAA(a);
      Na *= Srva(a);
    } // end if else for a >= max iteration
  } // end a loop
  
  return SBPR;
} // end get spr function


// Get Yield Per Recruit
template<class Type>
Type Get_YPR(Type F, // trial F value
             vector<Type> selex, // vector of age-specific selectivities
             Type M, // Female natural mortality
             vector<Type> waa, // vector of weight at ages
             vector<Type> ages // vector of ages
               ) {
  
  Type Na = 1; // initialize population
  Type YPR = ((selex(0) * F) / (M + selex(0) * F)) * // Get YPR for the first age class
             Na * (Type(1) - exp(-(M + selex(0) * F))) * waa(0);
  // Finish filling in quantities w/ age loop
  for(int a = 1; a < ages.size() * 4; a++) {
    if(a >= ages.size()) { // if max iteration
      Na *= exp(-(M + selex(ages.size() - 1) * F));
      YPR += (selex(ages.size() - 1) * F) / (M + selex(ages.size() - 1) * F) * 
             Na * (Type(1) - exp(-(M + selex(ages.size() - 1) * F))) * waa(ages.size() - 1);
    } else{
      Na *= exp(-(M + selex(a) * F));
      YPR += (selex(a) * F) / (M + selex(a) * F) * 
             Na * (Type(1) - exp(-(M + selex(a) * F))) * waa(a);
    } // end if else for max iteration loop
  } // end a loop
  
  return YPR;
} // end Get_YPR function
  
  
// Get equilibrium SSB for Bmsy
template<class Type>
Type Get_SSBe(Type F, // trial F value
              vector<Type> N_init, // vector of inital numbers at age
              vector<Type> selex, // vector of age-specific selectivities
              Type M, // Female natural mortality
              vector<Type> waa, // vector of weight at ages
              vector<Type> MatAA, // vector of maturity at ages
              vector<Type> ages // vector of ages
                ) {
  // Set up
  Type ssbe = 0; // equilibrium ssb
  vector<Type> Fa = F * selex; // fishing mortality at age
  vector<Type> Za = Fa + M; // total mortality at age
  vector<Type> Sa = exp(-Za); // survival at age
  array<Type> N_equilibrium(ages.size(), ages.size() * 2); // set up array 
  N_equilibrium.col(0) = N_init; // input initial numbers at age into the first year
  
  // Run annual cycle 
  for(int y = 1; y < (ages.size() * 2); y++) {
    // Initial recruitment
    N_equilibrium(0, y) = N_init(0);
    for(int a = 1; a < ages.size(); a++) 
      N_equilibrium(a, y) = N_equilibrium(a - 1, y - 1) * Sa(a - 1); // project population forward (not plus group)
      N_equilibrium(ages.size() - 1, y) = N_equilibrium(ages.size()  - 2, y - 1) * Sa(ages.size()  - 2) +
                                          N_equilibrium(ages.size()  - 1, y - 1) * Sa(ages.size()  - 1); // plus group calculations
  }
  
  // Calculate SSB equilibrium
  for(int a = 0; a < ages.size(); a++) {
    ssbe += N_equilibrium(a, (ages.size() * 2) - 1) * exp(-Za(a)) * MatAA(a) * waa(a);
  } // end a loop

  return(ssbe);
} // end Get SSBe function
