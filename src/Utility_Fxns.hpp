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
  Type ln_BH_tmp_rec = log(BH_first_part/BH_sec_part); // Deterministic BH
  return ln_BH_tmp_rec;
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

