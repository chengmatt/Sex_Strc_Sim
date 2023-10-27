# Sex Structure Simulations

## Purpose
To evaluate good practices for the treatment of sex-structure in integrated stock assessments, with several questions in mind: 
1)	Treatment of composition proportions (i.e., proportions sum to 1 within each sex, or sum to 1 across sexes), 
2)	Utility of sex-specific catch data, 
3)	Consequences of estimating an aggregated M, when sex-specific M is present,
4)	Estimability of initial sex ratios at recruitment,
5)	Consequences of mis-specfiying initial sex ratios


### Repository Structure
| Folder  | Items |
| --------| --------|
|R| Contains R scripts for running simulations, utility functions, etc. |
|figs| Contains general figures and model outputs |
|input| Excel files for simulation inputs |
|output| Output files for simulations |
|src| Contains source code for TMB model, which is compiled through R|

### Operating Model Options
The OM is a single area, age-structure population model with  several fleets (fishery and survey) coded in, but is set up for only single fleet observations (can be modified as needed). Additionally, the OM can generate both age and length composition data, and is able to accomodate multiple sexes (gonochoristic).
| OM Component  | Options |
| --------| --------|
|Recruitment| Beverton Holt Recruitment |
|Fishing mortality pattern| Contrast |
|Selectivity| Logistic |
|Compositional Likelihoods| Multinomial |

### Estimation Model Options
The EM is able to be specified for multiple fishery and survey fleets, as well as multiple sexes. 
| EM Component  | Options |
| --------| --------|
|Recruitment| Beverton Holt Recruitment |
|Selectivity| Logistic |
|Natural Mortality| Sex-specific M estiamted as an offset, or single aggregated M |
|Compositional Proportions| Proportions sum to 1 across sexes, without a sex-ratio nLL, or Proportions sum to 1 within sex, while fitting to a sex-ratio nLL (sex-ratios can be fit to as within a given age or year, or within a given year)|
|Initial Sex Ratios| Can be fixed or estimated in logit space, bounded between 0 and 1 |
|Compositional Likelihoods| Multinomial |

