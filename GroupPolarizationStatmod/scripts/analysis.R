##
# Two functions. First, calculate a tibble that calculates the family-wise
# error rate (FWER) and false detection rate (FDR) for three significance values
# d^*=0.2,0.5,0.8 at three different levels of aggregation (experimental condition, 
# journal article, and all). The significance values correspond to Cohen's
# (1988; pp. ~25-30) recommendation. 
#
# Second, visualize the tibble at the different aggregation levels for building
# into a composite figure.
# 
# Author: Matthew A. Turner <maturner01@gmail.com>
# Date: 2024-11-05
#


source("src/false_detection_rate.R")


fdr <- function(fwer, br, p) {

  num = (1 - br) * fwer
  denom = (num + (p * br))

  return (num / denom)
}


calc_fdr_vs_significance_tibble = function(
  probit_fits_dir = "data/probit_fits", base_rate = 0.1, power = 0.8, 
  agglevel = "experiment", significance_vals = c(0.2, 0.5, 0.8),
  n_experiments = 54
) {

  # Initialize return tibble.
  ret_tbl = tibble(StudyID = character(n_experiments),
               ExperimentID = character(n_experiments),
               FWER = numeric(n_experiments),
               FDR = numeric(n_experiments),
               SigVal = numeric(n_experiments),
               Power = numeric(n_experiments),
               BaseRate = numeric(n_experiments))

  probit_fits_data = load_probit_data()
  
  # Calculate FWER for each significance value.
  ii = 1
  while (ii <= length(significance_vals)) {
    
    # Calculate row indices to be updated in return tibble.
    low_rng = ((ii-1) * n_experiments + 1)
    high_rng = (ii * n_experiments)
    update_rng = low_rng:high_rng

    # Extract this significance value for readability.
    significance_val = significance_vals[ii]

    # Get tibble with StudyID, ExperimentID, and FWER columns...
    fwer_tbl = calculate_fwer(probit_fits_data, significance_val)
    # ...to calculate FDR and create a combined table.
    fdr_col = tibble(FDR = fdr(fwer_tbl$FWER, base_rate, power))
    fwer_fdr = bind_cols(fwer_tbl, fdr, col)
                
    # Update the return tibble with this significance_val.
    ret_tbl[update_rng, ] = 
      bind_cols(
        fwer_fdr,
        # Three metadata columns.
        tibble(SigVal = rep(significance_vals[ii], n_experiments),
               Power = rep(power, n_experiments),
               BaseRate = rep(base_rate, n_experiments))
      )
    
    ii = ii + 1
  }

  return (ret_tbl)
}


sigval_for_low_fwer = function(fdr_vs_sig_tbl, target_fwer = 0.05) {

  return(

    fdr_vs_sig_tbl %>%
      filter(tbl, FWER <= 0.05) %>% 
      group_by(StudyID, ExperimentID) %>% 
      filter(SigVal == min(SigVal))  %>% 
      arrange(desc(SigVal)) 

  )
}


add_default_fwer_rows = function(fdr_vs_sig_tbl, default_fwer = 0.05) {

  # Figure out the number of sig values used.
  n_sigvals = length(unique(fdr_vs_sig_tbl$SigVal))

  # Extract power and base rate.
  power = fdr_vs_sig_tbl[1, Power]
  
  # Use default FWER, FDR for all SigVals for non-plausibly-identified studies.
   
}


aggregate_fdr_vs_sig = function(fdr_vs_sig_tbl, default_fwer = 0.05) {
  
  ret = 
    add_default_fwer_rows(fdr_vs_sig_tbl, default_fwer) %>%
      group_by(StudyID, SigVal) %>%
      summarise(FWER = mean(FWER), FDR = mean(FDR))
}


NON_IDENTIFIED_STUDIES_EXPERIMENTS = 
  tibble(StudyID = c("Schkade2010", "Schkade2010", "Myers1975"),
         ExperimentID = c("COSprings-GlobalWarming", "Boulder-Aff.Act.", "Bad-Experimental"))

