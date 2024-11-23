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


# b for base rate, W for power. See 
fdr = function(fwer, b, W) {

  num = (1 - b) * fwer
  denom = (num + (b * W))

  return (num / denom)
}


calc_fdr_vs_significance = function(
  probit_fits_dir = "data/probit_fits", base_rate = 0.1, power = 0.8, 
  significance_vals = c(0.2, 0.5, 0.8)
) {

  probit_fits_data = load_probit_data()
  n_experiments = length(unique(probit_fits_data$ExperimentID))

  # Initialize return tibble.
  ret_tbl = tibble(StudyID = character(n_experiments),
                   ExperimentID = character(n_experiments),
                   FWER = numeric(n_experiments),
                   FDR = numeric(n_experiments),
                   SigVal = numeric(n_experiments),
                   Power = numeric(n_experiments),
                   BaseRate = numeric(n_experiments))
  
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
    fwer_fdr = bind_cols(fwer_tbl, fdr_col)
                
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


add_default_rows = function(fdr_vs_sig_tbl, default_fwer = 0.05) {

  # Figure out the number of sig values used.
  sigvals = fdr_vs_sig_tbl$SigVal
  n_sigvals = length(unique(sigvals))

  # Extract power and base rate.
  W = fdr_vs_sig_tbl[1, Power]
  b = fdr_vs_sig_tbl[1, BaseRate]
  
  # Load and use non-identified plausible.
  non_id = NON_IDENTIFIED

  n_non_id = n_distinct(non_id)
  n_default_rows = n_non_id * n_sigvals

  # Set default_fwer to be half of the least fwer calculated for that article
  # by first finding the minimum fwer for the study.
  study_ids = non_id$StudyID

  default_fwers_tbl = 
    fdr_vs_sig_tbl %>% 
    filter(StudyID %in% study_ids) %>%
    group_by(StudyID) %>%
    summarize(default_fwer = min(FWER))

  # Then create a vector that fits with default_rows.
  default_fwers_vec = map(default_frews$DEFAULT, 
                          \(frew) rep(frew / 2.0, n_sigvals))

  default_fdr_vec = fdr(default_fwers_vec, b, W)

  default_rows = tibble(StudyID = non_id$StudyID,
                        ExperimentID = non_id$ExperimentID,
                        FWER = default_fwers_vec,
                        FDR = default_fdr_vec,
                        SigVal = rep(sigvals, n_non_id),
                        Power = rep(W, n_default_rows),
                        BaseRate = rep(b, n_default_rows))
  
  # Use default FWER, FDR for all SigVals for non-plausibly-identified studies.
  return (bind_rows(fdr_vs_sig_tbl, default_rows))
}


aggregate_fdr_vs_sig = function(fdr_vs_sig_tbl, default_fwer = 0.05) {
  
  ret = 
    add_default_rows(fdr_vs_sig_tbl, default_fwer) %>%
      group_by(StudyID, SigVal) %>%
      summarise(FWER = mean(FWER), FDR = mean(FDR))
}


NON_IDENTIFIED = 
  tibble(
    StudyID = c("Schkade2010", "Schkade2010", "Myers1975"),
    ExperimentID = c("COSprings-GlobalWarming", "Boulder-Aff.Act.", "Bad-Experimental"),
    DefaultFWER = c(0.5, 0.5, 0.2)
  )


##
# Given fdr_vs_sig_tbl (which includes an FWER column), 
# calculate the smallest significance value that achieves an FWER of 
# target_fwer, set to a default main text value of 0.05.
#
sigval_for_low_fwer = function(fdr_vs_sig_tbl, target_fwer = 0.05) {

  return(

    fdr_vs_sig_tbl %>%
      filter(FWER <= target_fwer) %>% 
      group_by(StudyID, ExperimentID) %>% 
      filter(SigVal == min(SigVal))  %>% 
      mutate(SigVal = as_factor(SigVal)) %>%
      arrange(desc(SigVal)) 

  )
}
