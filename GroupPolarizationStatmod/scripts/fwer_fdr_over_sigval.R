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

calc_tibble = function(probit_fits_dir = "data/probit_fits",
                       base_rate = 0.1, power = 0.8, agglevel = "experiment", 
                       significance_vals = c(0.2, 0.5, 0.8)
                       n_experiments = 54) {

  # Initialize full experiments tibble.
  full_tbl = tibble(StudyID = character(n_experiments),
               ExperimentID = character(n_experiments),
               FWER = numeric(n_experiments),
               SigVal = numeric(n_experiments),
               Power = numeric(n_experiments),
               BaseRate = numeric(n_experiments))
  
  # Calculate FWER for each significance value.
  ii = 1
  while (ii <= 3) {
    
    # This reloads the data every time, but it's small data, so forget it.
    low_rng = ((ii-1) * n_experiments + 1)
    high_rng = (ii * n_experiments)
    update_rng = low_rng:high_rng

    tbl[update_rng, ] = 
      col_bind(
        calculate_fwer(probit_fits_dir, significance_vals[ii]),
        rep(significance_vals[ii], n_experiments),
        rep(power, n_experiments),
        rep(base_rate, n_experiments)
      )

    ii = ii + 1
  }

  return (tbl)
}



