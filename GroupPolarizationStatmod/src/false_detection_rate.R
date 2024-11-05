##
# Calculate the false detection rate from bayesian fits data.
#
# Author: Matthew A. Turner
# Date: 2024-10-30
#

library(tidyverse)
library(fs)


cohens_d <- function(pre_mean_estimate, post_mean_estimate, pre_sd, post_sd) {
  
  numerator <- post_mean_estimate - pre_mean_estimate
  
  denominator = sqrt((pre_sd**2 + post_sd**2) / 2.0)
  
  return (numerator / denominator)
}


load_fdr_data = function(probit_fits_dir = "data/probit_fits", 
                         sync_file = "data/probit_fits/all.csv",
                         overwrite = FALSE) {
   
  if (!file.exists(sync_file) || overwrite) {

    ret <- dir_ls(probit_fits_dir, glob = "*.csv") %>%
      read_csv() %>% unite(ExperimentID, ArticleTag, TreatmentTag)

    write_csv(ret, sync_file)

  } else {

    ret <- read_csv(sync_file)
  }

  return (ret)
}


false_detection_rate_by_condition <- 
    function(probit_fits_dir = "data/probit_fits", sigval = 0.8) {
  
  fdr_data = load_fdr_data(probit_fits_dir) 

  unique_treatments = unique(fdr_data$TreatmentTag)

  n_plausible = length(unique_treatments)

  Treatment = unlist(
    map(unique_treatments, \(treatment) rep(treatment, sigvals))
  )

  alpha_v_sigval = tibble(Alpha = rep(NA, n_plausible*length(sigvals)),
                          Sigval = rep(sigvals, n_plausible),
                          Treatment = Treatment)

  for (this_sigval in sigvals) {
    gb = load_fdr_data(fdr_data_file) %>%
           select(Sigval == this_sigval) %>%
           group_by(TreatmentTag)

    # Calculate alpha from this gb and put in alpha_v_sigval.
    # Or somehow need to save one for each treatment or sig val then
    # rbind all the results.
    
  }
}


##
# Family-wise error rate calculation using whatever ExperimentIDs are
#
#
calculate_fwer = function (probit_fits_dir = "data/probit_fits", sigval = 0.8) {

  fdr_data = load_fdr_data(probit_fits_dir)
     
  fdr_data$Cohens_d = cohens_d(fdr_data$LatentMeanPrePosteriorMean, 
                                fdr_data$LatentMeanPostPosteriorMean,
                                fdr_data$LatentMeanPrePosteriorSD, 
                                fdr_data$LatentMeanPostPosteriorSD)

  fdr_data$Significant = ifelse(fdr_data$Cohens_d >= sigval, 1, 0)

  # return (fdr_data)
  ret = fdr_data %>%
    group_by(ExperimentID) %>% 
    summarise(FWER = mean(Significant)) %>%
    separate_wider_delim(ExperimentID, delim = "_", 
                         names = c("StudyID", "ExperimentID"))

  return (ret)
}


false_detection_rate = function(by = "condition") {
  if (by == "condition") {
    
  }
}


ret = calculate_fwer()

