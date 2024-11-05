##
# Calculate the false detection rate from bayesian fits data.
#
# Author: Matthew A. Turner
# Date: 2024-10-30
#

# library(tidyverse)
library(dplyr)
library(fs)


cohens_d <- function(pre_mean_estimate, post_mean_estimate, pre_sd, post_sd) {
  
  numerator <- post_mean_estimate - pre_mean_estimate
  
  denominator = sqrt((pre_sd**2 + post_sd**2) / 2.0)
  
  return (numerator / denominator)
}


load_probit_data = function(probit_fits_dir = "data/probit_fits", 
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


##
# Family-wise error rate calculation using whatever ExperimentIDs are present.
# Note probit_fits may be a tibble for avoiding reading data every calculation.
#
calculate_fwer = function (probit_fits = "data/probit_fits", sigval = 0.8) {

  # Load probit data if necessary.
  if (typeof(probit_fits) == "character") {
    probit_data = load_probit_data(probit_fits)
  } else {
    probit_data = probit_fits
  }
     
  probit_data$Cohens_d = cohens_d(probit_data$LatentMeanPrePosteriorMean, 
                                  probit_data$LatentMeanPostPosteriorMean,
                                  probit_data$LatentMeanPrePosteriorSD, 
                                  probit_data$LatentMeanPostPosteriorSD)

  probit_data$Significant = ifelse(probit_data$Cohens_d >= sigval, 1, 0)

  # return (probit_data)
  ret = probit_data %>%
    group_by(ExperimentID) %>% 
    summarise(FWER = mean(Significant)) %>%
    separate_wider_delim(ExperimentID, delim = "_", 
                         names = c("StudyID", "ExperimentID")) %>%
    arrange(desc(FWER))

  return (ret)
}


