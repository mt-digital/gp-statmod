##
# Utility functions used in other src and scripts files.
#
# Author: Matthew A. Turner
# Date: 2024-12-27
#

library(dplyr)
library(readr)
library(tidyr)
library(fs)

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


