##
# 
#

load_fdr_data = function(probit_fits_dir = "data/probit_fits", 
                         sync_file = "data/probit_fits/all.csv",
                         overwrite = FALSE) {
   
  if (file_exists(sync_file) && !overwrite) {

    ret <- dir_ls(ordinal_data_dir, glob = "*.csv") %>%
      read_csv() %>% unite(ExperimentID, ArticleTag, TreatmentTag)

    write_csv(df, sync_file)

  } else {

    ret <- read_csv(sync_file)
  }

  return (ret)
}


false_detection_rate_by_condition = function(probit_fits_dir = "data/probit_fits",
                                             sigvals = seq(0.2, 2.0, 0.15)) {
  
  fdr_data = read_csv(load_fdr_data(probit_fits_dir)) 

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
    # rbind all the results
    
  }
         

}

