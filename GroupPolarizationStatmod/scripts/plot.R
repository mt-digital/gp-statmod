##
#
# Plot model illustrations and analysis of family-wise error rate simulations.
# 
# Currently model illustrations are in the root project directory and will
# soon be ported here.
# 
# Author: Matthew A. Turner <maturner01@gmail.com>
# Date: 2024-11-05
#

library(forcats)
library(ggplot2)

source("scripts/analysis.R")


plot_sigval_for_low_fwer = function(target_fwer = 0.05, 
                                    probit_fits_dir = "data/probit_fits",
                                    base_rate = 0.1, power = 0.8, 
                                    significance_vals = seq(0.2, 3.0, 0.01)) {

  
  sigval_for_low_fwer(calc_fdr_vs_significance(significance_vals = 
                                               significance_vals)) %>%
    mutate(ExperimentID = paste(StudyID, ExperimentID, "_")) %>%
    ggplot(aes(x = SigVal, y = ExperimentID, color = ))
    

}


plot_fwer_fdr = function(probit_fits_dir = "data/probit_fits",
                         base_rate = 0.1, power = 0.8, 
                         significance_vals = c(0.2, 0.5, 0.8),
                         medium_sig_val = 0.5) {

  # Load data.
  dat = calc_fdr_vs_significance(significance_vals = significance_vals) %>%
    mutate(ExperimentID = paste(StudyID, ExperimentID, sep = "_")) %>%
    select(!c(StudyID)) %>% 
    mutate(SigVal = as_factor(SigVal))

  # Filter to keep just "medium" effect sizes and sort factors.
  dat_med = 
    filter(dat, SigVal == medium_sig_val) %>% 
    mutate(ExperimentID = fct_reorder(ExperimentID, FDR))
    
  p = 
    ggplot(dat_med, aes(x = FDR, y = ExperimentID)) + 
      geom_bar(stat = "identity", fill = "red", color = "red", width = 0.15)  +
      geom_point(aes(x = FWER, y = ExperimentID), 
                 dat %>% filter(SigVal == medium_sig_val), size = 2) + 
      geom_point(aes(x = FDR, y = ExperimentID, shape = SigVal), dat) +
      xlim(c(0, 1))

      


  # %>%
  #   # Concatenate ExperimentID name.
  #   mutate(ExperimentID = paste(StudyID, ExperimentID, sep = "_"), 
  #          SigVal = as_factor(SigVal)) %>%
  #   # arrange(desc(FWER)) %>%
  #   ggplot(aes(x = FWER, y = ExperimentID, shape = SigVal)) + geom_point()
    # Tall-ify table to have an FDR and FWER for each significance value.
    # melt(id.vars = c("StudyID", "ExperimentID", "SigVal", "Power", "BaseRate"),
    #      value.name = "Value", variable.name = "Measure") %>%
    # Plot...
    # ggplot(aes(x = Value, y = ExperimentID, color = SigVal, shape = Measure)) +
    #   geom_point()

  ggsave("test.pdf", p)

  return (p)
}


plot_fwer_fdr()
