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
                         significance_vals = c(0.2, 0.5, 0.8)) {

  # Load data.
  p = calc_fdr_vs_significance(significance_vals = significance_vals) %>%
    # Concatenate ExperimentID name.
    mutate(ExperimentID = paste(StudyID, ExperimentID, sep = "_"), 
           SigVal = as_factor(SigVal)) %>%
    arrange(desc(FWER)) %>%
    ggplot(aes(x = FWER, y = ExperimentID, shape = SigVal)) + geom_point()
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
