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
library(latex2exp)

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


plot_fwer_fdr = function(save_path = "paper/Figures/fwer_fdr_experiment.pdf", 
                         aggregate = FALSE,
                         probit_fits_dir = "data/probit_fits",
                         base_rate = 0.1, power = 0.8, 
                         significance_vals = c(0.2, 0.5, 0.8),
                         medium_sig_val = 0.5) {

  # Load data.
  dat = calc_fdr_vs_significance(significance_vals = significance_vals) %>%
    mutate(ExperimentID = paste(StudyID, ExperimentID, sep = "_")) %>%
    # select(!c(StudyID)) %>% 
    mutate(SigVal = as_factor(SigVal))

  # Filter to keep just "medium" effect sizes and sort factors.
  dat_med = 
    filter(dat, SigVal == medium_sig_val) %>% 
    mutate(ExperimentID = fct_reorder(ExperimentID, FDR))
    
  # Redefine default shape palette
  scale_shape_discrete = function(...) {
    scale_shape_manual(values = c(6, 1, 0))
  }

  p = 
    ggplot(dat_med, aes(x = FDR, y = ExperimentID)) + 
      geom_bar(aes(fill = StudyID, color = StudyID), stat = "identity", width = 0.15)  +
      labs(fill = "Study ID", color = "Study ID") +
      # geom_point(aes(x = FWER, y = ExperimentID, color = StudyID, fill = StudyID), 
      #            dat %>% filter(SigVal == medium_sig_val), size = 2, pch=21) + 
      geom_point(aes(x = FDR, y = ExperimentID, shape = SigVal, color = StudyID), 
                 dat, size = 2.25) + 
      labs(shape = TeX('Significance, $d^*$')) +
      xlim(c(0, 1)) + 
      xlab("False detection rate (FDR)") + 
      ylab("Experiment ID")

  ggsave(save_path, p)

  return (p)
}


plot_fwer_fdr()
