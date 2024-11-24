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
                                    save_path = "paper/Figures/sigval_for_low_fwer.pdf",
                                    probit_fits_dir = "data/probit_fits",
                                    base_rate = 0.1, power = 0.8, 
                                    significance_vals = seq(0.2, 3.0, 0.01)) {
  
  # Going to create ggplot, p.
  p = 
    # First load and pre-process calculation of FWER ≤ target_fwer by experiment.
    sigval_for_low_fwer(calc_fdr_vs_significance(significance_vals = 
                                                 significance_vals)) %>%
      mutate(ExperimentID = paste(StudyID, ExperimentID, "_"),
             SigVal = as.numeric(SigVal)) %>%

    # Then plot the result.
    ggplot(aes(x = SigVal, y = ExperimentID, fill = StudyID, color = StudyID)) +
      geom_bar(stat = "identity", width = 0.15)  +
      xlab(TeX("Significance $d^*$ for $\\alpha \\leq 0.05")) + 
      ylab("Experiment ID")
        
     ggsave(save_path, p) 

   return (p)
}


plot_fdr = function(save_path = "paper/Figures/fdr_vs_experiment.pdf", 
                    aggregate = FALSE, use_non_identified = FALSE,
                    probit_fits_dir = "data/probit_fits",
                    default_fwer = 0.05,
                    base_rate = 0.1, power = 0.8, 
                    significance_vals = c(0.2, 0.5, 0.8),
                    medium_sig_val = 0.5) {

  # Load data.
  tbl = calc_fdr_vs_significance(significance_vals = significance_vals) %>%
    mutate(ExperimentID = paste(StudyID, ExperimentID, sep = "_")) %>%
    # select(!c(StudyID)) %>% 
    mutate(SigVal = as_factor(SigVal), Power = power, BaseRate = base_rate)

  if (aggregate)
    p = plot_fdr_aggregated(tbl, medium_sig_val, save_path, default_fwer, use_non_identified)
  else
    p = plot_fdr_nonaggregated(tbl, medium_sig_val, save_path)

  return (p)
}

# Redefine default shape palette
scale_shape_discrete = function(...) {
  scale_shape_manual(values = c(6, 1, 0))
}

plot_fdr_aggregated = function(tbl, medium_sig_val, 
                               save_path, default_fwer = 0.05,
                               use_non_identified = FALSE) {

  print(use_non_identified)
  
  tbl = as_tibble(aggregate_fdr_vs_sig(tbl, use_non_identified, default_fwer))

  tbl_med = 
    filter(tbl, SigVal == medium_sig_val) %>%
    select(!SigVal) %>%
    mutate(StudyID = fct_reorder(StudyID, FDR))

  cat(
    paste0("\n\n***** CALCULATED MEAN FDR = ", 
           mean(tbl_med$FDR), 
           ifelse(use_non_identified, " with", " without"),
           " default FWER", 
           ifelse(use_non_identified, paste0(" = ", default_fwer), ""),
           " *****\n\n"
    )
  )

  print(tbl_med)

  p = 
    ggplot(tbl_med, aes(x = FDR, y = StudyID)) +
      geom_bar(aes(fill = StudyID, color = StudyID), stat = "identity", width = 0.15) +
      labs(fill = "Study ID", color = "Study ID") +
      geom_point(aes(x = FDR, y = StudyID, shape = SigVal, color = StudyID),
                 tbl, size = 5) + 
      labs(shape = TeX("Significance, $d^*$")) +
      xlim(c(0, 1)) +
      xlab("False detection rate (FDR)") +
      ylab("Study ID") +
      guides(color = guide_legend(reverse = TRUE), fill = guide_legend(reverse = TRUE))

  ggsave(save_path, p)

  return (p)
}

plot_fdr_nonaggregated = function(tbl, medium_sig_val, save_path) {
  
  # Filter to keep just "medium" effect sizes and sort factors.
  tbl_med = 
    filter(tbl, SigVal == medium_sig_val) %>% 
    mutate(ExperimentID = fct_reorder(ExperimentID, FDR))

  p = 
    ggplot(tbl_med, aes(x = FDR, y = ExperimentID)) + 
      geom_bar(aes(fill = StudyID, color = StudyID), stat = "identity", width = 0.15)  +
      labs(fill = "Study ID", color = "Study ID") +
      geom_point(aes(x = FDR, y = ExperimentID, shape = SigVal, color = StudyID), 
                 tbl, size = 2.25) + 
      labs(shape = TeX('Significance, $d^*$')) +
      xlim(c(0, 1)) + 
      xlab("False detection rate (FDR)") + 
      ylab("Experiment ID") +
      guides(color = guide_legend(reverse = TRUE), fill = guide_legend(reverse = TRUE))

  ggsave(save_path, p)

  return (p)
}

plot_fdr()

plot_fdr(aggregate = TRUE, save_path = "paper/Figures/fdr_vs_study.pdf")

plot_fdr(aggregate = TRUE, 
         save_path = "paper/Figures/fdr_vs_study_w_nonid.pdf", 
         use_non_identified = TRUE)
# plot_sigval_for_low_fwer()
