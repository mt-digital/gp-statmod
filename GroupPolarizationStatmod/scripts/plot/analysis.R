##
#
# Plot analysis of family-wise error rate simulations.
# 
# Currently model illustrations are in the root project directory and will
# soon be ported here.
# 
# Author: Matthew A. Turner <maturner01@gmail.com>
# Date: 2024-12-07
#

library(forcats)
library(ggplot2)
library(latex2exp)
library(emojifont)

source("src/fwer.R")
source("scripts/analysis.R")

mytheme = theme(
  axis.line = element_line(), legend.key=element_rect(fill = NA),
  text = element_text(size=18),# family = 'PT Sans'),
  # legend.key.width = unit(2, 'cm'),
  # legend.key.size = unit(1.5, 'lines'),
  panel.background = element_rect(fill = "white"),
  panel.grid.major.x = element_line(color = "lightgrey", linewidth = 0.1, linetype = 2)
)

rhg_cols <- rev(c("#771C19", "#AA3929", "#E25033", "#F27314", "#F8A31B", 
              "#E2C59F", "#B6C5CC", "#8E9CA3", "#556670", "#556AAA", "#000000"))

rhg_cols_10 <- rev(c("#771C19", "#AA3929", "#E25033", "#F27314", "#F8A31B", 
              "#E2C59F", "#B6C5CC", "#8E9CA3", "#556670", "#556AAA"))

plot_ordinal_cohens <- function(ordinal_data_dir = "data/probit_fits",
                                sync_file = "data/probit_fits/all.csv",
                                output_file = 
                                  "paper/Figures/Analysis/ordinal_cohens.pdf",
                                overwrite = FALSE) {
  
  df <- load_probit_data(probit_fits_dir, sync_file, overwrite)
  
  df$Cohens_d <- cohens_d(df$LatentMeanPrePosteriorMean, 
                          df$LatentMeanPostPosteriorMean,
                          df$LatentMeanPrePosteriorSD, 
                          df$LatentMeanPostPosteriorSD)
  
  ggplot(df, mapping = aes(x = Cohens_d, y = ExperimentID)) + 
    geom_vline(xintercept = 0, color = "red") + 
    geom_vline(xintercept = c(-1, 1), color="red", linetype="dotted") + 
    geom_boxplot() 
  
  ggsave(output_file)
}




plot_fdr = function(save_path = "paper/Figures/Analysis/fdr_vs_experiment.pdf", 
                    aggregate = FALSE, use_non_identified = FALSE,
                    default_fwer = 0.05,
                    base_rate = 0.1, power = 0.8, 
                    significance_vals = c(0.2, 0.5, 0.8),
                    focal_sig_val = 0.5) {

  # Load data.
  tbl = calc_fdr_vs_significance(significance_vals = significance_vals) %>%
    mutate(ExperimentID = paste(StudyID, ExperimentID, sep = "_")) %>%
    mutate(SigVal = as_factor(SigVal), Power = power, BaseRate = base_rate)

  if (aggregate)
    p = plot_fdr_aggregated(tbl, focal_sig_val, save_path, default_fwer, use_non_identified)
  else
    p = plot_fdr_nonaggregated(tbl, focal_sig_val, save_path)

  return (p)
}


plot_fdr_aggregated = function(tbl, focal_sig_val, 
                               save_path, default_fwer = 0.05,
                               use_non_identified = FALSE) {
  
  tbl = 
    as_tibble(aggregate_fdr_vs_sig(tbl, use_non_identified, default_fwer))

  tbl_foc = 
    filter(tbl, SigVal == focal_sig_val) %>%
      select(!SigVal) %>%
      mutate(StudyID = fct_reorder(StudyID, FWER))

  cat(
    paste0("\n\n***** CALCULATED MEAN FDR = ", 
           mean(tbl_foc$FDR), 
           ifelse(use_non_identified, " with", " without"),
           " default FWER", 
           ifelse(use_non_identified, paste0(" = ", default_fwer), ""),
           " *****\n\n"
    )
  )

  tbl_all_foc = tibble(StudyID = as.factor(c("All")), 
                       FWER = c(mean(tbl_foc$FWER)), 
                       FDR = c(mean(tbl_foc$FDR)))

  tbl_all = 
    group_by(tbl, SigVal) %>% 
      summarize(StudyID = as.factor("All"), FWER = mean(FWER), FDR = mean(FDR))

  tbl_foc = rbind(tbl_all_foc, tbl_foc)
  tbl = rbind(tbl_all, tbl)

  p = 
    ggplot(tbl_foc, aes(x = FWER, y = StudyID, fill = StudyID)) +
      geom_hline(yintercept = 1.5, linetype = "dashed", linewidth = 0.25) +
      geom_bar(aes(fill = StudyID, color = StudyID), stat = "identity", width = 0.1) +
      geom_point(aes(x = FWER, y = StudyID, shape = SigVal, color = StudyID), 
                 tbl, size = 6) + 
      geom_point(aes(x = FDR, y = StudyID, shape = SigVal, color = StudyID), 
                 fill = "transparent", 
                 tbl, size = 6) + 

      xlim(c(0, 1)) +
      xlab("Family-wise error rate (\U25CF); false discovery rate (\U25CB)") + 
      ylab("Study ID") +
      labs(shape = TeX("Significance, $d^*$")) +
      labs(fill = "Study ID", color = "Study ID") +
      scale_shape_manual(values = c(24, 22, 21), 
                         labels = c("Low, 0.2", "Medium, 0.5", "High, 0.8"))  +
      guides(shape = guide_legend(reverse = TRUE), 
             color = guide_legend(reverse = TRUE), 
             fill = guide_legend(reverse = TRUE))  +
      scale_fill_manual(values = rhg_cols) + 
      scale_color_manual(values = rhg_cols) + 
      mytheme

  ggsave(save_path, p, width = 9.5, height = 10)

  return (tbl_all)
}


plot_fdr_nonaggregated = function(tbl, focal_sig_val, save_path) {
  
  # Filter to keep just "medium" effect sizes and sort factors.
  tbl_foc = 
    filter(tbl, SigVal == focal_sig_val) %>% 
    mutate(StudyID = fct_reorder(StudyID, FWER, mean), 
           ExperimentID = fct_reorder(ExperimentID, FWER, max))

  p = 
    ggplot(tbl_foc, aes(x = FWER, y = ExperimentID, fill = StudyID)) + 
      geom_bar(aes(fill = StudyID, color = StudyID), stat = "identity", width = 0.15)  +
      labs(fill = "Study ID", color = "Study ID") +
      geom_point(aes(x = FWER, y = ExperimentID, shape = SigVal, color = StudyID), tbl, size = 3.35) + 
      geom_point(aes(x = FDR, y = ExperimentID, shape = SigVal, color = StudyID), fill = "transparent", 
                 tbl, size = 3.35) + 
      labs(shape = TeX('Significance, $d^*$')) +
      xlim(c(0, 1)) + 
      # xlab(TeX("Family-wise error rate ($\\bullet $); false discovery rate ($\\omicron$)")) + 
      xlab("Family-wise error rate (\U25CF), False discovery rate (\U25CB)") + 
      ylab("Experiment ID") +
      scale_shape_manual(values = c(24, 22, 21), labels = c("Low, 0.2", "Medium, 0.5", "High, 0.8")) +
      guides(shape = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE), 
             fill = guide_legend(reverse = TRUE)) +
      scale_fill_manual(values = rhg_cols_10) + scale_color_manual(values = rhg_cols_10) + 
      mytheme

  ggsave(save_path, p, width = 9.5, height = 10)

  return (p)
}


plot_sigval_for_low_fwer = function(target_fwer = 0.05, 
                                    save_path = "paper/Figures/Analysis/sigval_for_low_fwer.pdf",
                                    probit_fits_dir = "data/probit_fits",
                                    focal_sig_val = 0.5,
                                    base_rate = 0.1, power = 0.8, 
                                    significance_vals = seq(0.5, 3.0, 0.01)) {
  
  # Going to create ggplot, p.
  # First load and pre-process calculation of FWER ≤ target_fwer by experiment.
  pre_tbl = calc_fdr_vs_significance(significance_vals = significance_vals)
  pre_tbl_focal = filter(pre_tbl, SigVal == focal_sig_val) %>%
    mutate(StudyID = fct_reorder(StudyID, FWER, mean),
           ExperimentID = fct_reorder(ExperimentID, FWER, max))

  print(pre_tbl_focal)
  study_id_fct = levels(pre_tbl_focal$StudyID)

  tbl = sigval_for_low_fwer(pre_tbl) %>%
    mutate(ExperimentID = as.factor(paste(StudyID, ExperimentID, sep = "_")),
           SigVal = as.numeric(as.character(SigVal))) 
  
  # Then plot the result.
  p = ggplot(tbl, aes(x = SigVal, y = ExperimentID, fill = StudyID, color = StudyID)) +
      geom_vline(xintercept = c(0.2, 0.5, 0.8), linetype = "dashed", 
                 linewidth = 0.25, color = "darkgrey") + 
      geom_bar(aes(fill = StudyID, color = StudyID), 
               stat = "identity", width = 0.15)  +
      geom_point(aes(fill = StudyID, color = StudyID), shape = 18, size = 5) +
      xlab(TeX("Significance $d^*$ for $\\alpha \\leq 0.05")) +
      ylab("Experiment ID") +
      scale_x_continuous(breaks = c(0.2, 0.5, 0.8, 1.1, 1.5, 2.0, 2.5)) +
      scale_fill_manual(breaks = study_id_fct, values = rhg_cols_10) +
      scale_color_manual(breaks = study_id_fct, values = rhg_cols_10) +
      guides(shape = guide_legend(reverse = TRUE), 
             color = guide_legend(reverse = TRUE),
             fill = guide_legend(reverse = TRUE)) +
      mytheme
        
   ggsave(save_path, p, width = 9.5, height = 10)

   return (p)
}


## MAKE MAIN PUBLICATION PLOTS ##
# plot_ordinal_cohens()

plot_fdr()

plot_fdr(aggregate = TRUE, save_path = "paper/Figures/Analysis/fdr_vs_study.pdf")

plot_fdr(aggregate = TRUE,
         save_path = "paper/Figures/Analysis/fdr_vs_study_w_nonid.pdf",
         use_non_identified = TRUE)

plot_sigval_for_low_fwer()

