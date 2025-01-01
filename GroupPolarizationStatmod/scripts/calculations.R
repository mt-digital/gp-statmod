##
# Calculations that are included in paper prose.
#
# Author: Matthew A. Turner
# Date: 2024-12-27
#


source("scripts/analysis.R")


# Create a CSV list of quantile calculations that may be included in publications.
quantiles_tbl = function(tex_output_file = "paper/tex/quantile_table.tex", 
                         q_resolution_or_probs = 0.05, 
                         significance_vals = c(0.2, 0.5, 0.8),
                         base_rate = 0.1, power = 0.8
                         ) {
  
  # Hacky but just need this to proceed for now.
  if (length(q_resolution_or_probs) > 1) {
    probs <- q_resolution_or_probs
  } else {
    probs <- seq(0.0, 1.0, q_resolution_or_probs)
  }

  initialized = FALSE
  for (sigval in significance_vals) {

    tbl <- calc_fdr_vs_significance(significance_vals = sigval, 
                                    base_rate = base_rate, power = power)

    for (measure in c("FWER", "FDR")) {

      # Get quantile for this measure, base rate, and power.
      q <- quantile(as_vector(tbl[, measure]), probs = probs)

      if (!initialized) {
        # Load quantile named list into tbl with quantiles as colnames sans "%".
        q_tbl <- tibble(!!!q, .name_repair = ~ str_replace(., "%", ""))
        q_tbl$BaseRate <- base_rate
        q_tbl$Power <- power
        q_tbl$SigVal <- sigval
        q_tbl$Measure <- measure

        initialized = TRUE
      }
      else {
        # Create new row.
        new_row <- tibble(!!!q, .name_repair = ~ str_replace(., "%", ""))
        new_row$BaseRate <- base_rate
        new_row$Power <- power
        new_row$SigVal <- sigval
        new_row$Measure <- measure

        # Bind new row to existing.
        q_tbl <- rbind(q_tbl, new_row)
      }
    }
  }

  # Reorder columns to be easier to read.
  q_tbl <- 
    q_tbl %>% 
      relocate(Measure, .before = `0`) %>% 
      relocate(SigVal, .after = Measure) %>% 
      relocate(BaseRate, .after = SigVal) %>% 
      relocate(Power, .after = BaseRate) %>%
      arrange(desc(Measure))

  # q_tbl_lim <- filter(SigVal == 0.5) %>% select(!c("SigVal", 
  q_tbl_lim <- q_tbl %>% select(!c("BaseRate", "Power"))
  
  if (!is_null(tex_output_file)) {
    print(xtable(q_tbl_lim, type = "latex"), file = tex_output_file, 
          include.rownames = FALSE, booktabs = TRUE, floating = FALSE)
  }

  return (q_tbl)
}


build_quantiles_across_contexts <- 
  function(pessimistic_b_W = c(0.05, 0.5), 
           middle_b_W = c(0.1, 0.8),
           optimistic_b_W = c(0.5, 0.95),
           quantile_probs = as.vector(c(0.0, 0.2, 0.5, 0.8, 1.0))) {

  pbw <- pessimistic_b_W
  mbw <- middle_b_W
  obw <- optimistic_b_W

  qt_p <- quantiles_tbl(quantile_probs, 
                                  base_rate = pbw[1], power = pbw[2],
                                  tex_output_file = NULL)

  qt_m <- quantiles_tbl(quantile_probs, 
                             base_rate = mbw[1], power = mbw[2],
                             tex_output_file = NULL)

  qt_o <- quantiles_tbl(quantile_probs, 
                                 base_rate = obw[1], power = obw[2],
                                 tex_output_file = NULL)

  # FWER is the same in each of the above.
  qt_fwer <- qt_p %>% filter(Measure == "FWER")

  qt_fdr = rbind(qt_p, qt_m, qt_o) %>% filter(Measure == "FDR")

  return (list(fwer = qt_fwer, fdr = qt_fdr))
}

