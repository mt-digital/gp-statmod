library(ggplot2)
library(readr)
library(dplyr)

source("model.R")  # for binProb function used in ordinal plot function below.


# mytheme = theme(
#   axis.line = element_line(), legend.key=element_rect(fill = NA),
#   text = element_text(size=16),# family = 'PT Sans'),
#   # legend.key.width = unit(2, 'cm'),
#   # legend.key.size = unit(1.5, 'lines'),
#   panel.background = element_rect(fill = "white"),
#   panel.grid.major.x = element_line(color = "lightgrey", linewidth = 0.1, linetype = 2)
# )
  mytheme = theme(
    panel.border = element_blank(), axis.line = element_line(), 
    text = element_text(size=16, family = 'PT Sans'), 
    panel.background = element_rect(fill = "white")
  )

plot_latent_pdf_integration = function(mu = 0, sd = 2.5, min_bin, max_bin, 
                                       bins, bin_colors, # xlim = c(-6.0, 6.0), 
                                       save_path = "test_latent_pdf.pdf", 
                                       max_prob = 1.0) {
  print(save_path)


  xmin = min_bin - 0.5
  xmax = max_bin + 0.5

  xlim = c(xmin, xmax)

  print(xlim)

  p = ggplot(data.frame(x = xlim), aes(x)) 

  p = p + stat_function(fun = dnorm,
                        args = c(mu, sd),
                        geom = "line",
                        xlim = c(xmin, xmax)) 

  # First bin from -infinity to min_bin + 0.5
  p = p + stat_function(fun = dnorm,
                        geom = "area",
                        args = c(mu, sd),
                        fill = bin_colors[1],
                        xlim = c(xmin, min_bin + 0.5))

  n_bins = length(bins)

  # Middle bins.
  for (idx in 2:(n_bins - 1)) {

    p = p + stat_function(fun = dnorm,
                          geom = "area",
                          args = c(mu, sd),
                          fill = bin_colors[idx],
                          xlim = c(bins[idx - 1] + 0.5, bins[idx] + 0.5))
  }

  # Final bin from last_bin - 0.5 to infinity.
  p = p + stat_function(fun = dnorm,
                        geom = "area",
                        args = c(mu, sd),
                        fill = bin_colors[n_bins],
                        xlim = c(xmax, max_bin - 0.5))

  p = p + xlab("Latent opinion") + ylab("Probability density\n\n") + 
  geom_vline(xintercept=c((bins[1:(length(bins) - 1)] + 0.5)), linetype="dotted") + 
  geom_vline(xintercept=mu, linetype="dashed") + 
  xlim(c(xmin, xmax)) + 
  ylim(c(0, max_prob)) +
  scale_x_continuous(breaks = bins) +
  mytheme
  
  ggsave(filename = save_path, device = cairo_pdf, p, width=5, height=3.5, units="in")
  
  return (p)
}


plot_ordinal_distribution = function(mu = 0, sd = 2.5, min_bin = -2, 
                                     max_bin = 2, bins = -2:2,
                                     bin_colors = c("lightpink", "lightblue", 
                                                      "#f7ae3d", "lightgreen", 
                                                      "#F6CFFF"),
                                     save_path = "test_ordinal_pdf.pdf",
                                     max_prob = 1.0) {
  
  # Use same parameters as latent distribution to calculate integration over
  # bin threshold limits \theta_{k-1}..\theta_k
  opinion_bin_vec <- seq(min_bin, max_bin)

  n_bins <- length(opinion_bin_vec)
  
  ord_dist_df <- data.frame(OpinionBin = seq(min_bin, max_bin), 
                            Density = rep(0, n_bins))
  
  K <- max_bin - min_bin + 1
  k_idxs <- bins - min_bin + 1
  
  for (ii in 1:length(k_idxs)) {
    bprob <- binProb(bins[ii], min_bin, K, mu, sd)
    ord_dist_df[ii, "Density"] <- bprob
  }
  
  mean_observed <- sum(ord_dist_df$OpinionBin * ord_dist_df$Density)
  
  print(paste("Mean observed for sd = ", sd, ": ", mean_observed, sep = ""))
  
  # Set x-axis limits to match continuous distro.
  xmin = min_bin - 0.5
  xmax = max_bin + 0.5

  xlim = c(xmin, xmax)

  print(xlim)

  p <- ggplot(data=ord_dist_df, aes(x=OpinionBin, y=Density)) +
         geom_bar(stat = "identity", fill = bin_colors) + 
         geom_vline(xintercept = mean_observed, linetype="dashed") +
         xlab("Ordinal opinion measurement") + ylab("Probability density\n\n") +
         ylim(c(0.0, max_prob)) + 
         scale_x_continuous(breaks = bins, limits = xlim) +
         mytheme
  
  ggsave(filename = save_path, device = cairo_pdf, p, 
         width=5, height=3.5, units="in")
  
  return (p)
}


# Consensus occurs when group opinion variance decreases; simple consensus when mu 
# remains constant, while group polarization is when the mean becomes more extreme.
make_distros_consensus_figure =
  function(mu = 0.8, sigma_pre = 2.0, sigma_post = 1.0, min_bin = -2, 
           max_bin = 2, bins = -2:2, bin_colors = c("lightpink", "lightblue", 
                                                    "#f7ae3d", "lightgreen", 
                                                    "#F6CFFF"),
           save_dir = 
             "~/workspace/gp-statmod/GroupPolarizationStatmod/paper/Figures/Model"
           , max_prob = 1.0) {

  sds = c(sigma_pre, sigma_post)
  for (sd in sds) {
    save_path = fs::path_join(c(save_dir, paste0("ordinal", "_mu=", mu, "_sd=", round(sd, 2), ".pdf")))

    plot_ordinal_distribution(mu, sd, min_bin, max_bin, bins, bin_colors, save_path = save_path, max_prob = max_prob)
  }
  
  for (sd in sds) {
    save_path = fs::path_join(c(save_dir, paste0("latent", "_mu=", mu, "_sd=", round(sd, 2), ".pdf")))

    print(save_path)

    plot_latent_pdf_integration(mu, sd, min_bin, max_bin, bins, bin_colors, 
                                save_path = save_path, max_prob = max_prob)
  }
  
}


# Run figure-making routines using data from CaseStudies.xlsx
# Get schkade2010 data
case_studies_tbl = read_csv("data/StudiesAnalysis.csv")

# XXX note that the analyzed study CSV uses column name "ArticleTag" not StudyID
# and "TreatmentTag" not ExperimentID. 
# TODO Not sure right now where this change is done, need to find out
# schkade2010 = filter(case_studies_tbl, ArticleTag == "Schkade2010" &
#                      TreatmentTag == "COSprings-CivilUnions")
                         
# mu = schkade2010$LatentMean
# sigma_pre = schkade2010$LatentSDPre
# sigma_post = schkade2010$LatentSDPost
# min_bin = schkade2010$MinBinValue
# max_bin = schkade2010$MaxBinValue
# bins = min_bin:max_bin


# Trying Moscovici "Americans" as example to match other model figure.
moscovici1969 = filter(case_studies_tbl, 
                     ArticleTag == "Moscovici1969" & 
                       TreatmentTag == "Americans")
                         
mu = moscovici1969$LatentMean
sigma_pre = moscovici1969$LatentSDPre
sigma_post = moscovici1969$LatentSDPost
min_bin = moscovici1969$MinBinValue
max_bin = moscovici1969$MaxBinValue
bins = min_bin:max_bin

rhg_cols_10 <- rev(c("#771C19", "#AA3929", "#E25033", "#F27314", "#F8A31B", 
              "#E2C59F", "#B6C5CC", "#8E9CA3", "#556670", "#556AAA"))

bin_colors = rhg_cols_10[1:7]

make_distros_consensus_figure(mu, sigma_pre, sigma_post, 
                              min_bin, max_bin, bins, bin_colors, max_prob = 0.4)
