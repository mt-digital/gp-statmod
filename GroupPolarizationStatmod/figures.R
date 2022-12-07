library(ggplot2)

plot_latent_pdf_integration = function(mu = 0, sd = 2.5, xlim = c(-6.0, 6.0), 
                                       bin1 = c(-0.5, 0.5), bin2 = c(1.5, 5.0),
                                       save_path = "test_latent_pdf.pdf") {

  mytheme = theme(
    panel.border = element_blank(), axis.line = element_line(), 
    text = element_text(size=16, family = 'PT Sans'), 
    panel.background = element_rect(fill = "white")
  )

  p = ggplot(data.frame(x = xlim), aes(x)) + 

    stat_function(fun = dnorm,
                  args = c(mu, sd),
                  geom = "line",
                  xlim = xlim + 1) +
    
    stat_function(fun = dnorm,
                  geom = "area",
                  args = c(mu, sd),
                  fill = "lightpink",
                  xlim = c(-5.0, -1.5)) +
    
    stat_function(fun = dnorm,
                  geom = "area",
                  args = c(mu, sd),
                  fill = "lightblue",
                  xlim = c(-1.5, -0.5)) +

    stat_function(fun = dnorm,
                  geom = "area",
                  args = c(mu, sd),
                  fill = "#ffcb58",
                  xlim = c(-0.5, 0.5)) +

    stat_function(fun = dnorm,
                  geom = "area",
                  args = c(mu, sd),
                  fill = "lightgreen",
                  xlim = c(0.5, 1.5)) +
    
    stat_function(fun = dnorm,
                  geom = "area",
                  args = c(mu, sd),
                  fill = "#F6CFFF",
                  xlim = c(1.5, 5)) +

    xlab("Latent opinion") + ylab("Probability density\n\n") + 
    geom_vline(xintercept=-1.5:1.5, linetype="dotted") + 
    geom_vline(xintercept=mu, linetype="dashed") + 
    xlim(c(-5, 5)) + 
    # ylim(0.01, NA) + 

    mytheme
  
  ggsave(filename = save_path, device = cairo_pdf, p, width=5, height=3.5, units="in")
  
  return (p)
}


plot_ordinal_distribution = function(mu = 0, sd = 2.5, min_bin = -2, 
                                     max_bin = 2, bins = -2:2,
                                     bin_colors = c("lightpink", "lightblue", 
                                                      "#f7ae3d", "lightgreen", 
                                                      "#F6CFFF"),
                                     save_path = "test_ordinal_pdf.pdf") {
  
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
  
  p <- ggplot(data=ord_dist_df, aes(x=OpinionBin, y=Density)) +
         geom_bar(stat = "identity", fill = bin_colors) + 
         geom_vline(xintercept = mean_observed, linetype="dashed") +
         xlab("Binned ordinal opinion") + ylab("Probability density\n\n") +
         mytheme
  
  ggsave(filename = save_path, device = cairo_pdf, p, 
         width=5, height=3.5, units="in")
  
  return (p)
}


# Consensus occurs when group opinion variance decreases; simple consensus when mu 
# remains constant, while group polarization is when the mean becomes more extreme.
make_distros_consensus_figure <- function(mu = 0.8, sds = c(2.0, 1.0), 
                                          save_dir = 
  "~/workspace/gp-statmod/GroupPolarizationStatmod/paper/ConsensusDistroIllustration/") {

  for (sd in sds) {
    save_path = paste(save_dir, "ordinal", "_mu=", mu, "_sd=", sd, ".pdf", 
                      sep = "")

    plot_ordinal_distribution(mu, sd, save_path = save_path)
  }
  
  for (sd in sds) {
    save_path = paste(save_dir, "latent", "_mu=", mu, "_sd=", sd, ".pdf", 
                      sep = "")

    plot_latent_pdf_integration(mu, sd, save_path = save_path)
  }
  
}

