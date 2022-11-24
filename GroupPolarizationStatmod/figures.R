library(ggplot2)

plot_latent_pdf_integration = function() {

  mytheme = theme(
    panel.border = element_blank(), axis.line = element_line(), 
    text = element_text(size=16, family = 'PT Sans'), 
    panel.background = element_rect(fill = "white")
  )

  p = ggplot(data.frame(x = c(-3, 3)), aes(x)) + 

    stat_function(fun = dnorm,
                  args = c(0, 2.5),
                  geom = "line",
                  xlim = c(-5, 5)) +
    stat_function(fun = dnorm,
                  geom = "area",
                  args = c(0, 2.5),
                  fill = "steelblue",
                  xlim = c(0.5, 1.5)) +
    stat_function(fun = dnorm,
                  geom = "area",
                  args = c(0, 2.5),
                  fill = "steelblue",
                  xlim = c(2.5, 5)) +

    xlab("Latent opnion") + ylab("Probability density") +

    mytheme
  
  return (p)
}


plot_ordinal_distribution = function() {
  # Use same parameters as latent distribution to calculate integration over
  # bin threshold limits \theta_{k-1}..\theta_k
}
