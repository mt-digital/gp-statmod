# Group Polarization Generative Modeling

## Background

This code implements numerical analyses, generative simulations and statistical
tests described in the paper "Many claims of group polarization are plausibly 
false due to inappropriate statistics". 

## Code

### Simulate simple consensus group polarization data

Given a static latent mean, a pre-deliberation and post-deliberation
standard deviation, we can generate simulated simple consensus data to be fit
using either a metric- or ordinal-data statistical model. The metric model is
a t-test, which is equivalent to a linear model. The ordinal model is the 
ordered probit model, which accounts for the binning process that experiment
participants do when reporting their ordinal opinions.

### Metric linear model of group polarization

We fit a linear model of the simulated simple consensus data drawn from
continous normal pre- and post-deliberation distributions and put in ordinal
bins. For simplicity we use the `t.test` function in R, which may be also be a 
common approach used for this type of inference. Specifically, in the
`simulate_metric_cohens_d` function in the file
[`experiments.R`](https://github.com/mt-digital/gp-statmod/blob/main/GroupPolarizationStatmod/experiments.R)
there is the following line that fits a `t.test` to simulated observations of
pre- and post-deliberation opinions

```R
t_result <- t.test(sim_obs_pre, sim_obs_post, paired = paired, var.equal = var.equal)
```

The linear estimates of means and standard deviations are then extracted in the
following lines:

```R
pre_mean_estimate <- t_result$estimate[[1]] 
post_mean_estimate <- t_result$estimate[[2]] 

pre_sd <- sd(sim_obs_pre)
post_sd <- sd(sim_obs_post)
```

Effect size is quantified as Cohen's _d_.  We calculate _d_ in 1000 trials using `make_metric_cohens_d_table` in
[`analysis.R`](https://github.com/mt-digital/gp-statmod/blob/main/GroupPolarizationStatmod/analysis.R).
The output is written to a csv file specified by the `output_file` argument.
To analyze the distribution of _d_ values we create a [boxplot using
`ggplot2`](http://www.sthda.com/english/wiki/ggplot2-box-plot-quick-start-guide-r-software-and-data-visualization).

### Ordered probit model of group polariztion


