# Generative modeling of group polarization experiments to control false detection rate

## Code

### Identify plausibly false detections

### Bayesian fits

( explain how `singleBayesianFitTrial` works in [`experiments.R`](/experiments.R) )


### Cluster script for Bayesian fits

[`bayesian_fit_trial.sh`](/bayesian_fit_trial.sh) uses
`singleBayesianFitTrial` [`in experiments.R`](/experiments.R)

Submit a job array to run 1000 trials like so

```
qsub --array=[1-1000] bayesian_fit_trial.sh
```

Each of the 1000 trials writes its output data 
to a randomly-named .csv into the `data/probit_fits`
directory, which must be created before submitting the job array.
Each .csv output has 55 rows: 54 rows of data, one for each  and one header row.


### False detection rate calculations

Cohen's *d* is calculated with `cohens_d` in
[`src/false_detection_rate.R`](/src/false_experiments.R). 
This asymmetric-variance version of Cohen's d is enumerated in Liddell and
Kruschke (2018) and used by them at line 884 in their script
[OrdinalScaleGroupJags.R](https://osf.io/5jrgz), hosted in the associated 
["Ordinal Data Analysis" OSF repository](https://osf.io/53ce9/).

XXXX EDIT XXXX For each plausible false detection experimental condition we calculate $\alpha =
\Pr(D|\not E)$ for a given significance value, $\alpha^*$, by calculating the
fraction of $d$ values that are greater than $\alpha^*$ (EVENTUALLY/SOON)
to be implemented in [`false_detection_rate.R`](/false_detection_rate.R).



See `src/false_detection_rate.R`, used by `scripts/analysis.R`, both used by
`scripts/plot.R`. The analysis of the probit fits is entirely here now, I believe.
Still need a lot of the other stuff, which should be cleaned up and organized
into `src` or `scripts` if needed, or discarded if no longer used.

