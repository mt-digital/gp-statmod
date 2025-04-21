---
editor_options: 
  markdown: 
    wrap: 80
---

# Group polarization replications may be marred by high false discovery rates

This repository provides code to analyze whether published replications of
*group polarization*, where groups of people theoretically become more extreme
in their opinions following deliberation with like-minded others. This code is
somewhat complex because of the needs of the project, which have evolved over
the past few years. The project has accumulated some technical debt that will
need to be paid up one day. For now, please excuse a bit of disorganization in
the code. Hopefully the notes here are sufficient to find functions performing
different parts of the data input and analysis pipeline.

The ultimate goal of the this code and the project is to identify which
published detections of group polarization are *plausibly spurious*, meaning
that we have identified opinion distributions for before and after deliberation
where average group extremism in terms of latent psychological opinions does not
change, but when before and after opinions are measured on an ordinal scale
(e.g., a Likert scale), the opinions appear to have shifted to become more
extreme. This is caused by ceiling effects induced by the ordinal scale that may
not capture extreme opinions becoming more centrist, but does capture centrist
opinions becoming more extreme. We use the term *simple consensus* for the case
where opinion variance in the group decreases but the mean opinion does not
change, to distinguish this case from *group polarization*, also a form of
consensus, but one in which the mean opinion becomes more extreme following
deliberation.

Here is an outline of the components gathered in this repository, with
references to figures, tables, analyses, etc., used in our analysis and
presentation:

1.  **Identify of plausibly spurious group polarization detections.**

    -   Shiny web app for data entry, researcher input of initial guesses, and
        evaluation of parameters identified to induce spurious group
        polarization in the *distributional model*, that identifies plausible
        spurious group polarization by integrating latent opinion distributions
        over ordinal scale bins.

    -   To identify plausibly spurious parameters, the researcher must first
        input identifying information for the research article (the "study"),
        the experimental condition within the study, the observed pre- and
        post-deliberation mean opinions, and a guess for what constant latent
        opinion mean and pre- and post-deliberation variances could generate the
        reported observations.

    -   The researcher then can run a hillclimbing-style algorithm
        (`solveForLatentSD` function in [`numerical.R`](numerical.R)) that
        identifies simple consensus latent distribution parameters that generate
        spurious group polarization observations. The researcher must inspect
        these results before deciding if it should count as a plausibly spurious
        detection. To finalize the analysis of an experimental condition the
        researcher must indicate that it should be included as a plausibly
        spurious detection of group polarization.

2.  **Use the simple consensus parameters found for all plausibly spurious
    detections to seed generative simulations of group polarization experiments
    to which ordered probit models are fit to estimate best-case scenario
    family-wise error rates and false discovery rates.**

    -    At the top level for generating the analysis here, see
        [`bayesian_fit_trial.sh`](/bayesian_fit_trial.sh) , a slurm submission
        script that requests resources and runs the function
        `singleBayesianFitTrial` [`in experiments.R`](/experiments.R).

    -   `singleBayesianFitTrial` in turn sets up a new write file in the output
        directory and runs `makeBayesianFitTable`, also in
        [experiments.R](experiments.R), which encapsulates the following steps:

        1.  Read latent simple consensus parameters that generated spurious
            group polarization from the file
            [data/StudiesAnalysis.csv](data/StudiesAnalysis.csv), which stores
            the identified parameters, and other data, using the web app as
            explained above.

        2.  Create simulated observed opinions pre- and post-deliberation by
            drawing latent opinions for $N$ simulated participants from the
            spurious group polarization distributions, then binning them as if
            reported by simulated participants on the same ordinal scale and $N$
            as used in the real-world experiments (`makeJAGSModelData` in
            [`experiments.R`](experiments.R)).

        3.  A normal continuous distribution is inferred from MCMC fitting of an
            ordered probit model (done in `calculateBayesian` in
            [`experiments.R`](experiments.R)) for both pre- and
            post-deliberation latent opinion distributions.

3.  **Estimate best-case scenario family-wise error rates (FWERs) and false
    discovery rates (FDRs) given ordered probit inferences from simulated
    datasets.**

4.  **Present results in a table and plots for individual experiments and
    aggregated to the study level and over all experiments in all studies.**

    -   Quantiles table made using source("scripts/calculations.R");
        \`quantiles_tbl\`

5.  **Make illustration of the large-N model of spurious group polarization from
    simple consensus.**

    -   

### Identify plausibly false detections

### Bayesian fits

( explain how `singleBayesianFitTrial` works in
[`experiments.R`](/experiments.R) )

### Cluster script for Bayesian fits

Submit a job array to run 1000 trials like so

```         
qsub --array=[1-1000] bayesian_fit_trial.sh
```

The `bayesian_fit_trial.sh` script sets up cluster parameters, requesting 4 CPUs
with 8 GB of memory for each trial.

Each of the 1000 trials writes its output data to a randomly-named .csv into the
`data/probit_fits` directory, which must be created before submitting the job
array. Each .csv output has 55 rows: 54 rows of data, one for each of the
experimental conditions indexed by $e$ in the paper, and one header row.

### False detection rate calculations

Cohen's *d* is calculated with `cohens_d` in
[`src/false_detection_rate.R`](src/false_detection_rate.R). This
asymmetric-variance version of Cohen's d is enumerated in Liddell and Kruschke
(2018) and used by them at line 884 in their script
[OrdinalScaleGroupJags.R](https://osf.io/5jrgz), hosted in the associated
["Ordinal Data Analysis" OSF repository](https://osf.io/53ce9/).

XXXX EDIT XXXX For each plausible false detection experimental condition we
calculate $\alpha =
\Pr(D|\neg E)$ for a given significance value, $\alpha^*$, by calculating the
fraction of $d$ values that are greater than $\alpha^*$ (EVENTUALLY/SOON) to be
implemented in [`false_detection_rate.R`](/false_detection_rate.R).

See `src/false_detection_rate.R`, used by `scripts/analysis.R`, both used by
`scripts/plot.R`. The analysis of the probit fits is entirely here now, I
believe. Still need a lot of the other stuff, which should be cleaned up and
organized into `src` or `scripts` if needed, or discarded if no longer used.
