# GP-stats: Spurious Group Polarization Analysis

This project implements the modeling, simulation, and analysis pipeline described in:

> **Turner & Smaldino (2025)**  
> *If the Null Fits, You Must Omit: Ubiquitous False Detections of Group Polarization*

---

## Overview

GP-stats estimates the **false detection rate (FDR)** of group polarization findings by simulating **null-consistent latent opinion distributions** measured on ordinal (Likert) scales. It includes:

- A **Shiny web app** for user-guided parameter fitting
- R scripts for **large-N analytical models**, **hillclimbing search**, **Bayesian ordered probit fitting**, and **parallel simulations**
- Tools to **aggregate and visualize** error rates across experiments and studies

---

## Code ↔ Paper Mapping

| Paper Section | Description | Implementation |
|---------------|-------------|----------------|
| **2.2 – How study design can induce false detections** | Latent–ordinal mapping, bin thresholds (Eq. 2–6), observed mean computation | `model.R` (prob. integrals), `numerical.R` (`meanObs`, `sdObs`) |
| **2.3 – Web app & hillclimbing algorithm** | Find σ_pre, σ_post s.t. μ constant and simulated means match empirical | `numerical.R` (`solveForLatentSD`), `experiments.R` (wrapper logic), `app.R` (UI), `util.R` (UI helpers) |
| **2.4 – Estimating & limiting FDR** | Simulate experiments using fitted β, fit ordered probit, compute FWER & FDR (Eq. 7–10) | `experiments.R` (`makeBayesianFitTable`), `analysis.R` (loop & aggregation), `model.R` (effect size calc) |
| **3.1 – Plausibly spurious detections** | Mark conditions after visual sanity-check in app | `app.R` / `util.R` interactive marking, results saved to CSV |
| **3.2 – FDR distribution** | Quantiles & summary tables (Table 1) | `analysis.R` (`makePlausibleFPTable`) |
| **3.3 – Required d* thresholds** | Inverse problem: find d* for FWER ≤ 0.05 | `experiments.R` (threshold scan), plotted in `plot.R` |

---

## Workflow Diagram

The diagram below shows the GP-stats data flow from empirical inputs, through the Shiny app and batch simulations, to statistical outputs and figures.

![](gp-stats-workflow.png){width=500}

---

## Data Flow

1. **Load empirical data**  
   `data/StudiesAnalysis.csv` → `analysis.R`  
   Contains pre/post means, N, and bin settings per experiment.

2. **Parameter fitting**  
   Shiny app or batch mode hillclimbing finds β = {μ, σ_pre, σ_post}.

3. **Simulation & fitting**  
   For each β, simulate `N_emp` trials → fit Bayesian ordered probit → compute d, α, and FDR.

4. **Aggregation & output**  
   Study-level and corpus-level summaries; saved under `data/output/`.

5. **Visualization**  
   `plot.R` renders figures for FWER/FDR distributions and d* thresholds.

---

## Requirements

- **R ≥ 4.0**
- Packages: `runjags`, `rjags`, `dplyr`, `tibble`, `parallel`, `uuid`, `shiny`
- **JAGS** installed and on PATH

---

## Quick Start

```r
# Launch web app for parameter fitting
shiny::runApp("app.R")

# Run full simulation analysis
source("analysis.R")
makePlausibleFPTable()

# Plot results
source("plot.R")
```


## License

Creative Commons Attribution 4.0 International (CC BY 4.0)
You are free to share and adapt this work for any purpose, even commercially, provided that you give appropriate credit to the original authors:

> Turner, M.A. & Smaldino, P.E. (2025). If the Null Fits, You Must Omit: Ubiquitous False Detections of Group Polarization.

A copy of the full license is available at: https://creativecommons.org/licenses/by/4.0/