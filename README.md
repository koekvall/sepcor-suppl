# Supplementary code for "Likelihood-Based Inference with Separable Correlation Matrices"

Code to reproduce the tables, figures, and data example in the paper. Each
script is self-contained and can be run on its own from the repository root.

## Prerequisites

```r
install.packages(c("ggplot2", "MASS", "splines", "dplyr", "foreach",
                   "doParallel", "doRNG"))
devtools::install_github("koekvall/sepcor")  # or install from a local copy
```

## Scripts by paper result

| Paper result | Script |
|--------------|--------|
| Figure 1 (data-example covariance heatmaps) | `scripts/fig_data_example.R` |
| Figure 2 (data-example variances by season)  | `scripts/fig_data_example.R` |
| Figure 3 (estimation error vs. sample size)  | `scripts/fig_error.R` |
| Figure 4 (test rejection rates)              | `scripts/fig_tests.R` |
| Table 1 (existence/uniqueness thresholds)    | `scripts/tab_boundary.R` |
| Table 2 (computing-time comparison)          | `scripts/tab_benchmark.R` |
| Coverage rates (Section 5 text)              | `scripts/tab_coverage.R` |
| Data example: estimates, standard errors, Wald test, bootstrap LRT, timings, multi-start diagnostic (Section 3) | `data_example/data_example.R` |

Figures are written to `Figures/` (git-ignored); the `tab_*` scripts print a
LaTeX table and cache their results in `scripts/*_results.rds`.

## Data

- `data_example/use_data_do.rds` — prepared dissolved oxygen data (Upper
  Mississippi River) used by the data-example scripts.

## Other scripts

- `counterexample.R` — numerical check of a conjectured sample-size condition;
  not used for any paper float.
- `scripts/benchmark.R`, `scripts/fig_benchmark.R`, `scripts/fig_boundary.R` —
  exploratory timing and threshold plots not included in the paper.

## Reproducing a single result

From the repository root, for example:

```bash
Rscript scripts/fig_data_example.R
Rscript data_example/data_example.R
```
