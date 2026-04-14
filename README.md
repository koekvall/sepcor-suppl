# Supplementary code for "Likelihood-Based Inference with Separable Correlation Matrices"

Code to reproduce the simulations, tables, figures, and data example in the paper.

## Prerequisites

```r
install.packages(c("ggplot2", "tidyr", "patchwork", "splines", "dplyr"))
devtools::install_github("koekvall/sepcor")  # or install from local copy
```

## Structure

- `scripts/run_all.R` — master script that sources the simulation scripts below
- `scripts/fig_error.R` — Figure 1 (estimation error vs. sample size)
- `scripts/fig_tests.R` — Figure 2 (bootstrap and asymptotic test rejection rates)
- `scripts/tab_boundary.R` — Table 1 (sample size thresholds for existence/uniqueness)
- `scripts/tab_benchmark.R` — Table 2 (computing time comparison)
- `scripts/benchmark.R` — timing benchmarks
- `scripts/fig_benchmark.R` — benchmark figures
- `scripts/fig_boundary.R` — boundary threshold figures
- `scripts/tab_coverage.R` — Wald interval coverage rates
- `data_example/data_example.R` — dissolved oxygen data example (Section 3)
- `data_example/use_data_do.rds` — prepared dissolved oxygen data
- `counterexample.R` — numerical exploration of sample size conditions

## Reproducing the simulations

From the repository root:

```bash
Rscript scripts/run_all.R
```

For the data example:

```bash
Rscript data_example/data_example.R
```
