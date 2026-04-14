# Master script: reproduce all simulation figures
#
# Prerequisites:
#   install.packages(c("ggplot2", "tidyr", "patchwork"))
#   devtools::install("~/GitHub/sepcor")  # or install from GitHub
#
# Usage: from the project root directory, run:
#   Rscript scripts/run_all.R
#
# Output: PDF figures in Figures/

cat("=== Generating Figure 1: estimation error ===\n")
source("scripts/fig_error.R")

cat("\n=== Generating Figure 2: test rejection rates ===\n")
source("scripts/fig_tests.R")

cat("\n=== Generating Table: convergence boundary thresholds ===\n")
source("scripts/tab_boundary.R")

cat("\nAll figures and tables generated.\n")
22