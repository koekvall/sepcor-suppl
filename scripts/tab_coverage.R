# Table: Coverage of Wald confidence intervals based on expected Fisher information
# DGP: separable correlation with r = c = 5, using AR(1) component matrices.
# For each off-diagonal entry of C_1 and C_2, build 95% CI as estimate +/- 1.96*SE
# where SE comes from sepcor_se() (expected Fisher information).
# Report average coverage across off-diagonal entries of C_1 and C_2 separately.
#
# Usage: Rscript scripts/tab_coverage.R
# Output: printed table and scripts/tab_coverage_results.rds

library(sepcor)
library(parallel)

set.seed(2029)
m       <- 500     # replications per setting
n_cores <- max(1L, detectCores() - 1L)
alpha   <- 0.05
z       <- qnorm(1 - alpha / 2)

# ---- DGP: AR(1) component matrices with rho = 0.6 --------------------------
r   <- 5L
cc  <- 5L
q   <- r * cc
rho <- 0.6

# AR(1) correlation matrix
ar1 <- function(d, rho) rho^abs(outer(seq_len(d), seq_len(d), "-"))
C1_true <- ar1(r,  rho)
C2_true <- ar1(cc, rho)

# True parameter values (off-diagonals of C1 and C2)
theta_C2_true <- C2_true[upper.tri(C2_true)]   # length r*(r-1)/2 = 10
theta_C1_true <- C1_true[upper.tri(C1_true)]   # length cc*(cc-1)/2 = 10

# Cholesky factor for data generation (D = I_q for sepcor DGP)
Sig_c <- kronecker(t(chol(C2_true)), t(chol(C1_true)))

# ---- Simulation -------------------------------------------------------------
n_vals <- c(50L, 100L, 320L)

results <- data.frame()

for (n_sim in n_vals) {
  res <- mclapply(seq_len(m), function(j) {
    E  <- Sig_c %*% matrix(rnorm(q * n_sim), ncol = n_sim)
    fc <- tryCatch(sepcor(E, r), error = function(e) list(info = -1L))
    if (fc$info != 0L) return(NULL)

    se_obj <- tryCatch(sepcor_se(fc, E, r),
                       error = function(e) NULL)
    if (is.null(se_obj)) return(NULL)

    # Off-diagonal entries of C_2 (upper tri, column-major)
    est_C2 <- fc$C2[upper.tri(fc$C2)]
    se_C2  <- se_obj$se_C2
    cov_C2 <- as.integer(
      abs(est_C2 - theta_C2_true) <= z * se_C2
    )

    # Off-diagonal entries of C_1 (upper tri, column-major)
    est_C1 <- fc$C1[upper.tri(fc$C1)]
    se_C1  <- se_obj$se_C1
    cov_C1 <- as.integer(
      abs(est_C1 - theta_C1_true) <= z * se_C1
    )

    list(cov_C2 = cov_C2, cov_C1 = cov_C1)
  }, mc.cores = n_cores)

  # Drop NULLs (non-converged or SE failures)
  res <- Filter(Negate(is.null), res)
  n_valid <- length(res)

  # Average coverage per entry, then across entries
  cov_C2_mat <- do.call(rbind, lapply(res, `[[`, "cov_C2"))
  cov_C1_mat <- do.call(rbind, lapply(res, `[[`, "cov_C1"))

  avg_cov_C2 <- round(mean(colMeans(cov_C2_mat)), 3)
  avg_cov_C1 <- round(mean(colMeans(cov_C1_mat)), 3)

  results <- rbind(results, data.frame(
    n       = n_sim,
    n_valid = n_valid,
    cov_C1  = avg_cov_C1,
    cov_C2  = avg_cov_C2
  ))
  cat(sprintf("n = %3d  valid = %d  cov(C1) = %.3f  cov(C2) = %.3f\n",
              n_sim, n_valid, avg_cov_C1, avg_cov_C2))
}

# ---- Print -----------------------------------------------------------------
cat("\n=== Coverage results ===\n")
print(results)

saveRDS(results, "scripts/tab_coverage_results.rds")
cat("Saved to scripts/tab_coverage_results.rds\n")
