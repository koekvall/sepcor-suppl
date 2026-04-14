# Data example: Dissolved oxygen concentration (Upper Mississippi River)
# Fits separable correlation model, computes standard errors via expected
# Fisher information, tests H0: C_2 = I_3 (no seasonal correlation) via a
# Wald test, and performs a parametric bootstrap LRT of H0: separable
# covariance vs H1: separable correlation (B = 10000).
#
# Usage: Rscript data_example/data_example.R

library(sepcor)
library(splines)
library(dplyr)

set.seed(2031)

# ---- Load and prepare data -------------------------------------------------
use_data_do <- readRDS("data_example/use_data_do.rds") %>%
  mutate(location = interaction(fs, strat)) %>%
  droplevels() %>%
  arrange(year, season, location)

n_r <- length(unique(use_data_do$location))  # 16 areas
n_c <- length(unique(use_data_do$season))    # 3 seasons
n   <- length(unique(use_data_do$year))
q   <- n_r * n_c                             # 48

# ---- Fit separable correlation and covariance models -----------------------
X  <- cbind(1, bs(unique(use_data_do$year), df = 5, degree = 3))
Y  <- matrix(use_data_do$DO, nrow = n, byrow = TRUE)
E0 <- t(Y - X %*% solve(crossprod(X), crossprod(X, Y)))

fit_cor <- sepcor(E0, n_r)
fit_cov <- sepcor(E0, n_r, sepcov = TRUE)
cat(sprintf("Sep. correlation converged in %d iterations (info = %d)\n",
            fit_cor$iter, fit_cor$info))
cat(sprintf("Sep. covariance  converged in %d iterations (info = %d)\n",
            fit_cov$iter, fit_cov$info))

# ---- Parametric bootstrap LRT: H0: sep. cov. vs H1: sep. cor. -------------
# Observed test statistic
T_obs <- 2 * (fit_cor$ll - fit_cov$ll)
B     <- 10000L

# Generate bootstrap samples under H0 (separable covariance)
Lc_null <- kronecker(t(chol(fit_cov$C2)), t(chol(fit_cov$C1)))
T_boot  <- numeric(B)
for (b in seq_len(B)) {
  E_b    <- Lc_null %*% matrix(rnorm(q * n), ncol = n)
  fc_b   <- tryCatch(sepcor(E_b, n_r), error = function(e) list(info = -1L))
  fv_b   <- tryCatch(sepcor(E_b, n_r, sepcov = TRUE),
                      error = function(e) list(info = -1L))
  if (fc_b$info == 0L && fv_b$info == 0L)
    T_boot[b] <- 2 * (fc_b$ll - fv_b$ll)
  else
    T_boot[b] <- NA_real_
}

p_boot <- mean(T_boot >= T_obs, na.rm = TRUE)
cat(sprintf("\nParametric bootstrap LRT: H0: sep. covariance vs H1: sep. correlation\n"))
cat(sprintf("  T_obs = %.2f, B = %d, p-value = %.4f\n", T_obs, B, p_boot))
cat(sprintf("  Bootstrap samples with convergence issues: %d / %d\n",
            sum(is.na(T_boot)), B))

# ---- Standard errors via expected Fisher information -----------------------
se <- sepcor_se(fit_cor, E0, n_r)

# ---- C_2 estimates and SEs -------------------------------------------------
# C_2 is n_c x n_c = 3x3 season correlation matrix
# upper.tri() entries in column-major order: (1,2), (1,3), (2,3)
n_U       <- n_c * (n_c - 1L) / 2L          # = 3
theta_C2  <- fit_cor$C2[upper.tri(fit_cor$C2)]      # estimated off-diagonals
Var_C2    <- se$vcov[seq_len(n_U), seq_len(n_U)]    # asymptotic variance

cat("\nC_2 (season correlation matrix):\n")
print(round(fit_cor$C2, 3))

cat("\nUpper-triangular entries of C_2:\n")
cat(sprintf("  entry  estimate     SE    z-value\n"))
for (k in seq_len(n_U)) {
  cat(sprintf("  (%d,%d)   %6.3f  %6.3f   %6.2f\n",
              which(upper.tri(fit_cor$C2), arr.ind = TRUE)[k, 1],
              which(upper.tri(fit_cor$C2), arr.ind = TRUE)[k, 2],
              theta_C2[k], se$se_C2[k], theta_C2[k] / se$se_C2[k]))
}

# ---- Joint Wald test: H0: C_2 = I_3 (all off-diagonals = 0) ---------------
W       <- as.numeric(t(theta_C2) %*% solve(Var_C2) %*% theta_C2)
p_wald  <- pchisq(W, df = n_U, lower.tail = FALSE)

cat(sprintf("\nWald test H0: C_2 = I_%d\n", n_c))
cat(sprintf("  W = %.2f  (df = %d)  p-value = %.4f\n", W, n_U, p_wald))
