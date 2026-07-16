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
t_boot  <- system.time(
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
)["elapsed"]

p_boot <- mean(T_boot >= T_obs, na.rm = TRUE)
cat(sprintf("\nParametric bootstrap LRT: H0: sep. covariance vs H1: sep. correlation\n"))
cat(sprintf("  T_obs = %.2f, B = %d, p-value = %.4f\n", T_obs, B, p_boot))
cat(sprintf("  Bootstrap samples with convergence issues: %d / %d\n",
            sum(is.na(T_boot)), B))
cat(sprintf("  Total bootstrap time (%d model fits): %.1f seconds\n", 2L * B, t_boot))

# ---- Timing: single fit, Algorithm 1 vs BFGS with analytic gradients -------
# Backs the runtime comparison in Section 3 of the manuscript. BFGS maximizes
# the same Kronecker-structured log-likelihood (prof_log_lik_sep) with the
# analytic score and the same relative tolerance (1e-8) as Algorithm 1.
# Parameters packed as [C2 upper tri | C1 upper tri | log(D)]. No random
# numbers are drawn here, so results elsewhere in the script are unaffected.
make_nll <- function(E, nr, nc) {
  n_U <- nc * (nc - 1L) / 2L
  n_V <- nr * (nr - 1L) / 2L
  q   <- nr * nc
  function(theta) {
    C2 <- diag(nc)
    C2[upper.tri(C2)] <- theta[seq_len(n_U)]
    C2[lower.tri(C2)] <- t(C2)[lower.tri(C2)]
    C1 <- diag(nr)
    C1[upper.tri(C1)] <- theta[n_U + seq_len(n_V)]
    C1[lower.tri(C1)] <- t(C1)[lower.tri(C1)]
    D_diag <- exp(theta[n_U + n_V + seq_len(q)])
    ll <- tryCatch(prof_log_lik_sep(C1, C2, D_diag, E),
                   error = function(e) -1e20)
    if (!is.finite(ll)) return(1e20)
    -ll
  }
}

# Analytic score. With Sigma = D (C2 %x% C1) D, Sbar the average residual
# outer product, G = Sigma^-1 Sbar Sigma^-1 - Sigma^-1, and H = D G D:
#   d ll / d [C1]_{ab}  = n [sum_{j,j'} [C2]_{j'j} H_{(j,j')}]_{ab},
#   d ll / d [C2]_{ab}  = n tr(H_{(a,b)} C1) (symmetrized),
#   d ll / d log d_j    = n ([Sigma^-1 Sbar]_{jj} - 1),
# where H_{(j,j')} is the (j,j')th r x r block of H.
make_gr <- function(E, nr, nc) {
  n_U   <- nc * (nc - 1L) / 2L
  n_V   <- nr * (nr - 1L) / 2L
  q     <- nr * nc
  n_obs <- ncol(E)
  Sbar  <- tcrossprod(E) / n_obs
  function(theta) {
    C2 <- diag(nc)
    C2[upper.tri(C2)] <- theta[seq_len(n_U)]
    C2[lower.tri(C2)] <- t(C2)[lower.tri(C2)]
    C1 <- diag(nr)
    C1[upper.tri(C1)] <- theta[n_U + seq_len(n_V)]
    C1[lower.tri(C1)] <- t(C1)[lower.tri(C1)]
    d  <- exp(theta[n_U + n_V + seq_len(q)])
    Sigma <- d * t(d * kronecker(C2, C1))
    Si  <- tryCatch(chol2inv(chol(Sigma)), error = function(e) NULL)
    if (is.null(Si)) return(rep(0, length(theta)))
    SiS <- Si %*% Sbar
    G   <- SiS %*% Si - Si
    H   <- d * t(d * G)
    M1  <- matrix(0, nr, nr)
    T2  <- matrix(0, nc, nc)
    for (j in seq_len(nc)) for (jp in seq_len(nc)) {
      blk <- H[((j - 1L) * nr + 1L):(j * nr), ((jp - 1L) * nr + 1L):(jp * nr)]
      M1  <- M1 + C2[jp, j] * blk
      T2[j, jp] <- sum(blk * C1)
    }
    g_C2 <- n_obs * ((T2 + t(T2)) / 2)[upper.tri(T2)]
    g_C1 <- n_obs * M1[upper.tri(M1)]
    g_ld <- n_obs * (diag(SiS) - 1)
    -c(g_C2, g_C1, g_ld)
  }
}

nll    <- make_nll(E0, n_r, n_c)
gr     <- make_gr(E0, n_r, n_c)
S0     <- tcrossprod(E0) / n
theta0 <- c(rep(0, n_c * (n_c - 1L) / 2L), rep(0, n_r * (n_r - 1L) / 2L),
            log(sqrt(diag(S0))))

loop <- 100L  # repeat to get stable sub-millisecond timings
t_alg1 <- system.time(
  for (j in seq_len(loop)) sepcor(E0, n_r)
)["elapsed"] / loop
t_bfgs <- system.time(
  for (j in seq_len(5L)) opt <- optim(theta0, nll, gr, method = "BFGS",
                                      control = list(maxit = 5000, reltol = 1e-8))
)["elapsed"] / 5

cat(sprintf("\nSingle-fit timing (r = %d, c = %d, n = %d):\n", n_r, n_c, n))
cat(sprintf("  Algorithm 1:              %.4f seconds\n", t_alg1))
cat(sprintf("  BFGS (analytic gradient): %.4f seconds\n", t_bfgs))
cat(sprintf("  Implied time for %d BFGS fits: %.1f minutes\n",
            2L * B, 2 * B * t_bfgs / 60))
cat(sprintf("  BFGS log-lik minus Algorithm 1 log-lik: %.2e\n",
            -opt$value - fit_cor$ll))

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

# ---- Multi-start diagnostic: 150 random starting values --------------------
# Perturbs the starting values that shape the trajectory: the correlation
# matrix that is not updated first, drawn as a rescaled Wishart with one more
# degree of freedom than its dimension, and the log-standard deviations, set
# to their sample values plus N(0, 0.5^2) noise. Same mechanism as the
# multi-start option in sepcor(). All starts should reach the maximizer found
# from the default start (fit_cor).
n_starts_diag <- 150L
S0    <- tcrossprod(E0) / n
par_0 <- c(fit_cor$C1[upper.tri(fit_cor$C1)],
           fit_cor$C2[upper.tri(fit_cor$C2)],
           as.vector(fit_cor$D))

rand_corr <- function(d){
  W  <- rWishart(1, d + 1, diag(d))[,,1]
  Di <- diag(1 / sqrt(diag(W)))
  Di %*% W %*% Di
}

ll_diff  <- rep(NA_real_, n_starts_diag)
par_diff <- rep(NA_real_, n_starts_diag)
for (s in seq_len(n_starts_diag)) {
  W_init <- sqrt(diag(S0)) * exp(rnorm(q, 0, 0.5))
  fit_s  <- sepcor:::sepcor_rcpp(E0, W_init, n_r, 1e-8, 1000L, FALSE, 0,
                                 rand_corr(n_r))
  if (fit_s$info == 0L) {
    par_s       <- c(fit_s$C1[upper.tri(fit_s$C1)],
                     fit_s$C2[upper.tri(fit_s$C2)],
                     as.vector(fit_s$D))
    ll_diff[s]  <- abs(fit_s$ll - fit_cor$ll)
    par_diff[s] <- max(abs(par_s - par_0))
  }
}

cat(sprintf("\nMulti-start diagnostic (%d random starts):\n", n_starts_diag))
cat(sprintf("  Converged: %d / %d\n", sum(!is.na(ll_diff)), n_starts_diag))
cat(sprintf("  Max |log-lik - default|:   %.2e\n", max(ll_diff, na.rm = TRUE)))
cat(sprintf("  Max |parameter - default|: %.2e\n", max(par_diff, na.rm = TRUE)))
