# Table: Computing time of sepcor vs. BFGS with analytic gradients
# Averages over 5 replications per setting.
#
# Usage: Rscript scripts/tab_benchmark.R
# Output: printed LaTeX table and scripts/tab_benchmark_results.rds

library(sepcor)
library(MASS)

set.seed(2030)
reps <- 5

# ---- Negative log-likelihood for optim (Kronecker-aware) ------------------
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

make_theta0 <- function(E, nr, nc) {
  n_U <- nc * (nc - 1L) / 2L
  n_V <- nr * (nr - 1L) / 2L
  q   <- nr * nc
  S   <- tcrossprod(E) / ncol(E)
  c(rep(0, n_U), rep(0, n_V), log(sqrt(diag(S))))
}

# ---- Analytic score for optim ----------------------------------------------
# With Sigma = D (C2 %x% C1) D, Sbar the average residual outer product,
# G = Sigma^-1 Sbar Sigma^-1 - Sigma^-1, and H = D G D:
#   d ll / d [C1]_{ab}  = n [sum_{j,j'} [C2]_{j'j} H_{(j,j')}]_{ab},
#   d ll / d [C2]_{ab}  = n tr(H_{(a,b)} C1) (symmetrized),
#   d ll / d log d_j    = n ([Sigma^-1 Sbar]_{jj} - 1),
# where H_{(j,j')} is the (j,j')th r x r block of H. Verified against
# central finite differences of prof_log_lik_sep.
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

# ---- Run one setting, averaged over reps ----------------------------------
run_one <- function(nr, nc, n) {
  C1    <- toeplitz(0.5^seq(0, nr - 1))
  C2    <- toeplitz(0.4^seq(0, nc - 1))
  D     <- diag(seq(0.5, 2, length.out = nr * nc))
  Sigma <- D %*% kronecker(C2, C1) %*% D

  t_custom <- t_optim <- ll_diff <- est_diff <- numeric(reps)
  loop     <- 20L   # repeat sepcor to get sub-ms precision
  n_U <- nc * (nc - 1L) / 2L
  n_V <- nr * (nr - 1L) / 2L
  for (k in seq_len(reps)) {
    E      <- t(mvrnorm(n, mu = rep(0, nr * nc), Sigma = Sigma))
    theta0 <- make_theta0(E, nr, nc)
    # Time sepcor via repeated calls; report in milliseconds
    t_custom[k] <- system.time(
      for (j in seq_len(loop)) fit <- sepcor(E, nr, tol = 1e-8)
    )["elapsed"] / loop * 1000
    # BFGS with the analytic score and the same relative tolerance
    t_optim[k]  <- system.time(
      opt <- optim(theta0, make_nll(E, nr, nc), make_gr(E, nr, nc),
                   method = "BFGS",
                   control = list(maxit = 5000, reltol = 1e-8))
    )["elapsed"] * 1000
    # Agreement between the two solutions
    C2_b <- diag(nc); C2_b[upper.tri(C2_b)] <- opt$par[seq_len(n_U)]
    C2_b[lower.tri(C2_b)] <- t(C2_b)[lower.tri(C2_b)]
    C1_b <- diag(nr); C1_b[upper.tri(C1_b)] <- opt$par[n_U + seq_len(n_V)]
    C1_b[lower.tri(C1_b)] <- t(C1_b)[lower.tri(C1_b)]
    d_b  <- exp(opt$par[n_U + n_V + seq_len(nr * nc)])
    ll_diff[k]  <- fit$ll - (-opt$value)
    est_diff[k] <- max(abs(c(C1_b[upper.tri(C1_b)] - fit$C1[upper.tri(fit$C1)],
                             C2_b[upper.tri(C2_b)] - fit$C2[upper.tri(fit$C2)],
                             d_b - as.vector(fit$D))))
  }

  data.frame(nr = nr, nc = nc, n = n,
             t_sepcor_ms  = mean(t_custom),
             t_bfgs_ms    = mean(t_optim),
             speedup      = mean(t_optim) / mean(t_custom),
             ll_diff_min  = min(ll_diff),
             est_diff_max = max(est_diff))
}

# ---- Settings -------------------------------------------------------------
settings <- list(
  c(nr = 5,  nc = 6,  n = 30),
  c(nr = 8,  nc = 8,  n = 50),
  c(nr = 6,  nc = 10, n = 40),
  c(nr = 10, nc = 10, n = 50),
  c(nr = 8,  nc = 16, n = 50),
  c(nr = 12, nc = 12, n = 50),
  c(nr = 16, nc = 16, n = 50)
)

results <- do.call(rbind, lapply(settings, function(s) {
  cat(sprintf("r=%d, c=%d, n=%d ... ", s["nr"], s["nc"], s["n"]))
  res <- run_one(s["nr"], s["nc"], s["n"])
  cat(sprintf("sepcor=%.2fms  BFGS=%.0fms  speedup=%.0fx  ll_diff_min=%.1e  est_diff_max=%.1e\n",
              res$t_sepcor_ms, res$t_bfgs_ms, res$speedup,
              res$ll_diff_min, res$est_diff_max))
  res
}))

# ---- Print LaTeX table ----------------------------------------------------
cat("\n=== LaTeX table ===\n")
cat("\\begin{tabular}{rr r rr r}\n")
cat("$r$ & $c$ & $n$ & sepcor (ms) & BFGS (ms) & Speedup \\\\\n")
cat("\\hline\n")
for (i in seq_len(nrow(results))) {
  cat(sprintf("%d & %d & %4d & %5.2f & %5.0f & %4.0f$\\times$ \\\\\n",
              results$nr[i], results$nc[i], results$n[i],
              results$t_sepcor_ms[i], results$t_bfgs_ms[i], results$speedup[i]))
}
cat("\\end{tabular}\n")

saveRDS(results, "scripts/tab_benchmark_results.rds")
cat("Saved to scripts/tab_benchmark_results.rds\n")
