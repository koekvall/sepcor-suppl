# Table: Computing time of sepcor vs. L-BFGS-B
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

# ---- Run one setting, averaged over reps ----------------------------------
run_one <- function(nr, nc, n) {
  C1    <- toeplitz(0.5^seq(0, nr - 1))
  C2    <- toeplitz(0.4^seq(0, nc - 1))
  D     <- diag(seq(0.5, 2, length.out = nr * nc))
  Sigma <- D %*% kronecker(C2, C1) %*% D

  t_custom <- t_optim <- numeric(reps)
  loop     <- 20L   # repeat sepcor to get sub-ms precision
  for (k in seq_len(reps)) {
    E      <- t(mvrnorm(n, mu = rep(0, nr * nc), Sigma = Sigma))
    theta0 <- make_theta0(E, nr, nc)
    # Time sepcor via repeated calls; report in milliseconds
    t_custom[k] <- system.time(
      for (j in seq_len(loop)) sepcor(E, nr, tol = 1e-8)
    )["elapsed"] / loop * 1000
    t_optim[k]  <- system.time(
      optim(theta0, make_nll(E, nr, nc), method = "BFGS",
            control = list(maxit = 5000, reltol = 1e-8))
    )["elapsed"]
  }

  data.frame(nr = nr, nc = nc, n = n,
             t_sepcor_ms = mean(t_custom),
             t_bfgs_s    = mean(t_optim),
             speedup     = mean(t_optim) / (mean(t_custom) / 1000))
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
  cat(sprintf("sepcor=%.2fms  BFGS=%.3fs  speedup=%.1fx\n",
              res$t_sepcor_ms, res$t_bfgs_s, res$speedup))
  res
}))

# ---- Print LaTeX table ----------------------------------------------------
cat("\n=== LaTeX table ===\n")
cat("\\begin{tabular}{rr r rr r}\n")
cat("$r$ & $c$ & $n$ & sepcor (ms) & BFGS (s) & Speedup \\\\\n")
cat("\\hline\n")
for (i in seq_len(nrow(results))) {
  cat(sprintf("%d & %d & %4d & %5.2f & %5.2f & %4.0f$\\times$ \\\\\n",
              results$nr[i], results$nc[i], results$n[i],
              results$t_sepcor_ms[i], results$t_bfgs_s[i], results$speedup[i]))
}
cat("\\end{tabular}\n")

saveRDS(results, "scripts/tab_benchmark_results.rds")
cat("Saved to scripts/tab_benchmark_results.rds\n")
