# Benchmark: sepcor block-coordinate descent vs. optim (BFGS)
#
# Both optimize the same profile log-likelihood. The custom algorithm
# exploits the Kronecker structure; optim treats it as a dense problem.

library(sepcor)
library(MASS)

# ---- Negative log-likelihood for optim --------------------------------
# Parameters packed as: [C2 upper tri | C1 upper tri | log(D)]
make_nll <- function(E, nr, nc) {
  q   <- nr * nc
  n_U <- nc * (nc - 1L) / 2L
  n_V <- nr * (nr - 1L) / 2L

  function(theta) {
    C2 <- diag(nc)
    C2[upper.tri(C2)] <- theta[seq_len(n_U)]
    C2[lower.tri(C2)] <- t(C2)[lower.tri(C2)]

    C1 <- diag(nr)
    C1[upper.tri(C1)] <- theta[n_U + seq_len(n_V)]
    C1[lower.tri(C1)] <- t(C1)[lower.tri(C1)]

    D_diag <- exp(theta[n_U + n_V + seq_len(q)])

    Sigma <- diag(D_diag) %*% kronecker(C2, C1) %*% diag(D_diag)
    L <- tryCatch(t(chol(Sigma)), error = function(e) NULL)
    if (is.null(L)) return(1e20)
    -prof_log_lik(L, E)
  }
}

# Initial theta: C1 = I, C2 = I, D = sqrt(diag(S))
make_theta0 <- function(E, nr, nc) {
  n_U <- nc * (nc - 1L) / 2L
  n_V <- nr * (nr - 1L) / 2L
  q   <- nr * nc
  S   <- tcrossprod(E) / ncol(E)
  c(rep(0, n_U), rep(0, n_V), log(sqrt(diag(S))))
}

# ---- Run benchmark for one (nr, nc, n) setting -------------------------
run_one <- function(nr, nc, n, seed = 1) {
  set.seed(seed)
  C1    <- toeplitz(0.5^seq(0, nr - 1))
  C2    <- toeplitz(0.4^seq(0, nc - 1))
  D     <- diag(seq(0.5, 2, length.out = nr * nc))
  Sigma <- D %*% kronecker(C2, C1) %*% D
  E     <- t(mvrnorm(n, mu = rep(0, nr * nc), Sigma = Sigma))

  nll     <- make_nll(E, nr, nc)
  theta0  <- make_theta0(E, nr, nc)

  t_custom <- system.time(fit <- sepcor(E, nr, tol = 1e-8))["elapsed"]
  t_optim  <- system.time(
    opt <- optim(theta0, nll, method = "BFGS",
                 control = list(maxit = 5000, reltol = 1e-8))
  )["elapsed"]

  data.frame(
    nr       = nr,
    nc       = nc,
    n        = n,
    t_custom = t_custom,
    t_optim  = t_optim,
    speedup  = t_optim / t_custom,
    ll_custom = fit$ll,
    ll_optim  = -opt$value,
    converged = opt$convergence == 0
  )
}

# ---- Settings ----------------------------------------------------------
settings <- list(
  c(nr = 3,  nc = 4,  n = 100),
  c(nr = 3,  nc = 4,  n = 500),
  c(nr = 5,  nc = 6,  n = 200),
  c(nr = 5,  nc = 6,  n = 1000),
  c(nr = 8,  nc = 10, n = 500)
)

results <- do.call(rbind, lapply(settings, function(s) {
  cat("Running nr =", s["nr"], ", nc =", s["nc"], ", n =", s["n"], "... ")
  res <- run_one(s["nr"], s["nc"], s["n"])
  cat("done\n")
  res
}))

print(results, digits = 3)
