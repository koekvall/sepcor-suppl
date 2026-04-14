# Figure: Computation time of sepcor vs. L-BFGS-B (optim)
#
# Usage: Rscript scripts/fig_benchmark.R
# Output: Figures/fig_benchmark.pdf

library(sepcor)
library(MASS)
library(ggplot2)

# ---- Negative log-likelihood for optim --------------------------------
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

make_theta0 <- function(E, nr, nc) {
  n_U <- nc * (nc - 1L) / 2L
  n_V <- nr * (nr - 1L) / 2L
  q   <- nr * nc
  S   <- tcrossprod(E) / ncol(E)
  c(rep(0, n_U), rep(0, n_V), log(sqrt(diag(S))))
}

# ---- Run benchmark for one setting, averaged over reps ----------------
run_one <- function(nr, nc, n, reps = 5, seed = 1) {
  set.seed(seed)
  C1    <- toeplitz(0.5^seq(0, nr - 1))
  C2    <- toeplitz(0.4^seq(0, nc - 1))
  D     <- diag(seq(0.5, 2, length.out = nr * nc))
  Sigma <- D %*% kronecker(C2, C1) %*% D

  t_custom <- t_optim <- numeric(reps)
  for (k in seq_len(reps)) {
    E      <- t(mvrnorm(n, mu = rep(0, nr * nc), Sigma = Sigma))
    nll    <- make_nll(E, nr, nc)
    theta0 <- make_theta0(E, nr, nc)
    t_custom[k] <- system.time(sepcor(E, nr, tol = 1e-8))["elapsed"]
    t_optim[k]  <- system.time(
      optim(theta0, nll, method = "BFGS",
            control = list(maxit = 5000, reltol = 1e-8))
    )["elapsed"]
  }

  label <- sprintf("r=%d, c=%d\nn=%d", nr, nc, n)
  rbind(
    data.frame(setting = label, method = "sepcor", time = mean(t_custom)),
    data.frame(setting = label, method = "BFGS",   time = mean(t_optim))
  )
}

# ---- Settings ----------------------------------------------------------
settings <- list(
  c(nr = 3, nc = 4,  n = 100),
  c(nr = 3, nc = 4,  n = 500),
  c(nr = 5, nc = 6,  n = 200),
  c(nr = 5, nc = 6,  n = 1000),
  c(nr = 8, nc = 10, n = 500)
)

results <- do.call(rbind, lapply(settings, function(s) {
  cat(sprintf("Running r=%d, c=%d, n=%d ...\n", s["nr"], s["nc"], s["n"]))
  run_one(s["nr"], s["nc"], s["n"])
}))

# Fix factor order so settings appear left to right
results$setting <- factor(results$setting, levels = unique(results$setting))
results$method  <- factor(results$method,  levels = c("sepcor", "BFGS"))

# ---- Plot --------------------------------------------------------------
p <- ggplot(results, aes(x = setting, y = time, color = method, shape = method)) +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  scale_y_log10(
    breaks = 10^(-3:2),
    labels = c("0.001", "0.01", "0.1", "1", "10", "100")
  ) +
  scale_color_manual(values = c("sepcor" = "#2166ac", "BFGS" = "#d6604d")) +
  scale_shape_manual(values = c("sepcor" = 16, "BFGS" = 17)) +
  labs(x = NULL, y = "Time (seconds, log scale)", color = "Method",
       shape = "Method") +
  theme_bw(base_size = 10) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

ggsave("Figures/fig_benchmark.pdf", p, width = 5, height = 3.5)
cat("Saved Figures/fig_benchmark.pdf\n")
