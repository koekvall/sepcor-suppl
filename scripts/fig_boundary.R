# Figure: Convergence and uniqueness boundary
# Fraction of datasets where algorithm converges/yields unique solution, vs n
# Vertical dotted lines at separable covariance existence threshold
#
# Usage: Rscript scripts/fig_boundary.R
# Output: Figures/fig_boundary.pdf

library(sepcor)
library(ggplot2)

set.seed(2028)

configs <- list(
  list(r = 2, c = 9),
  list(r = 3, c = 5),
  list(r = 5, c = 5),
  list(r = 5, c = 15),
  list(r = 3, c = 3)
)

n_datasets <- 100
n_starts <- 10

results <- data.frame()

for(cfg in configs){
  r <- cfg$r; cc <- cfg$c; q <- r * cc
  U0 <- 0.5^abs(outer(1:cc, 1:cc, "-"))
  V0 <- 0.5^abs(outer(1:r, 1:r, "-"))
  Sig_c <- kronecker(t(chol(U0)), t(chol(V0)))

  sepcov_bound <- max(r / cc, cc / r) + 1

  n_vals <- sort(unique(c(2:12, 15, 20, 30)))

  for(n in n_vals){
    n_converged <- 0
    n_unique <- 0
    for(d in 1:n_datasets){
      E <- Sig_c %*% matrix(rnorm(q * n), ncol = n)
      S <- tcrossprod(E) / n
      lls <- numeric(n_starts)
      all_conv <- TRUE
      for(s in 1:n_starts){
        if(s == 1){
          W_init <- diag(S)
        } else {
          W_init <- (sqrt(pmax(diag(S), 1e-6)) * exp(rnorm(q, 0, 0.5)))^2
        }
        fit <- suppressWarnings(
          sepcor:::sepcor_rcpp(E, W_init, r, 1e-12, 3000, FALSE, 0)
        )
        lls[s] <- fit$ll
        if(fit$info != 0) all_conv <- FALSE
      }
      if(all_conv){
        n_converged <- n_converged + 1
        if(max(lls) - min(lls) < 1e-3) n_unique <- n_unique + 1
      }
    }
    results <- rbind(results, data.frame(
      r = r, c = cc, n = n,
      frac_conv = n_converged / n_datasets,
      frac_unique = if(n_converged > 0) n_unique / n_converged else NA,
      sepcov_bound = sepcov_bound
    ))
    cat(sprintf("r=%d, c=%d, n=%d: conv=%.2f, unique=%s\n",
      r, cc, n, n_converged / n_datasets,
      if(n_converged > 0) sprintf("%.2f", n_unique / n_converged) else "NA"))
  }
}

results$facet <- paste0("r = ", results$r, ", c = ", results$c)

# Threshold lines per facet
thresholds <- unique(results[, c("facet", "sepcov_bound")])

p <- ggplot(results, aes(x = n)) +
  geom_line(aes(y = frac_conv, color = "Convergence"), linewidth = 0.7) +
  geom_point(aes(y = frac_conv, color = "Convergence"), size = 1.8) +
  geom_line(aes(y = frac_unique, color = "Uniqueness"), linewidth = 0.7,
            linetype = "dashed", na.rm = TRUE) +
  geom_point(aes(y = frac_unique, color = "Uniqueness"), size = 1.8,
             shape = 2, na.rm = TRUE) +
  geom_vline(data = thresholds, aes(xintercept = sepcov_bound),
             linetype = "dotted", color = "grey40", linewidth = 0.5) +
  facet_wrap(~ facet, nrow = 1) +
  scale_color_manual(values = c("Convergence" = "#2166AC",
                                "Uniqueness" = "#B2182B")) +
  ylim(0, 1.05) +
  labs(x = "Sample size (n)", y = "Fraction",
       color = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 9),
        panel.grid.minor = element_blank())

ggsave("Figures/fig_boundary.pdf", p, width = 8, height = 3)
cat("Saved Figures/fig_boundary.pdf\n")
