# Figure: Estimation error vs sample size
# Compares separable correlation, separable covariance, and unrestricted MLE
# DGPs are the fitted separable correlation and separable covariance models
# from the dissolved oxygen data example (r = 16 areas, c = 3 seasons).
#
# Usage: Rscript scripts/fig_error.R
# Output: Figures/fig_error.pdf

library(sepcor)
library(splines)
library(dplyr)
library(ggplot2)
library(patchwork)

set.seed(2026)

# ---- Fit models to dissolved oxygen data to obtain DGP parameters ----------
use_data_do <- readRDS("~/GDrive/Work/Consulting/USGS/Data/use_data_do.rds") %>%
  mutate(location = interaction(fs, strat)) %>%
  droplevels() %>%
  arrange(year, season, location)

n_r <- length(unique(use_data_do$location))  # 16 areas
n_c <- length(unique(use_data_do$season))    # 3 seasons
n   <- length(unique(use_data_do$year))
q   <- n_r * n_c

X  <- cbind(1, bs(unique(use_data_do$year), df = 5, degree = 3))
Y  <- matrix(use_data_do$DO, nrow = n, byrow = TRUE)
E0 <- t(Y - X %*% solve(crossprod(X), crossprod(X, Y)))

fit_cor <- sepcor(E0, n_r)
fit_cov <- sepcor(E0, n_r, sepcov = TRUE)

# Cholesky factors for data generation
Sig_c_cor <- diag(as.vector(fit_cor$D)) %*%
  kronecker(t(chol(fit_cor$C2)), t(chol(fit_cor$C1)))
Sig_c_cov <- kronecker(t(chol(fit_cov$C2)), t(chol(fit_cov$C1)))

Sig0_cor <- tcrossprod(Sig_c_cor)
Sig0_cov <- tcrossprod(Sig_c_cov)

# ---- Simulation ------------------------------------------------------------
m      <- 200  # replications per setting
n_vals <- c(15, 20, 25, 35, 50, 75, 100, 150, 200, 320, 500)

results <- data.frame()

for (dgp in c("sepcor", "sepcov")) {
  Sig_c <- if (dgp == "sepcor") Sig_c_cor else Sig_c_cov
  Sig0  <- if (dgp == "sepcor") Sig0_cor  else Sig0_cov

  for (n_sim in n_vals) {
    err_cor <- err_cov <- err_ur <- numeric(m)

    for (j in seq_len(m)) {
      E <- Sig_c %*% matrix(rnorm(q * n_sim), ncol = n_sim)

      # Separable correlation
      fc <- tryCatch(sepcor(E, n_r), error = function(e) list(info = -1))
      if (fc$info == 0) {
        Sh <- diag(as.vector(fc$D)) %*% kronecker(fc$C2, fc$C1) %*%
          diag(as.vector(fc$D))
        err_cor[j] <- norm(Sh - Sig0, type = "2")
      } else {
        err_cor[j] <- NA
      }

      # Separable covariance
      fv <- tryCatch(sepcor(E, n_r, sepcov = TRUE),
                     error = function(e) list(info = -1))
      if (fv$info == 0) {
        Sh <- kronecker(fv$C2, fv$C1)
        err_cov[j] <- norm(Sh - Sig0, type = "2")
      } else {
        err_cov[j] <- NA
      }

      # Unrestricted (sample covariance, only when n > q)
      if (n_sim > q) {
        S <- tcrossprod(E) / n_sim
        err_ur[j] <- norm(S - Sig0, type = "2")
      } else {
        err_ur[j] <- NA
      }
    }

    results <- rbind(results,
      data.frame(dgp = dgp, n = n_sim,
                 method = "Sep. correlation", err = mean(err_cor, na.rm = TRUE)),
      data.frame(dgp = dgp, n = n_sim,
                 method = "Sep. covariance",  err = mean(err_cov, na.rm = TRUE)),
      data.frame(dgp = dgp, n = n_sim,
                 method = "Unrestricted",     err = mean(err_ur,  na.rm = TRUE))
    )
    cat(sprintf("Done: dgp = %s, n = %d\n", dgp, n_sim))
  }
}

# ---- Plot ------------------------------------------------------------------
col_vals <- c("Sep. correlation" = "#2166ac",
              "Sep. covariance"  = "#d6604d",
              "Unrestricted"     = "#bababa")

common_theme <- list(
  theme_bw(base_size = 10),
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        plot.title       = element_text(size = 10))
)

make_panel <- function(data, title) {
  ggplot(data, aes(x = n, y = err, color = method, shape = method)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2) +
    scale_x_log10() +
    scale_color_manual(values = col_vals) +
    labs(x = "Sample size (n)", y = "Average spectral norm error",
         color = "Estimator", shape = "Estimator", title = title) +
    common_theme
}

p_a <- make_panel(results[results$dgp == "sepcor", ],
                  "(a) DGP: separable correlation")
p_b <- make_panel(results[results$dgp == "sepcov", ],
                  "(b) DGP: separable covariance")

p <- (p_a | p_b) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("Figures/fig_error.pdf", p, width = 10, height = 3.5)
cat("Saved Figures/fig_error.pdf\n")
