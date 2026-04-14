# Figure: Test rejection rates
# Test (a): H0: sep. covariance vs H1: sep. correlation
#   -- bootstrap and asymptotic LRT, size only (DGP: sep. covariance)
# Test (b): H0: sep. correlation vs H1: unrestricted
#   -- bootstrap and asymptotic LRT, size only (DGP: sep. correlation)
#   -- asymptotic df = 1005 for r=16, c=3
# DGPs from dissolved oxygen data example (r=16, c=3)
#
# Usage: Rscript scripts/fig_tests.R
# Output: Figures/fig_tests.pdf

library(sepcor)
library(splines)
library(dplyr)
library(ggplot2)
library(parallel)
library(patchwork)

set.seed(2027)
m      <- 200   # replications per setting
B      <- 199   # bootstrap samples per replication
n_cores <- max(1L, detectCores() - 1L)

# ---- Fit models to data to obtain DGP parameters ---------------------------
use_data_do <- readRDS("~/GDrive/Work/Consulting/USGS/Data/use_data_do.rds") %>%
  mutate(location = interaction(fs, strat)) %>%
  droplevels() %>%
  arrange(year, season, location)

n_r <- length(unique(use_data_do$location))  # 16
n_c <- length(unique(use_data_do$season))    # 3
n   <- length(unique(use_data_do$year))
q   <- n_r * n_c                             # 48

X  <- cbind(1, bs(unique(use_data_do$year), df = 5, degree = 3))
Y  <- matrix(use_data_do$DO, nrow = n, byrow = TRUE)
E0 <- t(Y - X %*% solve(crossprod(X), crossprod(X, Y)))

fit_cor <- sepcor(E0, n_r)
fit_cov <- sepcor(E0, n_r, sepcov = TRUE)

Sig_c_cor <- diag(as.vector(fit_cor$D)) %*%
  kronecker(t(chol(fit_cor$C2)), t(chol(fit_cor$C1)))
Sig_c_cov <- kronecker(t(chol(fit_cov$C2)), t(chol(fit_cov$C1)))

# Degrees of freedom for asymptotic tests
df_a <- (n_r - 1L) * (n_c - 1L)                          # = 30 for r=16, c=3
df_b <- q*(q+1L)/2L - (n_r*(n_r-1L)/2L + n_c*(n_c-1L)/2L + q)  # = 1005

# ---- One-dataset test functions --------------------------------------------

# Test (a): returns bootstrap and asymptotic reject indicators
run_test_a <- function(E) {
  n_obs <- ncol(E)
  fc <- tryCatch(sepcor(E, n_r), error = function(e) list(info = -1L))
  fv <- tryCatch(sepcor(E, n_r, sepcov = TRUE),
                 error = function(e) list(info = -1L))
  if (fc$info != 0L || fv$info != 0L)
    return(c(boot = NA_real_, asymp = NA_real_))

  T_obs    <- 2 * (fc$ll - fv$ll)
  p_asymp  <- pchisq(T_obs, df = df_a, lower.tail = FALSE)

  # Bootstrap from sepcov null fit
  Lc_null <- kronecker(t(chol(fv$C2)), t(chol(fv$C1)))
  T_boot  <- vapply(seq_len(B), function(b) {
    E_b  <- Lc_null %*% matrix(rnorm(q * n_obs), ncol = n_obs)
    fc_b <- tryCatch(sepcor(E_b, n_r), error = function(e) list(info = -1L))
    fv_b <- tryCatch(sepcor(E_b, n_r, sepcov = TRUE),
                     error = function(e) list(info = -1L))
    if (fc_b$info == 0L && fv_b$info == 0L)
      2 * (fc_b$ll - fv_b$ll)
    else
      NA_real_
  }, numeric(1L))

  p_boot <- mean(T_boot >= T_obs, na.rm = TRUE)
  c(boot = p_boot < 0.05, asymp = p_asymp < 0.05)
}

# Test (b): bootstrap and asymptotic LRT, returns both reject indicators
run_test_b <- function(E) {
  n_obs <- ncol(E)
  fc <- tryCatch(sepcor(E, n_r), error = function(e) list(info = -1L))
  if (fc$info != 0L) return(c(boot = NA_real_, asymp = NA_real_))

  S    <- tcrossprod(E) / n_obs
  L_ur <- tryCatch(t(chol(S)), error = function(e) NULL)
  if (is.null(L_ur)) return(c(boot = NA_real_, asymp = NA_real_))

  ll_ur <- prof_log_lik(L_ur, E)
  T_obs <- 2 * (ll_ur - fc$ll)
  p_asymp <- pchisq(T_obs, df = df_b, lower.tail = FALSE)

  # Bootstrap from sepcor null fit
  Lc_null <- diag(as.vector(fc$D)) %*%
    kronecker(t(chol(fc$C2)), t(chol(fc$C1)))
  T_boot <- vapply(seq_len(B), function(b) {
    E_b  <- Lc_null %*% matrix(rnorm(q * n_obs), ncol = n_obs)
    fc_b <- tryCatch(sepcor(E_b, n_r), error = function(e) list(info = -1L))
    if (fc_b$info != 0L) return(NA_real_)
    S_b  <- tcrossprod(E_b) / n_obs
    L_b  <- tryCatch(t(chol(S_b)), error = function(e) NULL)
    if (is.null(L_b)) return(NA_real_)
    2 * (prof_log_lik(L_b, E_b) - fc_b$ll)
  }, numeric(1L))

  c(boot = mean(T_boot >= T_obs, na.rm = TRUE) < 0.05, asymp = p_asymp < 0.05)
}

# ---- Run simulations -------------------------------------------------------
cache_file <- "scripts/fig_tests_cache.rds"

if (file.exists(cache_file)) {
  cache    <- readRDS(cache_file)
  results_a <- cache$results_a
  results_b <- cache$results_b
  cat("Loaded results from cache.\n")
} else {

  # Test (a)
  n_vals_a <- c(15, 20, 25, 35, 50, 75, 100)
  results_a <- data.frame()

  for (dgp in c("sepcov", "sepcor")) {
    Sig_c <- if (dgp == "sepcov") Sig_c_cov else Sig_c_cor
    for (n_sim in n_vals_a) {
      res <- mclapply(seq_len(m), function(j) {
        E <- Sig_c %*% matrix(rnorm(q * n_sim), ncol = n_sim)
        run_test_a(E)
      }, mc.cores = n_cores)
      rej_boot  <- mean(sapply(res, `[`, "boot"),  na.rm = TRUE)
      rej_asymp <- mean(sapply(res, `[`, "asymp"), na.rm = TRUE)
      results_a <- rbind(results_a,
        data.frame(dgp = dgp, n = n_sim, test = "Bootstrap LRT",
                   rej = rej_boot),
        data.frame(dgp = dgp, n = n_sim, test = "Asymptotic LRT",
                   rej = rej_asymp))
      cat(sprintf("(a) dgp=%-7s n=%3d  boot=%.2f  asymp=%.2f\n",
                  dgp, n_sim, rej_boot, rej_asymp))
    }
  }

  # Test (b): sepcor DGP (size) only
  n_vals_b <- c(50, 60, 70, 85, 100)
  results_b <- data.frame()

  for (n_sim in n_vals_b) {
    res <- mclapply(seq_len(m), function(j) {
      E <- Sig_c_cor %*% matrix(rnorm(q * n_sim), ncol = n_sim)
      run_test_b(E)
    }, mc.cores = n_cores)
    rej_boot  <- mean(sapply(res, `[`, "boot"),  na.rm = TRUE)
    rej_asymp <- mean(sapply(res, `[`, "asymp"), na.rm = TRUE)
    results_b <- rbind(results_b,
      data.frame(n = n_sim, test = "Bootstrap LRT",  rej = rej_boot),
      data.frame(n = n_sim, test = "Asymptotic LRT", rej = rej_asymp))
    cat(sprintf("(b) n=%3d  boot=%.2f  asymp=%.2f\n", n_sim, rej_boot, rej_asymp))
  }

  saveRDS(list(results_a = results_a, results_b = results_b), cache_file)
  cat("Saved results to cache.\n")
}

# ---- Plot ------------------------------------------------------------------
# Keep only the sepcov DGP (size); power was 1.0 for all n
results_a_size <- results_a[results_a$dgp == "sepcov", ]

common_theme <- list(
  theme_bw(base_size = 10),
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        plot.title       = element_text(size = 10))
)

col_vals  <- c("Asymptotic LRT" = "#d6604d", "Bootstrap LRT" = "#2166ac")
shp_vals  <- c("Asymptotic LRT" = 16,       "Bootstrap LRT" = 17)

p_a <- ggplot(results_a_size,
    aes(x = n, y = rej, color = test, shape = test, group = test)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.05, linetype = "dotted",
             color = "grey40", linewidth = 0.4) +
  scale_x_log10() +
  scale_color_manual(values = col_vals) +
  scale_shape_manual(values = shp_vals) +
  ylim(0, 1) +
  labs(x = "Sample size (n)", y = "Rejection rate",
       color = "Test", shape = "Test",
       title = "(a) H0: sep. covariance vs H1: sep. correlation") +
  common_theme

p_b <- ggplot(results_b,
    aes(x = n, y = rej, color = test, shape = test, group = test)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.05, linetype = "dotted",
             color = "grey40", linewidth = 0.4) +
  scale_x_log10() +
  scale_color_manual(values = col_vals) +
  scale_shape_manual(values = shp_vals) +
  labs(x = "Sample size (n)", y = "Rejection rate",
       color = "Test", shape = "Test",
       title = "(b) H0: sep. correlation vs H1: unrestricted") +
  common_theme

p <- (p_a | p_b) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("Figures/fig_tests.pdf", p, width = 10, height = 3.5)
cat("Saved Figures/fig_tests.pdf\n")
