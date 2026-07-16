# Figures 1 and 2: dissolved oxygen data example.
#   Figures/cov_plot_1.pdf : heatmap of the residual (sample) covariance matrix
#   Figures/cov_plot_2.pdf : heatmap of the fitted separable-correlation covariance
#   Figures/var_plot.pdf   : estimated variances by location and season,
#                            separable correlation (dots) vs covariance (triangles)
# Opacity encodes covariance magnitude, matching the manuscript figures. PDF
# output preserves the opacity (transparency) without requiring cairo/X11.
#
# Usage: Rscript scripts/fig_data_example.R   (from the repository root)

library(sepcor)
library(splines)
library(dplyr)
library(ggplot2)

# ---- Data and fits (same pipeline as data_example/data_example.R) ----------
use_data_do <- readRDS("data_example/use_data_do.rds") %>%
  mutate(location = interaction(fs, strat)) %>%
  droplevels() %>%
  arrange(year, season, location)

n_r <- length(unique(use_data_do$location))   # 16 locations
n_c <- length(unique(use_data_do$season))     # 3 seasons
n   <- length(unique(use_data_do$year))
q   <- n_r * n_c

X  <- cbind(1, bs(unique(use_data_do$year), df = 5, degree = 3))
Y  <- matrix(use_data_do$DO, nrow = n, byrow = TRUE)
E0 <- t(Y - X %*% solve(crossprod(X), crossprod(X, Y)))

fit_cor <- sepcor(E0, n_r)                 # separable correlation
fit_cov <- sepcor(E0, n_r, sepcov = TRUE)  # separable covariance

dir.create("Figures", showWarnings = FALSE)

# ---- Figure 1: covariance heatmaps -----------------------------------------
S_resid <- tcrossprod(E0) / n                                  # residual sample cov
d_cor   <- as.vector(fit_cor$D)
Sig_cor <- diag(d_cor) %*% kronecker(fit_cor$C2, fit_cor$C1) %*% diag(d_cor)

amax <- max(abs(S_resid), abs(Sig_cor))                        # shared opacity scale
heat_df <- function(M) data.frame(row = c(row(M)), col = c(col(M)), val = c(M))
heat_plot <- function(M) {
  ggplot(heat_df(M), aes(col, row, alpha = abs(val))) +
    geom_tile(fill = "black") +
    scale_alpha_continuous(range = c(0, 1), limits = c(0, amax), guide = "none") +
    scale_y_reverse() +
    coord_equal() +
    theme_void()
}

ggsave("Figures/cov_plot_1.pdf", heat_plot(S_resid),
       width = 7, height = 7, units = "cm")
ggsave("Figures/cov_plot_2.pdf", heat_plot(Sig_cor),
       width = 7, height = 7, units = "cm")

# ---- Figure 2: estimated variances by location and season ------------------
# Element (j1, j2) has column-major index (j2 - 1) * n_r + j1, so variances
# reshape to n_r x n_c (location x season).
var_cor <- d_cor^2                                             # separable correlation
var_cov <- as.vector(outer(diag(fit_cov$C1), diag(fit_cov$C2)))  # separable covariance

season_lab <- factor(rep(seq_len(n_c), each = n_r),
                     levels = seq_len(n_c), labels = c("Spring", "Summer", "Fall"))
var_df <- rbind(
  data.frame(location = rep(seq_len(n_r), n_c), season = season_lab,
             variance = var_cor, Model = "Separable correlation"),
  data.frame(location = rep(seq_len(n_r), n_c), season = season_lab,
             variance = var_cov, Model = "Separable covariance")
)

p_var <- ggplot(var_df, aes(location, variance, shape = Model)) +
  geom_point(size = 1.8) +
  facet_wrap(~ season) +
  scale_shape_manual(values = c("Separable correlation" = 16,
                                "Separable covariance"  = 17)) +
  labs(x = "Location", y = "Estimated variance", shape = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank())

# Saved at the same physical size as Figures 3-4 (inches) so that fonts and
# points scale down consistently when included at \includegraphics width.
ggsave("Figures/var_plot.pdf", p_var, width = 10, height = 3)

cat("Wrote Figures/cov_plot_1.pdf, Figures/cov_plot_2.pdf, Figures/var_plot.pdf\n")
