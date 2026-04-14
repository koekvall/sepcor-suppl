# Table: Convergence and uniqueness thresholds
# For each (r, c) with 2 <= r <= c <= 10, find:
#   n1: smallest n for which all datasets yield at least one converged solution
#       (empirical analogue of a.s. existence)
#   n2: smallest n for which all datasets yield a unique converged solution
#       (empirical analogue of a.s. uniqueness)
# Exact lower bounds from Derksen & Makam (2021), Theorem 1.3.
# Search starts at the Derksen lower bound; no need to search below it.
#
# Usage: Rscript scripts/tab_boundary.R
# Output: printed LaTeX table

library(sepcor)
library(parallel)

set.seed(2028)
n_datasets <- 50
n_starts   <- 20
n_cores    <- max(1L, detectCores() - 1L)

# ---- Derksen & Makam (2021) Theorem 1.3 ------------------------------------
gcd <- function(a, b) if (b == 0L) a else Recall(b, a %% b)

derksen_bounds <- function(r, cc) {
  p <- min(r, cc); q <- max(r, cc)
  if (p == q) {
    return(if (p == 1L) c(n1 = 1L, n2 = 1L) else c(n1 = 1L, n2 = 3L))
  }
  d      <- gcd(p, q)
  r_tilde <- (p^2 + q^2 - d^2) / (p * q)
  if (isTRUE(all.equal(r_tilde, round(r_tilde)))) {   # integer case
    n1 <- as.integer(round(r_tilde))
    n2 <- if (d == 1L) n1 else n1 + 1L
  } else {
    val <- ceiling((p^2 + q^2) / (p * q))
    n1  <- val; n2 <- val
  }
  c(n1 = n1, n2 = n2)
}

# ---- One (r, c, n) check ---------------------------------------------------
check_n <- function(r, cc, n) {
  q   <- r * cc
  U0  <- 0.5^abs(outer(seq_len(cc), seq_len(cc), "-"))
  V0  <- 0.5^abs(outer(seq_len(r),  seq_len(r),  "-"))
  Sig_c <- kronecker(t(chol(U0)), t(chol(V0)))

  res <- mclapply(seq_len(n_datasets), function(d) {
    E <- Sig_c %*% matrix(rnorm(q * n), ncol = n)
    S <- tcrossprod(E) / n
    params    <- vector("list", n_starts)
    any_conv  <- FALSE
    all_conv  <- TRUE
    for (s in seq_len(n_starts)) {
      W_init <- if (s == 1L) diag(S) else
        (sqrt(pmax(diag(S), 1e-6)) * exp(rnorm(q, 0, 0.5)))^2
      fit <- suppressWarnings(
        sepcor:::sepcor_rcpp(E, W_init, r, 1e-8, 3000L, FALSE, 0)
      )
      if (fit$info != 0L) {
        all_conv <- FALSE
      } else {
        any_conv <- TRUE
        # Store parameter vector (C1 upper tri, C2 upper tri, D diagonal)
        params[[s]] <- c(fit$C1[upper.tri(fit$C1)],
                         fit$C2[upper.tri(fit$C2)],
                         as.vector(fit$D))
      }
    }
    is_unique <- if (!all_conv) FALSE else {
      # All starts converged: check all parameter vectors agree
      ref <- params[[1L]]
      all(sapply(params[-1L], function(p) max(abs(p - ref)) < 1e-4))
    }
    c(any_conv = any_conv, all_conv = all_conv, unique = is_unique)
  }, mc.cores = n_cores)

  frac_any    <- mean(sapply(res, `[`, "any_conv"))
  frac_all    <- mean(sapply(res, `[`, "all_conv"))
  frac_unique <- mean(sapply(res, `[`, "unique"))
  c(frac_any = frac_any, frac_all = frac_all, frac_unique = frac_unique)
}

# ---- Selected (r, c) pairs: covers all Derksen formula cases, range of dims -
pairs <- list(
  list(r = 2L,  cc = 2L),
  list(r = 2L,  cc = 3L),
  list(r = 2L,  cc = 6L),
  list(r = 2L,  cc = 9L),
  list(r = 3L,  cc = 3L),
  list(r = 3L,  cc = 6L),
  list(r = 5L,  cc = 5L),
  list(r = 5L,  cc = 6L),
  list(r = 5L,  cc = 10L)
)

results <- data.frame()

for (cfg in pairs) {
  r  <- cfg$r
  cc <- cfg$cc
  dk <- derksen_bounds(r, cc)

  # Search from Derksen lower bound; upper bound is generous
  n_min <- max(1L, dk["n1"])
  n_max <- max(dk["n2"] * 3L + 5L, 20L)

  sepcor_n1 <- NA_integer_
  sepcor_n2 <- NA_integer_

  for (n in n_min:n_max) {
    res <- check_n(r, cc, n)
    cat(sprintf("r=%2d c=%2d n=%2d  any=%.2f  all=%.2f  unique=%.2f\n",
                r, cc, n, res["frac_any"], res["frac_all"],
                res["frac_unique"]))

    if (is.na(sepcor_n1) && res["frac_all"] == 1) sepcor_n1 <- n
    if (is.na(sepcor_n2) && res["frac_unique"] == 1) {
      sepcor_n2 <- n
      break
    }
  }

  results <- rbind(results, data.frame(
    r = r, c = cc,
    dk_n1 = dk["n1"], dk_n2 = dk["n2"],
    sepcor_n1 = sepcor_n1, sepcor_n2 = sepcor_n2
  ))
}

# ---- Print LaTeX table ------------------------------------------------------
cat("\n=== LaTeX table ===\n")
cat("\\begin{tabular}{rr cc cc}\n")
cat(" & & \\multicolumn{2}{c}{Sep.\\ covariance} & \\multicolumn{2}{c}{Sep.\\ correlation} \\\\\n")
cat("$r$ & $c$ & $n_1$ & $n_2$ & $\\hat{n}_1$ & $\\hat{n}_2$ \\\\\n")
cat("\\hline\n")
for (i in seq_len(nrow(results))) {
  cat(sprintf("%d & %d & %d & %d & %d & %d \\\\\n",
              results$r[i], results$c[i],
              results$dk_n1[i], results$dk_n2[i],
              results$sepcor_n1[i], results$sepcor_n2[i]))
}
cat("\\end{tabular}\n")

saveRDS(results, "scripts/tab_boundary_results.rds")
cat("Saved results to scripts/tab_boundary_results.rds\n")
