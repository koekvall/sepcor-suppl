# Test conjecture: n > max(r/c, c/r) suffices (strict inequality, Beta = 0).
# With p = 0, this is just n > max(r/c, c/r).

set.seed(42)
n_rep <- 50000

run_k0 <- function(n, r, c_) {
  # Beta = 0, so p = 0, Q_X = I_n, E_i = Y_i
  calE <- matrix(rnorm(n * r * c_), n, r * c_)
  S_n_diag <- colSums(calE^2)
  check_d <- (0 + sqrt(0 + n * 4 * S_n_diag)) / (2 * n)
  check_E_list <- lapply(1:n, function(i) {
    matrix(calE[i, ] / check_d, r, c_)
  })
  T1 <- Reduce("+", lapply(check_E_list, function(cEi) cEi %*% t(cEi)))
  T2 <- Reduce("+", lapply(check_E_list, function(cEi) t(cEi) %*% cEi))
  c(rank_T1 = qr(T1)$rank, rank_T2 = qr(T2)$rank)
}

# Test cases at the boundary of n > max(r/c, c/r)
# For each (r, c), test:
#   n = floor(max(r/c,c/r))     => equality or below, should FAIL
#   n = floor(max(r/c,c/r)) + 1 => just above if max is integer, AT if not
#   n = ceil(max(r/c,c/r)) + 1  => strictly above, should be OK
cases <- list(
  # r = c: max = 1, boundary at n = 1 vs n = 2
  list(r = 2, c_ = 2, label = "r=c=2, max=1"),
  list(r = 3, c_ = 3, label = "r=c=3, max=1"),
  list(r = 4, c_ = 4, label = "r=c=4, max=1"),
  # c = 2r: max = 2, boundary at n = 2 vs n = 3
  list(r = 2, c_ = 4, label = "r=2,c=4, max=2"),
  list(r = 3, c_ = 6, label = "r=3,c=6, max=2"),
  # c = 3r: max = 3, boundary at n = 3 vs n = 4
  list(r = 2, c_ = 6, label = "r=2,c=6, max=3"),
  # Non-integer max: c/r = 3/2, boundary at n = 1 (below) vs n = 2 (above)
  list(r = 2, c_ = 3, label = "r=2,c=3, max=1.5"),
  list(r = 3, c_ = 2, label = "r=3,c=2, max=1.5"),
  # Non-integer max: c/r = 5/2, boundary at n = 2 (below) vs n = 3 (above)
  list(r = 2, c_ = 5, label = "r=2,c=5, max=2.5"),
  # Non-integer max: c/r = 4/3
  list(r = 3, c_ = 4, label = "r=3,c=4, max=1.33"),
  # Non-integer max: c/r = 5/3
  list(r = 3, c_ = 5, label = "r=3,c=5, max=1.67")
)

cat(sprintf("%-22s %-6s %-6s %-8s %-8s %-8s\n",
            "Case", "n", "n>max?", "T1_fail", "T2_fail", "status"))
cat(strrep("-", 62), "\n")

for (case in cases) {
  r <- case$r
  c_ <- case$c_
  mx <- max(r / c_, c_ / r)

  # Test n values around the boundary
  n_vals <- sort(unique(c(
    max(1, floor(mx)),      # at or below
    ceiling(mx),            # at (if non-integer) or equal
    ceiling(mx) + 1         # strictly above
  )))

  for (n in n_vals) {
    above <- n > mx
    results <- replicate(n_rep, run_k0(n, r, c_))
    T1_fail <- mean(results["rank_T1", ] < r)
    T2_fail <- mean(results["rank_T2", ] < c_)
    status <- if (T1_fail > 0 || T2_fail > 0) "FAIL" else "ok"
    cat(sprintf("%-22s %-6d %-6s %-7.2f%% %-7.2f%% %-8s\n",
                case$label, n, ifelse(above, "yes", "NO"),
                T1_fail * 100, T2_fail * 100, status))
  }
  cat("\n")
}
