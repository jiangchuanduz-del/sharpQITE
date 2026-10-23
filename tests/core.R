library(sharpQITE)
# Fixed unit priorities make these exhaustive checks deterministic. They are
# assignment-independent because these tests keep the original unit order.
fixed_test <- function(z,y,...,tie_order=seq_along(z))
  quantile_test(z,y,...,tie_order=tie_order)
fixed_band <- function(z,y,...,tie_order=seq_along(z))
  quantile_band(z,y,...,tie_order=tie_order)

near <- function(a, b, tol = 1e-11) stopifnot(isTRUE(all.equal(a, b, tolerance = tol, check.attributes = FALSE)))
expect_error <- function(expr) stopifnot(inherits(tryCatch({force(expr); NULL}, error = identity), "error"))
assignments <- function(n, m) utils::combn(n, m, function(i) {z <- integer(n); z[i] <- 1L; z}, simplify = FALSE)

# Independent exhaustive rank-path enumeration, including unbalanced groups
# and thresholds just off the lattice, checks inclusive >= tail orientation.
for (n in 2:10) for (m in seq_len(n - 1L)) {
  zz <- assignments(n, m)
  vals <- vapply(zz, function(z) max(c(0, cumsum(1 - z) / (n - m) - cumsum(z) / m)), numeric(1))
  support <- unique(round(vals * m * (n - m))) / (m * (n - m))
  q <- sort(unique(c(-Inf, -0.01, 0, support, support + 1e-7, 1, 1.01, Inf)))
  truth <- vapply(q, function(a) mean(vals >= a - 1e-13), numeric(1))
  near(rks_tail(q, m, n - m), truth)
}
near(rks_tail(1, 100, 100), exp(-lchoose(200, 100)), tol = 1e-12)
near(rks_tail(1, 100, 100) / exp(-lchoose(200, 100)), 1, tol = 1e-11)
stopifnot(rks_tail(1, 100, 100) > 0) # no 1 - CDF cancellation

# Exact Wilcoxon and Stephenson distributions against independent enumeration.
for (stat in c("wilcoxon", "stephenson")) for (s in c(2, 3, 5)) {
  n <- 7; m <- 3; zlist <- assignments(n, m)
  score <- if (stat == "wilcoxon") 1:n else choose(0:(n-1), s-1) / choose(n-1, s-1)
  vals <- vapply(zlist, function(z) sum(score[z == 1]), numeric(1))
  d <- randomization_distribution(m, n-m, stat, s, method = "exact")
  q <- c(-Inf, unique(vals), Inf)
  expected <- vapply(q, function(a) mean(vals >= a - 1e-12), numeric(1))
  near(randomization_tail(d, q), expected)
}

# The least-favourable configuration equals a brute-force minimum over all
# subsets of at most m-k unconstrained treated effects, with fixed tie priority.
z <- c(1, 0, 1, 0, 1, 0, 1); y <- c(2, 2, 1, 0, 2, -1, 4)
tie <- c(7, 2, 5, 1, 6, 4, 3)
for (stat in c("rks", "wilcoxon", "stephenson")) for (k in 1:4) for (cc in c(-1, 0, 1, 2)) {
  d <- randomization_distribution(4, 3, stat, 3, "exact")
  observed <- fixed_test(z, y, k = k, c = cc, statistic = stat, s = 3,
                            tie_order = tie, distribution = d)$statistic
  options <- list(integer())
  if (4-k > 0) for (j in seq_len(4-k)) options <- c(options, utils::combn(which(z == 1), j, simplify = FALSE))
  vals <- vapply(options, function(idx) {
    yy <- y; yy[z == 1] <- yy[z == 1] - cc; yy[idx] <- -Inf
    rank_statistic(z, yy, stat, 3, tie)
  }, numeric(1))
  near(observed, min(vals))
}

# Inversion checks at the boundary and on both sides, including ties.
for (stat in c("rks", "wilcoxon", "stephenson")) {
  z <- rep(c(0, 1), 5); y <- c(-1, 0, 0, 1, 1, 3, 2, 4, 2, 6)
  b <- fixed_band(z, y, alpha = .2, statistic = stat, s = 3, target = "treated", method = "exact")
  for (k in 1:5) {
    grid <- sort(unique(c(as.vector(outer(y[z==1], y[z==0], "-")), -100, 100)))
    grid <- sort(unique(c(grid, grid - 1e-6, grid + 1e-6)))
    p <- fixed_test(z, y, k = k, c = grid, statistic = stat, s = 3, method = "exact")$p.value
    in_band <- grid > b$lower[k] | (grid == b$lower[k] & b$lower_closed[k])
    stopifnot(identical(in_band, p > .2), all(diff(p) >= -1e-12))
  }
}

# Exhaustive finite-population validity: the treated null is random, so the
# relevant error probability is P(reject AND null true), not P(reject | null).
n <- 8; m <- 4; zlist <- assignments(n, m)
y0 <- c(-2, -1, -1, 0, 1, 1, 2, 3); tau <- c(-1, 0, 2, 0, 1, 3, -2, 1)
covered <- function(b, truth) all((truth > b$lower | (truth == b$lower & b$lower_closed)) &
                                 (truth < b$upper | (truth == b$upper & b$upper_closed)))
for (stat in c("rks", "wilcoxon", "stephenson")) {
  errors <- matrix(FALSE, length(zlist), m)
  cover_t <- cover_all <- cover_two <- logical(length(zlist))
  for (j in seq_along(zlist)) {
    z <- zlist[[j]]; y <- y0 + z * tau
    truth <- sort(tau[z == 1])
    for (k in 1:m) errors[j, k] <- truth[k] <= 0 &&
      fixed_test(z, y, k = k, c = 0, statistic = stat, s = 3, method = "exact")$p.value <= .2
    bt <- fixed_band(z, y, .2, stat, 3, "treated", method = "exact")
    ba <- fixed_band(z, y, .2, stat, 3, "all", method = "exact")
    b2 <- fixed_band(z, y, .2, stat, 3, "all", "two.sided", method = "exact")
    cover_t[j] <- covered(bt, truth)
    cover_all[j] <- covered(ba, sort(tau))
    cover_two[j] <- covered(b2, sort(tau))
  }
  stopifnot(all(colMeans(errors) <= .2 + 1e-12), mean(cover_t) >= .8,
            mean(cover_all) >= .8, mean(cover_two) >= .8)
}

# Sign/label transformations and reflection of upper bands.
z <- zlist[[17]]; y <- y0 + z * tau
near(fixed_test(z, y, k = 3, target = "control")$p.value,
     fixed_test(1-z, -y, k = 3)$p.value)
near(fixed_test(z, y, k = 3, alternative = "less")$p.value,
     fixed_test(z, -y, k = 2)$p.value)
lower <- fixed_band(z, -y, target = "all")
upper <- fixed_band(z, y, target = "all", alternative = "less")
near(upper$upper, -rev(lower$lower))
stopifnot(identical(upper$upper_closed, rev(lower$lower_closed)))

# Kernel CDF agrees with a hand-weighted empirical CDF, is monotone, and has
# exact tails. Check minimization independently over all jump candidates.
x <- c(-1, 0, 1, -1, 0, 1); z <- c(0, 0, 0, 1, 1, 1); y <- c(-2, 0, 2, -1, 2, 4)
fit <- fit_kernel_cdf(x, z, y, bandwidth = 1)
pred <- predict_kernel_cdf(fit, c(0, 1), c(-Inf, -2, 0, 2, 4, Inf), arm = 0)
w <- exp(-x[1:3]^2/2); w <- w / sum(w)
near(pred[1, ], c(0, w[1], sum(w[1:2]), 1, 1, 1))
stopifnot(all(apply(pred, 1, diff) >= -1e-14))
near(predict_kernel_cdf(fit, 1e3, Inf), matrix(1, 1, 1))
effects <- seq(-5, 7, by = .5)
bound <- sharp_cdf_bound(fit, effects, newx = c(-.5, .5))
brute <- vapply(effects, function(e) {
  b <- sort(unique(c(y[z==1], y[z==0] + e, Inf)))
  f1 <- predict_kernel_cdf(fit, c(-.5, .5), b, 1)
  f0 <- predict_kernel_cdf(fit, c(-.5, .5), b-e, 0)
  1 + mean(apply(f1-f0, 1, min))
}, numeric(1))
near(bound, brute); stopifnot(all(diff(bound) >= -1e-12))
adj <- estimate_adjustment(fit, c(-.5, .5), .7)
grid <- sort(unique(as.vector(outer(y[z==1], y[z==0], "-"))))
near(adj$tau_lower, grid[which(sharp_cdf_bound(fit, grid) >= .7)[1]])
fit_point <- fit_kernel_cdf(rep(0, 4), c(0,0,1,1), c(0,0,2,2))
point <- estimate_adjustment(fit_point, c(0,1), .5)
near(point$tau_lower, 2); stopifnot(all(point$adjustment == Inf))

# Fixed infinite adjustments and degenerate splits must remain defined.
stopifnot(is.finite(fixed_test(c(0,0,1,1), c(1,2,3,4), k=2,
                                adjustment=rep(Inf,4))$p.value))
deg <- split_adjusted_band(1:4, c(0,0,1,1), 1:4, folds=c(1,1,2,2))
stopifnot(all(deg$lower == -Inf), all(deg$upper == Inf))
split_coverage <- vapply(zlist, function(z) {
  y <- y0 + z*tau
  b <- split_adjusted_band(seq_len(n), z, y, alpha=.2, folds=rep(1:2,4),tie_order=seq_len(n))
  covered(b, sort(tau))
}, logical(1))
stopifnot(mean(split_coverage) >= .8)

# Monte Carlo plus-one calibration, reproducibility, and local RNG state.
set.seed(91); before <- .Random.seed
d <- randomization_distribution(3, 4, "stephenson", 3, "monte_carlo", B=99, seed=8)
stopifnot(identical(before, .Random.seed))
d2 <- randomization_distribution(3, 4, "stephenson", 3, "monte_carlo", B=99, seed=8)
near(d$values, d2$values)
near(randomization_tail(d, max(d$values)+1), .01)
near(randomization_tail(d, min(d$values)), 1)
near(randomization_tail(d, d$values[60]), (1+sum(d$values >= d$values[60]-1e-12))/100)
# Reusing seeded Monte Carlo draws must give identical inversion decisions.
z <- rep(0:1, 6); y <- seq_len(12) + 2*z
bm <- fixed_band(z, y, alpha=.2, statistic="stephenson", s=3,
                     target="treated", method="monte_carlo", B=99, seed=8)
for (k in seq_len(6)) {
  grid <- seq(-5, 10, by=.5)
  pp <- fixed_test(z, y, k=k, c=grid, statistic="stephenson", s=3,
                      method="monte_carlo", B=99, seed=8)$p.value
  stopifnot(identical(pp > .2, grid > bm$lower[k] |
                      (grid == bm$lower[k] & bm$lower_closed[k])))
}
# s=2 is calibrated exactly even when assignment enumeration is infeasible.
z <- rep(0:1, 40); y <- seq_len(80) + z
near(fixed_test(z,y,beta=.9,statistic="stephenson",s=2)$p.value,
     fixed_test(z,y,beta=.9,statistic="wilcoxon")$p.value)
expect_error(fixed_test(c(0,1), c(0,1), k=2))
expect_error(fixed_test(c(0,1), c(0,1), k=1, beta=.5))
expect_error(fit_kernel_cdf(1:4, c(0,0,1,1), 1:4, bandwidth=0))
expect_error(randomization_distribution(30,30,"stephenson",3,"exact"))

# Simulation smoke test checks both analysis and the design-based summary.
sim <- simulate_comparison(n=12, nsim=3, betas=.8, stephenson_s=c(2,3),
                           training_n=16, B=49, scenarios=simulation_scenarios()[1:2], seed=9)
summary <- summarize_simulation(sim)
stopifnot(nrow(sim$results)==30, nrow(summary)==10,
          all(summary$joint_false_rejection <= summary$rejection_rate),
          all(summary$null_probability[summary$scenario=="null"]==1))
cat("All exhaustive inference, distribution, adjustment, and simulation checks passed.\n")
