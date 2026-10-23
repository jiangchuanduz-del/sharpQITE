.as_x <- function(x, p = NULL) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.numeric(x) && is.null(dim(x)))
    x <- if (!is.null(p) && p > 1L) matrix(x, nrow = 1L) else matrix(x, ncol = 1L)
  if (!is.matrix(x) || !is.numeric(x) || !nrow(x) || !ncol(x) || any(!is.finite(x)))
    stop("x must be a finite numeric vector or matrix; encode factors explicitly.")
  if (!is.null(p) && ncol(x) != p) stop("newx has the wrong number of columns.")
  x
}

fit_kernel_cdf <- function(x, z, y, bandwidth = NULL) {
  .check_data(z, y)
  x <- .as_x(x)
  if (nrow(x) != length(y)) stop("x must have one row per training unit.")
  p <- ncol(x)
  if (is.null(bandwidth)) {
    scale <- apply(x, 2, stats::sd)
    scale[!is.finite(scale) | scale == 0] <- 1
    bandwidth <- scale * min(sum(z), sum(1 - z))^(-1 / (p + 4))
  }
  if (length(bandwidth) == 1L) bandwidth <- rep(bandwidth, p)
  if (!is.numeric(bandwidth) || length(bandwidth) != p || any(!is.finite(bandwidth)) ||
      any(bandwidth <= 0)) stop("bandwidth must be positive, with one value per covariate.")
  structure(list(x = x, z = z, y = y, bandwidth = bandwidth,
                 kernel = "Gaussian product kernel", call = match.call()),
            class = "sharp_kernel_cdf")
}

.kernel_weights <- function(fit, newx, arm) {
  idx <- which(fit$z == arm)
  idx <- idx[order(fit$y[idx], idx)]
  train <- fit$x[idx, , drop = FALSE]
  distance <- matrix(0, nrow(newx), length(idx))
  for (j in seq_len(ncol(newx)))
    distance <- distance + (outer(newx[, j], train[, j], "-") / fit$bandwidth[j])^2
  if (any(!is.finite(distance))) stop("Kernel distances overflow; rescale x or bandwidth.")
  # Subtracting the row maximum log-weight prevents all-zero weights far from data.
  distance <- distance - apply(distance, 1, min)
  w <- exp(-distance / 2)
  w <- w / rowSums(w)
  list(weights = w, y = fit$y[idx], effective_n = 1 / rowSums(w^2))
}

.check_fit <- function(fit) {
  if (!inherits(fit, "sharp_kernel_cdf")) stop("fit must come from fit_kernel_cdf().")
}

predict_kernel_cdf <- function(fit, newx, thresholds, arm = 1L) {
  .check_fit(fit)
  newx <- .as_x(newx, ncol(fit$x))
  arm <- .integer_scalar(arm, "arm", 0, 1)
  if (!is.numeric(thresholds) || !length(thresholds) || anyNA(thresholds))
    stop("thresholds must be a nonempty numeric vector without missing values.")
  a <- .kernel_weights(fit, newx, arm)
  a$weights %*% outer(a$y, thresholds, "<=")
}

.kernel_pair <- function(fit, newx) {
  a <- .kernel_weights(fit, newx, 1L)
  b <- .kernel_weights(fit, newx, 0L)
  list(treated = a, control = b,
       control_cdf = b$weights %*% outer(b$y, b$y, "<="),
       differences = outer(a$y, b$y, "-"))
}

.conditional_min <- function(pair, effect) {
  # The difference can decrease only at a shifted control-outcome jump.
  # +Inf is an additional candidate with value zero, as in equation (36).
  delta <- pair$treated$weights %*% (pair$differences <= effect) - pair$control_cdf
  pick <- max.col(-delta, ties.method = "first")
  value <- delta[cbind(seq_len(nrow(delta)), pick)]
  use_infinity <- value >= -64 * .Machine$double.eps
  value <- pmax(-1, pmin(0, value))
  value[use_infinity] <- 0
  adjustment <- pair$control$y[pick] + effect
  adjustment[use_infinity] <- Inf
  list(minimum = value, adjustment = adjustment)
}

sharp_cdf_bound <- function(fit, effect, newx = fit$x) {
  .check_fit(fit)
  newx <- .as_x(newx, ncol(fit$x))
  if (!is.numeric(effect) || !length(effect) || any(!is.finite(effect)))
    stop("effect must contain finite thresholds.")
  pair <- .kernel_pair(fit, newx)
  vapply(effect, function(e) 1 + mean(.conditional_min(pair, e)$minimum), numeric(1))
}

estimate_adjustment <- function(fit, newx, beta, integration_x = fit$x,
                                max_breaks = 1000000L) {
  .check_fit(fit)
  beta <- .probability(beta, "beta")
  newx <- .as_x(newx, ncol(fit$x))
  integration_x <- .as_x(integration_x, ncol(fit$x))
  max_breaks <- .integer_scalar(max_breaks, "max_breaks")
  n1 <- sum(fit$z); n0 <- sum(1 - fit$z)
  if (as.double(n1) * n0 > max_breaks)
    stop("Training cross-arm grid exceeds max_breaks; use a smaller training sample or raise the cap.")
  pair <- .kernel_pair(fit, integration_x)
  grid <- sort(unique(as.vector(pair$differences)))
  if (any(!is.finite(grid))) stop("Outcome differences overflow; rescale training outcomes.")
  fbar <- function(e) 1 + mean(.conditional_min(pair, e)$minimum)
  lo <- 1L; hi <- length(grid)
  while (lo < hi) {
    mid <- floor((lo + hi) / 2)
    if (fbar(grid[mid]) >= beta) hi <- mid else lo <- mid + 1L
  }
  tau_lower <- grid[lo]
  pred_pair <- .kernel_pair(fit, newx)
  pred <- .conditional_min(pred_pair, tau_lower)
  structure(list(adjustment = pred$adjustment, tau_lower = tau_lower,
                 beta = beta, estimated_cdf = fbar(tau_lower),
                 conditional_minimum = pred$minimum,
                 effective_n_treated = pred_pair$treated$effective_n,
                 effective_n_control = pred_pair$control$effective_n,
                 bandwidth = fit$bandwidth, call = match.call()),
            class = "sharp_adjustment")
}

.trivial_band <- function(n, alpha, alternative) {
  .new_band(data.frame(k = seq_len(n), beta = seq_len(n) / n,
                       lower = -Inf, upper = Inf,
                       lower_closed = FALSE, upper_closed = FALSE),
            alpha, "all", "rks", alternative, alpha, "uninformative", TRUE)
}

split_adjusted_band <- function(x, z, y, beta_adjust = 0.8, alpha = 0.05,
                                alternative = c("greater", "less", "two.sided"),
                                folds = NULL, bandwidth = NULL, seed = NULL,
                                tie_order = NULL) {
  .check_data(z, y)
  x <- .as_x(x)
  if (nrow(x) != length(z)) stop("x must have one row per unit.")
  beta_adjust <- .probability(beta_adjust, "beta_adjust")
  alpha <- .probability(alpha, "alpha")
  alternative <- match.arg(alternative)
  n <- length(z)
  if (n < 4L) stop("Sample splitting requires at least four units.")
  .with_seed(seed, {
  if (is.null(folds)) folds <- sample(rep(1:2, length.out = n))
  if (length(folds) != n || anyNA(folds) || !all(folds %in% 1:2) ||
      length(unique(folds)) != 2L) stop("folds must contain both 1 and 2, one per unit.")
  tie_order <- .tie_order(tie_order, n)
  bands <- vector("list", 2L); adjustments <- vector("list", 2L)
  for (f in 1:2) {
    test <- which(folds == f); train <- which(folds != f)
    if (length(unique(z[test])) < 2L || length(unique(z[train])) < 2L) {
      # Never redraw a split based on the realized assignments. Returning an
      # uninformative band preserves unconditional design-based coverage.
      bands[[f]] <- .trivial_band(length(test), alpha / 2, alternative)
    } else {
      fit <- fit_kernel_cdf(x[train, , drop = FALSE], z[train], y[train], bandwidth)
      adj <- estimate_adjustment(fit, x[test, , drop = FALSE], beta_adjust)
      bands[[f]] <- quantile_band(z[test], y[test], alpha / 2, "rks",
                                  alternative = alternative, adjustment = adj$adjustment,
                                  tie_order = tie_order[test], method = "exact")
      adjustments[[f]] <- adj
    }
  }
  low <- .pool_lower(lapply(bands, function(b) data.frame(bound = b$lower, closed = b$lower_closed)))
  neg_up <- .pool_lower(lapply(bands, function(b) data.frame(bound = -b$upper, closed = b$upper_closed)))
  ans <- .new_band(data.frame(k = seq_len(n), beta = seq_len(n) / n,
                              lower = low$bound, upper = -rev(neg_up$bound),
                              lower_closed = low$closed, upper_closed = rev(neg_up$closed)),
                   alpha, "all", "split-adjusted rks", alternative,
                   alpha / if (alternative == "two.sided") 8 else 4,
                   "exact_dp", TRUE, tie_order = tie_order, seed = seed)
  attr(ans, "folds") <- folds
  attr(ans, "fold_bands") <- bands
  attr(ans, "adjustments") <- adjustments
  attr(ans, "beta_adjust") <- beta_adjust
  ans
  })
}
