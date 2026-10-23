simulation_scenarios <- function() {
  make <- function(kind) {
    force(kind)
    function(n) {
      x <- stats::rnorm(n)
      y0 <- 2 * x + stats::rnorm(n)
      tau <- switch(kind,
                    null = rep(0, n),
                    constant = rep(2.5, n),
                    sparse = 6 * as.numeric(x > stats::qnorm(0.75)),
                    heterogeneous = 1 + 1.5 * x)
      list(x = matrix(x, ncol = 1), y0 = y0, y1 = y0 + tau)
    }
  }
  stats::setNames(lapply(c("null", "constant", "sparse", "heterogeneous"), make),
           c("null", "constant", "sparse", "heterogeneous"))
}

.check_population <- function(pop, n) {
  if (!is.list(pop) || !all(c("x", "y0", "y1") %in% names(pop)) ||
      !is.numeric(pop$y0) || !is.numeric(pop$y1) ||
      length(pop$y0) != n || length(pop$y1) != n ||
      any(!is.finite(c(pop$y0, pop$y1))))
    stop("A scenario must return a list with finite numeric y0, y1 and covariates x.")
  pop$x <- .as_x(pop$x)
  if (nrow(pop$x) != n) stop("Scenario x has an incorrect number of rows.")
  pop
}

simulate_comparison <- function(n = 80L, n1 = floor(n / 2), nsim = 500L,
                                betas = c(0.5, 0.8, 0.9), c = 0,
                                stephenson_s = c(2L, 5L, 10L, 20L),
                                training_n = 240L, alpha = 0.05, B = 19999L,
                                scenarios = simulation_scenarios(), seed = 20261003L) {
  n <- .integer_scalar(n, "n", 4); n1 <- .integer_scalar(n1, "n1", 1, n - 1)
  nsim <- .integer_scalar(nsim, "nsim")
  training_n <- .integer_scalar(training_n, "training_n", 4)
  alpha <- .probability(alpha, "alpha")
  if (!is.numeric(betas) || !length(betas) || anyNA(betas) || any(betas <= 0 | betas >= 1))
    stop("betas must be a nonempty vector in (0,1).")
  if (!is.numeric(c) || length(c) != 1L || !is.finite(c)) stop("c must be one finite threshold.")
  if (!is.numeric(stephenson_s) || anyNA(stephenson_s) ||
      any(stephenson_s != floor(stephenson_s) | stephenson_s < 2 | stephenson_s > n))
    stop("stephenson_s must contain integers between 2 and n.")
  if (!is.list(scenarios) || !length(scenarios) || is.null(names(scenarios)) ||
      any(!nzchar(names(scenarios))) || anyDuplicated(names(scenarios)) ||
      !all(vapply(scenarios, is.function, logical(1)))) stop("scenarios must be a named list of functions.")
  .with_seed(seed, {
    specs <- c(list(list(name = "RKS", statistic = "rks", s = 2L),
                    list(name = "RKS adjusted", statistic = "rks", s = 2L),
                    list(name = "Wilcoxon", statistic = "wilcoxon", s = 2L)),
               lapply(unique(stephenson_s), function(s)
                 list(name = paste0("Stephenson s=", s), statistic = "stephenson", s = s)))
    nulls <- lapply(specs, function(spec) {
      randomization_distribution(n1, n - n1, spec$statistic, spec$s, B = B)
    })
    rows <- vector("list", length(scenarios) * nsim * length(betas) * length(specs))
    cursor <- 0L; populations <- list(); training <- list(); fitted <- list()
    for (scenario in names(scenarios)) {
      # Draw once, then hold both potential-outcome vectors fixed throughout
      # all assignments, matching the finite-population design-based target.
      pop <- .check_population(scenarios[[scenario]](n), n)
      train <- .check_population(scenarios[[scenario]](training_n), training_n)
      z_train <- integer(training_n); z_train[sample.int(training_n, floor(training_n / 2))] <- 1L
      y_train <- ifelse(z_train == 1, train$y1, train$y0)
      fit <- fit_kernel_cdf(train$x, z_train, y_train)
      adj <- lapply(betas, function(beta) estimate_adjustment(fit, pop$x, beta))
      populations[[scenario]] <- pop
      training[[scenario]] <- list(population = train, z = z_train, y = y_train)
      fitted[[scenario]] <- adj
      for (iteration in seq_len(nsim)) {
        z <- integer(n); z[sample.int(n, n1)] <- 1L
        y <- ifelse(z == 1, pop$y1, pop$y0)
        tau_t <- sort((pop$y1 - pop$y0)[z == 1])
        for (j in seq_along(betas)) {
          k <- ceiling(.snap_integer(n1 * betas[j])); truth <- tau_t[k]
          for (l in seq_along(specs)) {
            spec <- specs[[l]]
            st <- spec$statistic
            start <- proc.time()[[3L]]
            test <- quantile_test(z, y, k = k, c = c, statistic = st, s = spec$s,
                                  adjustment = if (spec$name == "RKS adjusted") adj[[j]]$adjustment else NULL,
                                  distribution = nulls[[l]])
            elapsed <- proc.time()[[3L]] - start
            cursor <- cursor + 1L
            rows[[cursor]] <- data.frame(scenario = scenario, iteration = iteration,
                                         method = spec$name, beta = betas[j], k = k,
                                         threshold = c, true_quantile = truth,
                                         null_true = truth <= c, p_value = test$p.value,
                                         reject = test$p.value <= alpha,
                                         seconds = elapsed, calibration = test$method)
          }
        }
      }
    }
    structure(list(results = do.call(rbind, rows), populations = populations,
                   training = training, adjustments = fitted,
                   settings = list(n = n, n1 = n1, nsim = nsim, betas = betas,
                                   c = c, alpha = alpha, B = B, training_n = training_n,
                                   seed = seed, design = "fixed populations; repeated complete randomization")),
              class = "sharp_simulation")
  })
}

summarize_simulation <- function(simulation, conf_level = 0.95) {
  .probability(conf_level, "conf_level")
  dat <- if (inherits(simulation, "sharp_simulation")) simulation$results else simulation
  if (!is.data.frame(dat) || !all(c("scenario", "method", "beta", "reject", "null_true", "seconds") %in% names(dat)))
    stop("Expected output from simulate_comparison() or its results data frame.")
  zcrit <- stats::qnorm((1 + conf_level) / 2)
  groups <- split(seq_len(nrow(dat)), interaction(dat$scenario, dat$method, dat$beta, drop = TRUE))
  ans <- lapply(groups, function(idx) {
    d <- dat[idx, , drop = FALSE]; N <- nrow(d); p <- mean(d$reject)
    center <- (p + zcrit^2 / (2 * N)) / (1 + zcrit^2 / N)
    half <- zcrit * sqrt(p * (1 - p) / N + zcrit^2 / (4 * N^2)) / (1 + zcrit^2 / N)
    data.frame(scenario = d$scenario[1], method = d$method[1], beta = d$beta[1],
                nsim = N, rejection_rate = p, mcse = sqrt(p * (1 - p) / N),
                ci_low = max(0, center - half), ci_high = min(1, center + half),
                null_probability = mean(d$null_true),
                joint_false_rejection = mean(d$reject & d$null_true),
                power_when_null_false = if (all(d$null_true)) NA_real_ else mean(d$reject[!d$null_true]),
                mean_seconds = mean(d$seconds))
  })
  ans <- do.call(rbind, ans)
  ans <- ans[order(ans$scenario, ans$beta, ans$method), ]
  rownames(ans) <- NULL
  ans
}
