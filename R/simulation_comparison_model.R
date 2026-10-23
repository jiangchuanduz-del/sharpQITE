gaussian_covariate_model <- function(
    mu0 = function(x) 3 * x, mu1 = function(x) 4 * x,
    sigma0 = function(x) exp(-0.25 * x),
    sigma1 = function(x) 3 * exp(0.25 * x), xlim = c(-1, 1)) {
  if (!all(vapply(list(mu0, mu1, sigma0, sigma1), is.function, logical(1))))
    stop("mu0, mu1, sigma0 and sigma1 must be vectorized functions.")
  if (!is.numeric(xlim) || length(xlim) != 2L || any(!is.finite(xlim)) || diff(xlim) <= 0)
    stop("xlim must contain two increasing finite endpoints.")
  model <- list(mu0 = mu0, mu1 = mu1, sigma0 = sigma0, sigma1 = sigma1, xlim = xlim)
  .model_parameters(model, seq(xlim[1], xlim[2], length.out = 101))
  structure(model, class = "sharp_gaussian_model")
}

.model_parameters <- function(model, x) {
  out <- lapply(model[c("mu0", "mu1", "sigma0", "sigma1")], function(f) {
    a <- f(x)
    if (length(a) == 1L) a <- rep(a, length(x))
    if (!is.numeric(a) || length(a) != length(x) || any(!is.finite(a)))
      stop("Model functions must return finite numeric vectors of the input length.")
    a
  })
  if (any(out$sigma0 <= 0) || any(out$sigma1 <= 0)) stop("Conditional SDs must be positive.")
  out
}

.draw_gaussian_population <- function(model, n) {
  x <- stats::runif(n, model$xlim[1], model$xlim[2])
  p <- .model_parameters(model, x)
  list(x = x, y0 = p$mu0 + p$sigma0 * stats::rnorm(n),
       y1 = p$mu1 + p$sigma1 * stats::rnorm(n))
}

.gauss_uniform <- function(n, limits) {
  k <- seq_len(n - 1L)
  a <- matrix(0, n, n)
  a[cbind(k, k+1L)] <- a[cbind(k+1L, k)] <- k / sqrt(4*k^2 - 1)
  e <- eigen(a, symmetric = TRUE)
  ord <- order(e$values)
  list(x = mean(limits) + diff(limits)/2 * e$values[ord],
       w = e$vectors[1, ord]^2)
}

# All stationary points of F1(b|x)-F0(b-effect|x), plus the zero tail limit.
.normal_cdf_minimum <- function(effect, p) {
  a <- .5 * (1/p$sigma0^2 - 1/p$sigma1^2)
  b <- p$mu1/p$sigma1^2 - (effect+p$mu0)/p$sigma0^2
  cc <- .5*((effect+p$mu0)^2/p$sigma0^2 - p$mu1^2/p$sigma1^2) +
    log(p$sigma0/p$sigma1)
  roots1 <- roots2 <- rep(Inf, length(a))
  quadratic <- abs(a) > 1e-12 * pmax(1/p$sigma0^2, 1/p$sigma1^2)
  disc <- pmax(0, b[quadratic]^2 - 4*a[quadratic]*cc[quadratic])
  # Stable quadratic roots, including a repeated root at zero.
  q <- -.5 * (b[quadratic] + ifelse(b[quadratic] < 0, -1, 1) * sqrt(disc))
  roots1[quadratic] <- q / a[quadratic]
  roots2[quadratic] <- ifelse(q == 0, -b[quadratic]/(2*a[quadratic]), cc[quadratic]/q)
  linear <- !quadratic & abs(b) > 1e-14
  roots1[linear] <- -cc[linear]/b[linear]
  difference <- function(v) stats::pnorm((v-p$mu1)/p$sigma1) -
    stats::pnorm((v-effect-p$mu0)/p$sigma0)
  d1 <- difference(roots1); d2 <- difference(roots2)
  list(minimum = pmin(0, d1, d2),
       adjustment = ifelse(pmin(d1,d2) >= 0, Inf, ifelse(d1 <= d2, roots1, roots2)))
}

.mixture_quantile_bound <- function(p, weights, betas, outcome_grid_n) {
  yrange <- range(c(p$mu0 - 11*p$sigma0, p$mu0 + 11*p$sigma0,
                    p$mu1 - 11*p$sigma1, p$mu1 + 11*p$sigma1))
  ygrid <- seq(yrange[1], yrange[2], length.out = outcome_grid_n)
  cdf0 <- function(t) vapply(t, function(a) sum(weights * stats::pnorm((a-p$mu0)/p$sigma0)), numeric(1))
  cdf1 <- function(t) vapply(t, function(a) sum(weights * stats::pnorm((a-p$mu1)/p$sigma1)), numeric(1))
  q0 <- stats::approxfun(cdf0(ygrid), ygrid, ties = "ordered", rule = 2)
  q1 <- stats::approxfun(cdf1(ygrid), ygrid, ties = "ordered", rule = 2)
  vapply(betas, function(beta) {
    # Search all grid-local maxima, then refine each; no unimodality assumption.
    u <- beta * seq(1e-7, 1-1e-7, length.out = 501)
    value <- q1(u) - q0(1+u-beta)
    peaks <- which(value[2:500] >= value[1:499] & value[2:500] >= value[3:501]) + 1L
    peaks <- unique(c(peaks, which.max(value)))
    max(vapply(peaks, function(j) {
      stats::optimize(function(v) q1(v)-q0(1+v-beta),
                       interval = u[c(max(1,j-1),min(501,j+1))], maximum = TRUE,
                       tol = 1e-10)$objective
    }, numeric(1)))
  }, numeric(1))
}

superpopulation_bounds <- function(model = gaussian_covariate_model(),
                                    betas = seq(.005, .995, by = .005),
                                    quadrature_n = 161L, outcome_grid_n = 8001L,
                                    fixed_beta = .5) {
  if (!inherits(model, "sharp_gaussian_model")) stop("Use gaussian_covariate_model().")
  if (!is.numeric(betas) || !length(betas) || anyNA(betas) || any(betas <= 0 | betas >= 1))
    stop("betas must be in (0,1).")
  fixed_beta <- .probability(fixed_beta, "fixed_beta")
  quadrature_n <- .integer_scalar(quadrature_n, "quadrature_n", 8)
  outcome_grid_n <- .integer_scalar(outcome_grid_n, "outcome_grid_n", 101)
  quadrature <- .gauss_uniform(quadrature_n, model$xlim)
  p <- .model_parameters(model, quadrature$x)
  marginal <- .mixture_quantile_bound(p, quadrature$w, betas, outcome_grid_n)
  fbar <- function(e) 1 + sum(quadrature$w * .normal_cdf_minimum(e, p)$minimum)
  yrange <- range(c(p$mu0 - 11*p$sigma0, p$mu0 + 11*p$sigma0,
                    p$mu1 - 11*p$sigma1, p$mu1 + 11*p$sigma1))
  erange <- c(-diff(yrange), diff(yrange))
  covariate <- vapply(betas, function(beta)
    stats::uniroot(function(e) fbar(e)-beta, erange, tol=1e-8, extendInt="upX")$root,
    numeric(1))
  fixed_effect <- stats::uniroot(function(e) fbar(e)-fixed_beta, erange,
                                tol=1e-8, extendInt="upX")$root
  shift <- .normal_cdf_minimum(fixed_effect,p)$adjustment
  # When only the zero tail limit minimizes the objective, use a finite
  # 11-SD tail proxy. This case does not occur in the default study model.
  nonfinite <- !is.finite(shift)
  shift[nonfinite] <- pmin(p$mu1-11*p$sigma1,
                          p$mu0+fixed_effect-11*p$sigma0)[nonfinite]
  adjusted <- p
  adjusted$mu0 <- p$mu0-shift
  adjusted$mu1 <- p$mu1-shift
  fixed <- .mixture_quantile_bound(adjusted,quadrature$w,betas,outcome_grid_n)
  ans <- data.frame(beta = betas, marginal = marginal, covariate = covariate,
                    covariate_fixed = fixed)
  attr(ans, "quadrature_n") <- quadrature_n
  attr(ans, "outcome_grid_n") <- outcome_grid_n
  attr(ans, "fixed_beta") <- fixed_beta
  attr(ans, "fixed_effect") <- fixed_effect
  attr(ans, "tail_proxy_nodes") <- sum(nonfinite)
  ans
}
