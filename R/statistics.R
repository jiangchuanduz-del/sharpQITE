# Resolve a random tie priority once, or reuse a supplied independent priority.
rank_statistic <- function(z, y, statistic = "rks", s = 2L, tie_order = NULL,
                           seed = NULL) {
  .check_data(z, y, allow_infinite = TRUE)
  statistic <- .statistic_name(statistic)
  .with_seed(seed, {
    tie_order <- .tie_order(tie_order, length(z))
    ord <- order(y, tie_order)
    structure(.stat_order(z[ord], statistic, s), tie_order = tie_order, seed = seed)
  })
}

.stat_order <- function(z_ordered, statistic, s) {
  n <- length(z_ordered)
  m <- as.double(sum(z_ordered))
  if (statistic == "rks") {
    # Integer lattice avoids cancellation in empirical CDF differences.
    return(max(c(0, cumsum(1 - z_ordered) * m -
                   cumsum(z_ordered) * (n - m))) / (m * (n - m)))
  }
  sum(.scores(n, statistic, s)[z_ordered == 1])
}

# Probability of ever reaching the inclusive upper boundary on a uniform
# rank-order path. Backwards, diagonal recursion: O(n1*n0) time, O(n1) space.
.rks_boundary_tail <- function(h, n1, n0) {
  n1 <- as.double(n1); n0 <- as.double(n0)
  if (h <= 0) return(1)
  if (h > n1 * n0) return(0)
  if (min(n1,n0)==1) return((n1+n0-ceiling(h))/(n1+n0))
  if (lchoose(n1+n0,n1) < log(2^53)-1e-12) {
    # Count paths using exactly representable integers, then divide once.
    # This avoids accumulating roundoff at attainable small-sample tails.
    all_next <- hit_next <- numeric(n1+1L)
    all_next[n1+1L] <- 1
    n <- n1+n0
    for(total in seq.int(n-1L,0L)) {
      i <- seq.int(max(0,total-n0),min(n1,total)); j <- total-i
      all_current <- hit_current <- numeric(n1+1L)
      take_t <- i<n1; take_c <- j<n0
      all_current[i[take_t]+1L] <- all_next[i[take_t]+2L]
      hit_current[i[take_t]+1L] <- hit_next[i[take_t]+2L]
      idx <- i[take_c]+1L
      all_current[idx] <- all_current[idx]+all_next[idx]
      hit_current[idx] <- hit_current[idx]+hit_next[idx]
      hit <- i[j*n1-i*n0>=h]+1L
      hit_current[hit] <- all_current[hit]
      all_next <- all_current; hit_next <- hit_current
    }
    return(hit_next[1L]/all_next[1L])
  }
  if (n1 == n0) {
    # Reflection principle for a balanced random-walk bridge. Equivalent to
    # the DP: choose(2m,m+r)/choose(2m,m), r = ceiling(h/m).
    r <- ceiling(.snap_integer(h / n1))
    j <- seq_len(r)
    return(exp(sum(log1p(-(2 * j - 1) / (n1 + j)))))
  }
  next_prob <- numeric(n1 + 1L)
  n <- n1 + n0
  for (total in seq.int(n - 1L, 0L)) {
    i <- seq.int(max(0L, total - n0), min(n1, total))
    j <- total - i
    current <- numeric(n1 + 1L)
    hit <- j * n1 - i * n0 >= h
    current[i[hit] + 1L] <- 1
    u <- i[!hit]
    v <- j[!hit]
    if (length(u)) {
      p <- numeric(length(u))
      can_treat <- u < n1
      p[can_treat] <- (n1 - u[can_treat]) / (n - total) *
        next_prob[u[can_treat] + 2L]
      can_control <- v < n0
      p[can_control] <- p[can_control] +
        (n0 - v[can_control]) / (n - total) * next_prob[u[can_control] + 1L]
      current[u + 1L] <- p
    }
    next_prob <- current
  }
  min(1, max(0, next_prob[1L]))
}

rks_tail <- function(q, n1, n0) {
  n1 <- .integer_scalar(n1, "n1")
  n0 <- .integer_scalar(n0, "n0")
  if (!is.numeric(q) || anyNA(q)) stop("q must be numeric without missing values.")
  q <- as.double(q)
  vapply(q, function(a) {
    if (a <= 0) return(1)
    if (a > 1) return(0)
    h <- ceiling(.snap_integer(a * n1 * n0))
    .rks_boundary_tail(h, n1, n0)
  }, numeric(1))
}

# Smallest integer numerator whose exact inclusive RKS tail rejects.
.rks_critical_height <- function(n1, n0, alpha) {
  n1 <- .integer_scalar(n1,"n1"); n0 <- .integer_scalar(n0,"n0")
  alpha <- .probability(alpha,"alpha")
  lo <- 1; hi <- as.double(n1)*n0
  if (.rks_boundary_tail(hi,n1,n0)>alpha) return(Inf)
  while(lo<hi) {
    mid <- floor((lo+hi)/2)
    if (.rks_boundary_tail(mid,n1,n0)<=alpha) hi <- mid else lo <- mid+1
  }
  lo
}

# Reject for a statistic strictly above ordered reference value j. The observed
# statistic contributes the plus one, including when no rejection is possible.
.mc_critical_index <- function(B, alpha) {
  B <- .integer_scalar(B,"B"); alpha <- .probability(alpha,"alpha")
  allowed <- floor(alpha*(B+1))
  # Compare the actual represented tail, not a snapped product. Near a grid
  # point, multiplication and division can round in different directions.
  while (allowed>0 && allowed/(B+1)>alpha) allowed <- allowed-1
  while (allowed<B+1 && (allowed+1)/(B+1)<=alpha) allowed <- allowed+1
  B - allowed + 1
}

randomization_distribution <- function(n1, n0, statistic = "rks", s = 2L,
                                       method = c("auto", "exact", "monte_carlo"),
                                       B = 9999L, seed = NULL,
                                       max_enumerate = 200000L) {
  n1 <- .integer_scalar(n1, "n1")
  n0 <- .integer_scalar(n0, "n0")
  statistic <- .statistic_name(statistic)
  method <- match.arg(method)
  n <- n1 + n0
  if (statistic == "stephenson") s <- .integer_scalar(s, "s", 2, n)
  max_enumerate <- .integer_scalar(max_enumerate, "max_enumerate")
  ans <- list(n1 = n1, n0 = n0, statistic = statistic, s = s,
              cache = new.env(parent = emptyenv()), values = NULL, B = 0L)
  if (method != "monte_carlo" && statistic == "rks") {
    ans$method <- "exact_dp"
  } else if (method != "monte_carlo" && (statistic == "wilcoxon" ||
                                        (statistic == "stephenson" && s == 2L))) {
    ans$method <- "exact_wilcoxon"
  } else {
    possible <- lchoose(n, n1) <= log(max_enumerate) + 1e-12
    if (method == "exact" && !possible)
      stop("Exact enumeration exceeds max_enumerate; choose monte_carlo or increase the cap.")
    exact <- method != "monte_carlo" && possible
    scores <- if (statistic == "rks") NULL else .scores(n, statistic, s)
    eval_subset <- function(idx) {
      if (statistic != "rks") return(sum(scores[idx]))
      zz <- integer(n); zz[idx] <- 1L
      .stat_order(zz, "rks", s)
    }
    if (exact) {
      ans$values <- sort(utils::combn(n, n1, FUN = eval_subset))
      ans$method <- "exact_enumeration"
    } else {
      B <- .integer_scalar(B, "B")
      ans$values <- sort(.with_seed(seed, replicate(B, eval_subset(sample.int(n, n1)))))
      ans$method <- "monte_carlo"
      ans$B <- B
    }
  }
  structure(ans, class = "sharp_null")
}

randomization_tail <- function(distribution, q) {
  if (!inherits(distribution, "sharp_null")) stop("Expected a sharp_null distribution.")
  if (!is.numeric(q) || anyNA(q)) stop("q must be numeric without missing values.")
  q <- as.double(q)
  d <- distribution
  if (d$method == "exact_wilcoxon") {
    if (d$statistic == "stephenson") q <- q * (d$n1 + d$n0 - 1) + d$n1
    u <- q - d$n1 * (d$n1 + 1.0) / 2
    return(stats::pwilcox(ceiling(.snap_integer(u)) - 1,
                         d$n1, d$n0, lower.tail = FALSE))
  }
  if (d$method == "exact_dp") {
    return(vapply(q, function(a) {
      if (a <= 0) return(1)
      if (a > 1) return(0)
      h <- ceiling(.snap_integer(a * d$n1 * d$n0))
      key <- sprintf("%.0f", h)
      if (!exists(key, d$cache, inherits = FALSE))
        assign(key, .rks_boundary_tail(h, d$n1, d$n0), d$cache)
      get(key, d$cache, inherits = FALSE)
    }, numeric(1)))
  }
  vapply(q, function(a) {
    if (a == -Inf) return(1)
    if (a == Inf) return(if (d$method == "monte_carlo") 1 / (d$B + 1) else 0)
    # Including near-ties is conservative and avoids platform rounding effects.
    cut <- a - 32 * .Machine$double.eps * max(1, abs(a))
    exceed <- length(d$values) - findInterval(cut, d$values, left.open = TRUE)
    if (d$method == "monte_carlo") (1 + exceed) / (d$B + 1)
    else exceed / length(d$values)
  }, numeric(1))
}

.get_null <- function(distribution, z, statistic, s, method, B, seed) {
  if (is.null(distribution))
    return(randomization_distribution(sum(z), length(z) - sum(z), statistic,
                                      s, method, B, seed))
  if (!inherits(distribution, "sharp_null") || distribution$n1 != sum(z) ||
      distribution$n0 != length(z) - sum(z) || distribution$statistic != statistic ||
      (statistic == "stephenson" && distribution$s != s))
    stop("distribution does not match the group sizes and statistic.")
  distribution
}
