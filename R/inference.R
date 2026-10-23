# Compare cross-arm differences directly. This makes a breakpoint c=Yt-Yc
# an exact tie, avoiding subtract-then-compare cancellation during inversion.
.worst_evaluator <- function(z, y, k, statistic, s, tie_order) {
  ti <- which(z == 1); ci <- which(z == 0)
  m <- length(ti); n0 <- length(ci)
  ord_t <- order(y[ti], tie_order[ti])
  yt <- y[ti]
  if (k < m) yt[ord_t[seq.int(k + 1L, m)]] <- -Inf
  yc <- y[ci]
  rt <- integer(m); rt[order(yt, tie_order[ti])] <- seq_len(m)
  rc <- integer(n0); rc[order(yc, tie_order[ci])] <- seq_len(n0)
  finite <- outer(is.finite(yt), is.finite(yc), "&")
  differences <- outer(yt, yc, "-")
  tie_before <- outer(tie_order[ti], tie_order[ci], ">")
  constant <- outer(yt, yc, ">") | (outer(yt, yc, "==") & tie_before)
  scores <- if (statistic == "rks") NULL else .scores(m + n0, statistic, s)
  evaluate <- function(c, side = "at") {
    before <- constant
    if (side == "left") {
      before[finite] <- differences[finite] >= c
    } else if (side == "right") {
      before[finite] <- differences[finite] > c
    } else {
      before[finite] <- differences[finite] > c |
        (differences[finite] == c & tie_before[finite])
    }
    ranks_t <- rt + rowSums(before)
    if (statistic != "rks") return(sum(scores[ranks_t]))
    zz <- integer(m + n0); zz[ranks_t] <- 1L
    .stat_order(zz, "rks", s)
  }
  list(evaluate = evaluate,
       breaks = sort(unique(differences[finite & is.finite(differences)])))
}

quantile_test <- function(z, y, k = NULL, beta = NULL, c = 0,
                          statistic = "rks", s = 2L,
                          target = c("treated", "control"),
                          alternative = c("greater", "less", "two.sided"),
                          adjustment = NULL, tie_order = NULL,
                          distribution = NULL,
                          method = c("auto", "exact", "monte_carlo"),
                          B = 9999L, seed = NULL) {
  .check_data(z, y)
  call <- match.call()
  target <- match.arg(target); alternative <- match.arg(alternative)
  statistic <- .statistic_name(statistic); method <- match.arg(method)
  yy <- .adjust_y(y, adjustment)
  if (target == "control") { z <- 1 - z; yy <- -yy }
  m <- sum(z)
  if (is.null(k) == is.null(beta)) stop("Supply exactly one of k and beta.")
  if (!is.null(beta)) k <- ceiling(.snap_integer(m * .probability(beta, "beta", TRUE)))
  k <- .integer_scalar(k, "k", 1, m)
  if (!is.numeric(c) || !length(c) || any(!is.finite(c)))
    stop("c must be a nonempty vector of finite thresholds.")
  .with_seed(seed, {
  tie_order <- .tie_order(tie_order, length(z))
  # Consume tie randomness before reference draws, without resetting the seed.
  d <- .get_null(distribution, z, statistic, s, method, B, NULL)
  lower_test <- function(yy, kk, cc) {
    f <- .worst_evaluator(z, yy, kk, statistic, s, tie_order)$evaluate
    value <- vapply(cc, f, numeric(1))
    list(value = value, p = randomization_tail(d, value))
  }
  if (alternative != "less") a <- lower_test(yy, k, c)
  if (alternative != "greater") b <- lower_test(-yy, m + 1L - k, -c)
  if (alternative == "two.sided") {
    p <- pmin(1, 2 * pmin(a$p, b$p))
    value <- cbind(greater = a$value, less = b$value)
  } else {
    result <- if (alternative == "greater") a else b
    p <- result$p; value <- result$value
  }
  structure(list(p.value = p, statistic = value, k = k, beta = k / m,
                 c = c, target = target, alternative = alternative,
                 statistic_name = statistic, s = if (statistic == "stephenson") s else NULL,
                 method = d$method, B = d$B, adjusted = !is.null(adjustment),
                 tie_order = tie_order, seed = seed, call = call), class = "sharp_test")
  })
}

print.sharp_test <- function(x, ...) {
  cat("Quantile bounded-null test:", x$statistic_name, "(", x$method, ")\n")
  cat("Target:", x$target, "units; rank", x$k, "; alternative:", x$alternative, "\n")
  print(data.frame(threshold = x$c, p.value = x$p.value), row.names = FALSE)
  invisible(x)
}

.invert_one <- function(z, y, k, alpha, statistic, s, tie_order, d) {
  f <- .worst_evaluator(z, y, k, statistic, s, tie_order)
  accepted <- function(c, side) randomization_tail(d, f$evaluate(c, side)) > alpha
  breaks <- f$breaks
  if (!length(breaks)) {
    return(list(bound = if (accepted(0, "at")) -Inf else Inf, closed = FALSE))
  }
  if (accepted(breaks[1], "left")) return(list(bound = -Inf, closed = FALSE))
  if (!accepted(breaks[length(breaks)], "right"))
    return(list(bound = Inf, closed = FALSE))
  lo <- 1L; hi <- length(breaks)
  while (lo < hi) {
    mid <- floor((lo + hi) / 2)
    if (accepted(breaks[mid], "right")) hi <- mid else lo <- mid + 1L
  }
  list(bound = breaks[lo], closed = accepted(breaks[lo], "at"))
}

.group_lower <- function(z, y, alpha, statistic, s, tie_order, d) {
  ans <- lapply(seq_len(sum(z)), function(k)
    .invert_one(z, y, k, alpha, statistic, s, tie_order, d))
  data.frame(bound = vapply(ans, `[[`, numeric(1), "bound"),
             closed = vapply(ans, `[[`, logical(1), "closed"))
}

.pool_lower <- function(parts) {
  a <- do.call(rbind, parts)
  # At the same boundary [b,Inf) contains (b,Inf), hence closed sorts first.
  a[order(a$bound, !a$closed), , drop = FALSE]
}

quantile_band <- function(z, y, alpha = 0.05, statistic = "rks", s = 2L,
                          target = c("all", "treated", "control"),
                          alternative = c("greater", "less", "two.sided"),
                          adjustment = NULL, tie_order = NULL,
                          method = c("auto", "exact", "monte_carlo"),
                          B = 9999L, seed = NULL) {
  .check_data(z, y)
  alpha <- .probability(alpha, "alpha")
  target <- match.arg(target); alternative <- match.arg(alternative)
  statistic <- .statistic_name(statistic); method <- match.arg(method)
  yy <- .adjust_y(y, adjustment)
  group_count <- if (target == "all") 2 else 1
  side_count <- if (alternative == "two.sided") 2 else 1
  local_alpha <- alpha / group_count / side_count
  .with_seed(seed, {
    tie_order <- .tie_order(tie_order, length(z))
    groups <- if (target == "all") c("treated", "control") else target
    ds <- lapply(groups, function(g) {
      zz <- if (g == "treated") z else 1 - z
      randomization_distribution(sum(zz), length(zz) - sum(zz), statistic, s, method, B)
    })
    lower <- function(sign) {
      parts <- lapply(seq_along(groups), function(i) {
        is_t <- groups[i] == "treated"
        zz <- if (is_t) z else 1 - z
        yyy <- sign * if (is_t) yy else -yy
        .group_lower(zz, yyy, local_alpha, statistic, s, tie_order, ds[[i]])
      })
      .pool_lower(parts)
    }
    size <- if (target == "all") length(z) else if (target == "treated") sum(z) else sum(1 - z)
    l <- data.frame(bound = rep(-Inf, size), closed = FALSE)
    u <- data.frame(bound = rep(Inf, size), closed = FALSE)
    if (alternative != "less") l <- lower(1)
    if (alternative != "greater") {
      neg <- lower(-1)
      u <- data.frame(bound = -rev(neg$bound), closed = rev(neg$closed))
    }
    ans <- data.frame(k = seq_len(size), beta = seq_len(size) / size,
                      lower = l$bound, upper = u$bound,
                      lower_closed = l$closed, upper_closed = u$closed)
    .new_band(ans, alpha, target, statistic, alternative, local_alpha,
              vapply(ds, `[[`, character(1), "method"), !is.null(adjustment),
              tie_order = tie_order, seed = seed)
  })
}

.new_band <- function(ans, alpha, target, statistic, alternative, local_alpha,
                      methods, adjusted, tie_order = NULL, seed = NULL) {
  structure(ans, class = c("sharp_band", "data.frame"), alpha = alpha,
            target = target, statistic = statistic, alternative = alternative,
            local_alpha = local_alpha, methods = methods, adjusted = adjusted,
            simultaneous = TRUE, tie_order = tie_order, seed = seed)
}

print.sharp_band <- function(x, ...) {
  cat(sprintf("%.1f%% simultaneous %s band for %s-unit ITE quantiles\n",
              100 * (1 - attr(x, "alpha")), attr(x, "statistic"), attr(x, "target")))
  print.data.frame(x, ...)
  invisible(x)
}

plot.sharp_band <- function(x, ..., xlab = "Quantile level", ylab = "ITE bound") {
  finite <- c(x$lower[is.finite(x$lower)], x$upper[is.finite(x$upper)])
  ylim <- if (length(finite)) range(finite) else c(-1, 1)
  if (diff(ylim) == 0) ylim <- ylim + c(-0.5, 0.5)
  graphics::plot(x$beta, rep(NA_real_, nrow(x)), type = "n", ylim = ylim,
                 xlab = xlab, ylab = ylab, ...)
  graphics::lines(x$beta, replace(x$lower, !is.finite(x$lower), NA_real_), type = "s", col = "#176B87")
  graphics::lines(x$beta, replace(x$upper, !is.finite(x$upper), NA_real_), type = "s", col = "#B65C32")
  invisible(x)
}
