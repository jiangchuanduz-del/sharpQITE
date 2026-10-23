.rks_critical_balanced <- function(m, alpha) {
  i <- seq_len(m)
  hit <- which(exp(cumsum(log1p(-(2*i-1)/(m+i)))) <= alpha)
  h <- if (length(hit)) hit[1] else m+1
  # Refine the fast candidate with the same tail returned by rks_tail().
  # Small laws use integer counts; sum() and cumsum() can also differ by an ulp.
  while (h<=m && .rks_boundary_tail(h*m,m,m)>alpha) h <- h+1
  while (h>1 && .rks_boundary_tail((h-1)*m,m,m)<=alpha) h <- h-1
  if (h<=m) h else Inf
}

simulation_comparison_calibration <- function(n, alpha = .05, s = c(6L,10L), B = 9999L, seed = NULL) {
  n <- .integer_scalar(n, "n", 4)
  if (n %% 2) stop("Simulation comparison calibration requires even n and equal arm sizes.")
  alpha <- .probability(alpha, "alpha")
  B <- .integer_scalar(B, "B")
  s <- vapply(s, .integer_scalar, integer(1), name="s", lower=2, upper=n)
  if (!length(s) || anyDuplicated(s)) stop("s must be a nonempty vector of distinct subset sizes.")
  j <- .mc_critical_index(B,alpha)
  if (j > B) stop("B is too small to reject at this alpha.")
  scores <- c(list(Wilcoxon = as.numeric(seq_len(n))),
              stats::setNames(lapply(s, function(a) .scores(n,"stephenson",a)), paste0("s=",s)))
  draws <- .with_seed(seed, {
    out <- matrix(NA_real_,B,length(scores))
    for (b in seq_len(B)) {
      idx <- sample.int(n,n/2)
      out[b,] <- vapply(scores,function(v) sum(v[idx]),numeric(1))
    }
    out
  })
  crit <- apply(draws,2,function(v) sort(v,partial=j)[j])
  names(crit) <- names(scores)
  structure(list(n=n,alpha=alpha,B=B,s=s,scores=scores,critical=crit,
                 reference=draws,seed=seed,rks_h=.rks_critical_balanced(n/2,alpha)),
            class="sharp_simulation_comparison_calibration")
}

.rks_lower_sorted <- function(yt,yc,k,h) {
  m <- length(yt)
  if (!is.finite(h) || k < h) return(-Inf)
  i <- seq_len(k-h+1L)
  max(yt[i] - yc[i+m-k+h-1L])
}

.rank_lower_sorted <- function(yt,yc,k,score,critical,tol) {
  m <- length(yt); n <- m+length(yc); removed <- m-k
  const <- if (removed) sum(score[seq_len(removed)]) else 0
  ranks_max <- seq.int(n-k+1L,n)
  cushion <- 64*.Machine$double.eps*max(1,abs(critical))
  if (const+sum(score[ranks_max]) <= critical+cushion) return(-Inf)
  keep <- yt[seq_len(k)]
  bases <- removed+seq_len(k)
  rejects <- function(c) {
    count <- findInterval(keep-c,yc,left.open=TRUE)
    const+sum(score[bases+count]) > critical+cushion
  }
  lo <- keep[1]-yc[length(yc)]-1
  hi <- keep[k]-yc[1]+1
  if (!is.finite(lo) || !is.finite(hi)) stop("Rescale outcomes to avoid overflow.")
  while (hi-lo > tol) {
    mid <- .midpoint(lo,hi)
    if (mid<=lo || mid>=hi) {
      warning('Requested rank-inversion tolerance is below floating-point resolution; returning the conservative lower bracket.',call.=FALSE)
      break
    }
    if (rejects(mid)) lo <- mid else hi <- mid
  }
  # The rejected side is no larger than the exact interval infimum; including
  # this endpoint gives a conservative numerical lower bound. If adjacent
  # floating-point endpoints prevent attaining tol, the warning reports it.
  lo
}

treated_lower_bounds <- function(z,y,betas=seq(.02,.98,by=.02),statistic="rks",
                                  s=6L,alpha=.05,calibration=NULL,tol=.002) {
  .check_data(z,y)
  m <- sum(z)
  if (2*m != length(z)) stop("This scalable routine requires equal treatment/control sizes.")
  if (!is.numeric(betas) || !length(betas) || anyNA(betas) || any(betas<=0 | betas>1))
    stop("betas must be in (0,1].")
  alpha <- .probability(alpha,"alpha")
  if (!is.numeric(tol) || length(tol)!=1 || !is.finite(tol) || tol<=0) stop("tol must be positive.")
  statistic <- .statistic_name(statistic)
  k <- ceiling(.snap_integer(m*betas))
  yt <- sort(y[z==1]); yc <- sort(y[z==0])
  if (statistic=="rks") {
    h <- .rks_critical_balanced(m,alpha)
    bounds <- vapply(k,function(kk) .rks_lower_sorted(yt,yc,kk,h),numeric(1))
  } else {
    if (!inherits(calibration,"sharp_simulation_comparison_calibration") ||
        calibration$n != length(z) || calibration$alpha != alpha)
      stop("Provide a matching simulation_comparison_calibration for rank-sum inference.")
    name <- if (statistic=="wilcoxon") "Wilcoxon" else paste0("s=",s)
    if (!name %in% names(calibration$scores)) stop("This s is absent from calibration.")
    bounds <- vapply(k,function(kk) .rank_lower_sorted(yt,yc,kk,
      calibration$scores[[name]],calibration$critical[[name]],tol),numeric(1))
  }
  data.frame(k=k,beta=k/m,lower=bounds)
}
