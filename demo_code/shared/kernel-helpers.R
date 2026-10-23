# Scalable versions of the existing conditional-CDF and band calculations.
# These helpers are checked against sharpQITE's exact small-sample routines.
kernel_pair_fast <- function(fit,newx) {
  newx <- getFromNamespace('.as_x','sharpQITE')(newx,ncol(fit$x))
  weights <- getFromNamespace('.kernel_weights','sharpQITE')
  a <- weights(fit,newx,1L); b <- weights(fit,newx,0L)
  # Aggregate ties, then take prefix sums instead of multiplying by a dense
  # triangular indicator matrix. This is the same weighted empirical CDF.
  collapse <- function(part) {
    yy <- unique(part$y)
    group <- match(part$y,yy)
    w <- t(rowsum(t(part$weights),group,reorder=FALSE))
    cdf <- t(matrix(apply(w,1,cumsum),nrow=ncol(w),ncol=nrow(w)))
    list(y=yy,cdf=cdf,effective_n=part$effective_n)
  }
  list(treated=collapse(a),control=collapse(b))
}
conditional_min_fast <- function(pair,effect) {
  yt <- pair$treated$y; yc <- pair$control$y
  # Direct differences preserve exact comparisons at training breakpoints.
  cut <- vapply(yc,function(y0)sum(yt-y0<=effect),integer(1))
  f1 <- pair$treated$cdf[,pmax(1L,cut),drop=FALSE]
  f1[,cut==0L] <- 0
  delta <- f1-pair$control$cdf
  pick <- max.col(-delta,ties.method='first')
  value <- delta[cbind(seq_len(nrow(delta)),pick)]
  at_infinity <- value>=-64*.Machine$double.eps
  value <- pmax(-1,pmin(0,value)); value[at_infinity] <- 0
  adjustment <- yc[pick]+effect
  adjustment[at_infinity] <- Inf
  list(minimum=value,adjustment=adjustment)
}
estimate_adjustment_fast <- function(fit,newx,beta=.5) {
  pair <- kernel_pair_fast(fit,fit$x)
  grid <- sort(unique(as.vector(outer(pair$treated$y,pair$control$y,'-'))))
  fbar <- function(e)1+mean(conditional_min_fast(pair,e)$minimum)
  lo <- 1L; hi <- length(grid)
  while(lo<hi) {
    mid <- floor((lo+hi)/2)
    if(fbar(grid[mid])>=beta) hi <- mid else lo <- mid+1L
  }
  tau <- grid[lo]
  pred_pair <- kernel_pair_fast(fit,newx)
  pred <- conditional_min_fast(pred_pair,tau)
  fallback <- max(fit$y[fit$z==1],fit$y[fit$z==0]+tau)+1
  n_infinite <- sum(!is.finite(pred$adjustment))
  # Above both training supports the CDF difference is exactly zero. This
  # finite minimizer represents the same objective as a +Inf argmin, while
  # keeping the scalable rank-statistic inversion on finite residual outcomes.
  pred$adjustment[!is.finite(pred$adjustment)] <- fallback
  list(adjustment=pred$adjustment,tau_lower=tau,beta=beta,
    estimated_cdf=fbar(tau),conditional_minimum=pred$minimum,
    effective_n_treated=pred_pair$treated$effective_n,
    effective_n_control=pred_pair$control$effective_n,
    bandwidth=fit$bandwidth,infinite_argmins=n_infinite,finite_fallback=fallback)
}

kernel_bandwidth <- function(x,z,multiplier) {
  scale <- apply(x,2,sd)
  scale[!is.finite(scale)|scale==0] <- 1
  multiplier*scale*min(sum(z),sum(1-z))^(-1/(ncol(x)+4))
}

kernel_distance <- function(train,newx,bandwidth) {
  distance <- matrix(0,nrow(newx),nrow(train))
  for(j in seq_len(ncol(train)))
    distance <- distance+(outer(newx[,j],train[,j],'-')/bandwidth[j])^2
  distance
}

kernel_normalize <- function(distance,multiplier=1) {
  logw <- -distance/(2*multiplier^2)
  w <- exp(logw-apply(logw,1,max))
  w/rowSums(w)
}

weighted_crps <- function(weights,support,observed) {
  stopifnot(nrow(weights)==length(observed),ncol(weights)==length(support),
    all(diff(support)>=0))
  # CRPS = E|V-y| - .5 E|V-V'|, with the pair term evaluated by
  # integrating F(t)(1-F(t)) over the sorted empirical support.
  absolute_error <- rowSums(weights*abs(outer(observed,support,'-')))
  # A point mass has no pairwise spread. Preserve one score per query row.
  if (length(support)==1L) return(absolute_error)
  cumulative <- t(apply(weights,1,cumsum))
  cc <- cumulative[,-ncol(cumulative),drop=FALSE]
  absolute_error -
    as.vector((cc*(1-cc)) %*% diff(support))
}

