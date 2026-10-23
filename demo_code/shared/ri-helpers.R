# Shared calibration and inversion for the retained real-data analyses.
rks_critical_unbalanced <- function(m,n0,alpha) {
  getFromNamespace('.rks_critical_height','sharpQITE')(m,n0,alpha)
}
rks_lower_unbalanced <- function(yt,yc,k,h) {
  if(!is.finite(h)) return(-Inf)
  m <- length(yt); n0 <- length(yc)
  i <- seq_len(k)
  j <- ceiling((h+(m-k+i-1)*n0)/m)
  keep <- j>=1 & j<=n0
  if(!any(keep)) return(-Inf)
  max(yt[i[keep]]-yc[j[keep]])
}
calibrate_rank_statistics <- function(m,n0,alpha=.05,B=100000L,seed=1L) {
  n <- m+n0
  score <- list(Wilcoxon=as.numeric(seq_len(n)),
    'Stephenson (s=6)'=getFromNamespace('.scores','sharpQITE')(n,'stephenson',6L),
    'Stephenson (s=10)'=getFromNamespace('.scores','sharpQITE')(n,'stephenson',10L))
  set.seed(seed)
  draws <- matrix(NA_real_,B,3L,dimnames=list(NULL,names(score)))
  for(b in seq_len(B)) {
    idx <- sample.int(n,m)
    draws[b,] <- vapply(score,function(s)sum(s[idx]),numeric(1))
  }
  # The plus-one Monte Carlo tail rejects only strictly above this cutoff.
  j <- getFromNamespace('.mc_critical_index','sharpQITE')(B,alpha)
  critical <- if(j>B)setNames(rep(Inf,length(score)),names(score)) else
    apply(draws,2,function(x)sort(x,partial=j)[j])
  list(n1=m,n0=n0,alpha=alpha,B=B,seed=seed,scores=score,
    critical=critical,draws=draws,monte_carlo_correction='plus_one',critical_index=j,rks_h=rks_critical_unbalanced(m,n0,alpha))
}
reverse_calibration <- function(cal) {
  reverse <- cal
  reverse$n1 <- cal$n0; reverse$n0 <- cal$n1
  reverse$draws <- sweep(-cal$draws,2,vapply(cal$scores,sum,numeric(1)),'+')
  j <- getFromNamespace('.mc_critical_index','sharpQITE')(cal$B,cal$alpha)
  reverse$critical <- if(j>cal$B)setNames(rep(Inf,length(cal$scores)),names(cal$scores)) else
    apply(reverse$draws,2,function(x)sort(x,partial=j)[j])
  reverse$critical_index <- j
  reverse$monte_carlo_correction <- 'plus_one' 
  # Reversing rank order and swapping arms preserves the one-sided KS law.
  reverse
}
fast_arm_lower <- function(z,y,label,cal,tol=1e-6) {
  yt <- sort(y[z==1]); yc <- sort(y[z==0]); m <- length(yt)
  if(label %in% c('RKS','RKS adjusted (beta=0.5)')) {
    lower <- vapply(seq_len(m),function(k)rks_lower_unbalanced(yt,yc,k,cal$rks_h),numeric(1))
  } else {
    inverse <- getFromNamespace('.rank_lower_sorted','sharpQITE')
    lower <- vapply(seq_len(m),function(k)
      inverse(yt,yc,k,cal$scores[[label]],cal$critical[[label]],tol),numeric(1))
  }
  # Closed intervals are conservative at discontinuities. Rank-sum endpoints
  # are rounded down within tol; no claim of exact open/closed inversion here.
  data.frame(bound=lower,closed=is.finite(lower))
}
