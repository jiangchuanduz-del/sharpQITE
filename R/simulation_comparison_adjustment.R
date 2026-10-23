.weighted_bins <- function(group, value, n) {
  out <- numeric(n)
  a <- rowsum(value, group, reorder = FALSE)
  out[as.integer(rownames(a))] <- a[,1]
  out
}

fit_kernel_cdf_grid <- function(x,z,y,x_grid=NULL,outcome_grid_n=NULL,bandwidth=NULL) {
  .check_data(z,y)
  if (!is.numeric(x) || length(x)!=length(y) || any(!is.finite(x)))
    stop("This scalable estimator requires one finite numeric covariate per unit.")
  n <- length(y)
  if (is.null(x_grid)) x_grid <- seq(min(x),max(x),length.out=max(41L,2L*ceiling(n^.3)+1L))
  if (!is.numeric(x_grid) || length(x_grid)<2 || any(!is.finite(x_grid)) ||
      any(diff(x_grid)<=0)) stop("x_grid must be finite and strictly increasing.")
  if (is.null(outcome_grid_n)) outcome_grid_n <- max(101L,2L*ceiling(sqrt(n))+1L)
  outcome_grid_n <- .integer_scalar(outcome_grid_n,"outcome_grid_n",10)
  if (is.null(bandwidth)) {
    sx <- stats::sd(x); if (!is.finite(sx) || sx==0) sx <- 1
    bandwidth <- sx*min(sum(z),sum(1-z))^(-1/5)
  }
  if (!is.numeric(bandwidth) || length(bandwidth)!=1 || !is.finite(bandwidth) || bandwidth<=0)
    stop("bandwidth must be positive.")
  knots <- unique(as.numeric(stats::quantile(y,seq(0,1,length.out=outcome_grid_n),type=1)))
  pad <- max(1,stats::sd(y))
  knots <- c(min(y)-pad,knots,max(y)+pad)
  fit_arm <- function(arm) {
    idx <- which(z==arm); idx <- idx[order(y[idx])]
    position <- findInterval(knots,y[idx]) + 1L
    cdfs <- matrix(0,length(x_grid),length(knots))
    ess <- numeric(length(x_grid))
    for (j in seq_along(x_grid)) {
      d <- ((x[idx]-x_grid[j])/bandwidth)^2
      if (any(!is.finite(d))) stop("Kernel distances overflow; rescale x.")
      w <- exp(-.5*(d-min(d))); w <- w/sum(w)
      cdfs[j,] <- c(0,cumsum(w))[position]
      ess[j] <- 1/sum(w^2)
    }
    cdfs[,1] <- 0; cdfs[,ncol(cdfs)] <- 1
    list(cdf=cdfs,ess=ess)
  }
  f0 <- fit_arm(0); f1 <- fit_arm(1)
  g <- pmin(length(x_grid)-1L,pmax(1L,findInterval(x,x_grid)))
  fraction <- pmin(1,pmax(0,(x-x_grid[g])/(x_grid[g+1L]-x_grid[g])))
  integration_weights <- (.weighted_bins(g,1-fraction,length(x_grid))+
                           .weighted_bins(g+1L,fraction,length(x_grid)))/n
  structure(list(x_grid=x_grid,outcome_grid=knots,F0=f0$cdf,F1=f1$cdf,
                 integration_weights=integration_weights,bandwidth=bandwidth,
                 ess0=f0$ess,ess1=f1$ess,n=n,
                 representation="kernel-weighted empirical CDF at knots; linear interpolation"),
            class="sharp_kernel_grid")
}

.interpolate_cdf_grid <- function(values,knots,q) {
  j <- pmax(1L,pmin(length(knots)-1L,findInterval(q,knots)))
  f <- pmax(0,pmin(1,(q-knots[j])/(knots[j+1L]-knots[j])))
  values[,j,drop=FALSE] * rep(1-f,each=nrow(values)) +
    values[,j+1L,drop=FALSE] * rep(f,each=nrow(values))
}

.kernel_grid_minimum <- function(fit,effect) {
  b <- fit$outcome_grid
  d1 <- fit$F1-.interpolate_cdf_grid(fit$F0,b,b-effect)
  d0 <- .interpolate_cdf_grid(fit$F1,b,b+effect)-fit$F0
  j1 <- max.col(-d1,ties.method="first"); j0 <- max.col(-d0,ties.method="first")
  i <- seq_len(nrow(d1))
  v1 <- d1[cbind(i,j1)]; v0 <- d0[cbind(i,j0)]
  use1 <- v1<=v0
  value <- pmin(0,v1,v0)
  adjustment <- ifelse(use1,b[j1],b[j0]+effect)
  # With the finite-support interpolants a zero minimum is attained at a
  # finite point below both supports. Avoid artificial +/-Inf outcome ties.
  adjustment[value >= -64*.Machine$double.eps] <- min(b)+min(0,effect)-1
  list(minimum=value,adjustment=adjustment,
       Fbar=1+sum(fit$integration_weights*value))
}

estimate_grid_adjustment <- function(fit,newx,betas=seq(.02,.98,by=.02),tol=.005) {
  if (!inherits(fit,"sharp_kernel_grid")) stop("Use fit_kernel_cdf_grid().")
  if (!is.numeric(newx) || any(!is.finite(newx))) stop("newx must be finite numeric values.")
  if (!is.numeric(betas) || !length(betas) || anyNA(betas) || any(betas<=0 | betas>=1))
    stop("betas must be in (0,1).")
  if (!is.numeric(tol) || length(tol)!=1 || !is.finite(tol) || tol<=0) stop("tol must be positive.")
  cache <- new.env(parent=emptyenv())
  evaluate <- function(e) {
    key <- sprintf("%.17g",e)
    if (!exists(key,cache,inherits=FALSE)) assign(key,.kernel_grid_minimum(fit,e),cache)
    get(key,cache,inherits=FALSE)
  }
  span <- diff(range(fit$outcome_grid))
  if (!is.finite(span)) stop('Rescale outcomes to avoid overflow in adjustment inversion.')
  lo <- rep(-span,length(betas)); hi <- rep(span,length(betas))
  stalled <- logical(length(betas))
  while (any(hi-lo>tol & !stalled)) {
    mids <- .midpoint(lo,hi)
    stalled <- stalled | ((mids<=lo | mids>=hi) & hi-lo>tol)
    move <- mids>lo & mids<hi & !stalled
    if (!any(move)) break
    unique_mid <- unique(mids)
    vals <- vapply(unique_mid,function(e) evaluate(e)$Fbar,numeric(1))[match(mids,unique_mid)]
    accept <- vals>=betas
    hi[accept & move] <- mids[accept & move]; lo[!accept & move] <- mids[!accept & move]
  }
  if (any(stalled)) warning('Requested adjustment-inversion tolerance is below floating-point resolution; returning the upper brackets.',call.=FALSE)
  at_grid <- vapply(hi,function(e) evaluate(e)$adjustment,numeric(length(fit$x_grid)))
  estimates <- matrix(vapply(seq_along(betas),function(j)
    stats::approx(fit$x_grid,at_grid[,j],xout=newx,rule=2)$y,numeric(length(newx))),
    nrow=length(newx),ncol=length(betas))
  list(adjustment=estimates,tau_lower=hi,betas=betas,x_grid=fit$x_grid,
       adjustment_grid=at_grid,bandwidth=fit$bandwidth,
       median_ess0=stats::median(fit$ess0),median_ess1=stats::median(fit$ess1),
       tolerance=tol,representation=fit$representation)
}
