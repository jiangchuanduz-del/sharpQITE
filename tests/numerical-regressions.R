library(sharpQITE)
inside <- function(name)getFromNamespace(name,'sharpQITE')
near <- function(a,b,tol=1e-12)
  stopifnot(isTRUE(all.equal(a,b,tolerance=tol,check.attributes=FALSE)))

# At a singleton-arm upper tail, 1/20 must be represented as 0.05, not an
# accumulation slightly above it. Exercise raw tests, band endpoints and cuts.
for(sizes in list(c(1L,19L),c(19L,1L))) {
  m <- sizes[1]; n0 <- sizes[2]
  stopifnot(identical(rks_tail(1,m,n0),.05),
    inside('.rks_critical_height')(m,n0,.05)==19,
    is.infinite(inside('.rks_critical_height')(m,n0,.05-1e-10)))
  z <- c(rep(0,n0),rep(1,m)); y <- c(rep(0,n0),rep(1,m))
  test <- quantile_test(z,y,k=m,seed=11)
  band <- quantile_band(z,y,target='treated',seed=11)
  stopifnot(test$p.value==.05,test$p.value<=.05,band$lower[m]==1)
}
# Compare every attainable small-sample tail to independently enumerated paths.
for(m in 1:5) for(n0 in 1:5) {
  heights <- combn(m+n0,m,FUN=function(i) {
    z<-integer(m+n0);z[i]<-1L
    max(c(0,cumsum(1-z)*m-cumsum(z)*n0))
  })
  for(h in sort(unique(heights))) {
    p <- mean(heights>=h)
    stopifnot(inside('.rks_boundary_tail')(h,m,n0)==p)
    if(h>0 && p<1) stopifnot(inside('.rks_critical_height')(m,n0,p)<=h)
  }
}
for(alpha in .05+c(-1,0,1)*.Machine$double.eps*.05) {
  z <- c(rep(0,19),1); y <- z
  p <- quantile_test(z,y,k=1,seed=5)$p.value
  band <- quantile_band(z,y,target='treated',alpha=alpha,seed=5)
  rejects_zero <- band$lower[1]>0 || (band$lower[1]==0 && !band$lower_closed[1])
  stopifnot(identical(rejects_zero,as.logical(p<=alpha)))
  stopifnot(is.finite(inside('.rks_critical_height')(1,19,alpha))==(p<=alpha))
}

# Balanced cutoffs must use the identical represented law as public tails,
# both at and immediately beside attainable probabilities.
for(m in c(1:12,30,100,1000)) {
  p <- rks_tail(seq_len(m)/m,m,m)
  for(base in unique(p[p>0 & p<1])[seq_len(min(10,sum(p>0 & p<1)))]) {
    for(alpha in base+c(-1,0,1)*.Machine$double.eps*base) {
      if(alpha<=0 || alpha>=1) next
      reject <- which(p<=alpha)
      expected <- if(length(reject))reject[1] else Inf
      stopifnot(inside('.rks_critical_balanced')(m,alpha)==expected)
    }
  }
}

# A Monte Carlo cutoff must agree with (1+exceedances)/(B+1) <= alpha,
# even when alpha is a representable value just below a probability grid point.
for(B in c(19L,99L,9999L,100000L)) for(count in c(1,floor(.05*(B+1)),B)) {
  p <- count/(B+1)
  for(alpha in p+c(-1,0,1)*.Machine$double.eps*p) {
    if(alpha<=0 || alpha>=1) next
    j <- inside('.mc_critical_index')(B,alpha)
    allowed <- B-j+1
    stopifnot(allowed/(B+1)<=alpha,(allowed+1)/(B+1)>alpha)
    stopifnot((count<=allowed)==(p<=alpha))
  }
}

sim <- simulate_comparison(n=200,n1=100,nsim=1,betas=c(.07,.14,.28,.55,.56),
  training_n=8,stephenson_s=numeric(),B=19,
  scenarios=simulation_scenarios()['null'],seed=23)
stopifnot(identical(unique(sim$results$k),c(7L,14L,28L,55L,56L)) ||
          identical(unique(sim$results$k),c(7,14,28,55,56)))
stopifnot(all(sim$results$true_quantile==0))

x <- seq(-1,1,length.out=12); z <- rep(0:1,each=6)
fit <- fit_kernel_cdf_grid(x,z,seq_len(12))
for(newx in list(.2,c(-.5,.2,.6))) for(beta in list(.5,c(.25,.5,.75))) {
  a <- estimate_grid_adjustment(fit,newx,betas=beta)
  stopifnot(identical(dim(a$adjustment),c(length(newx),length(beta))))
  for(i in seq_along(newx)) near(a$adjustment[i,],
    as.numeric(estimate_grid_adjustment(fit,newx[i],betas=beta)$adjustment))
}

# The old loop hangs between adjacent doubles 1e16-2 and 1e16. Require both
# termination and a conservative bracket, with an explicit precision warning.
warned <- FALSE
bound <- withCallingHandlers(inside('.rank_lower_sorted')(
  1e16+c(0,4),c(0,2),2,1:4,5,1e-6),warning=function(w) {
    stopifnot(grepl('floating-point resolution',conditionMessage(w)))
    warned <<- TRUE; invokeRestart('muffleWarning')
  })
stopifnot(warned,is.finite(bound),bound<=1e16,1e16-bound<=4)
near(inside('.midpoint')(-1e308,1e308),0)

large <- fit_kernel_cdf_grid(x,z,c(rep(0,6),rep(1e16,6)))
warned <- FALSE
a <- withCallingHandlers(estimate_grid_adjustment(large,.2,betas=c(.25,.5,.75),tol=1e-8),
  warning=function(w) {
    stopifnot(grepl('floating-point resolution',conditionMessage(w)))
    warned <<- TRUE; invokeRestart('muffleWarning')
  })
stopifnot(warned,all(is.finite(a$adjustment)),identical(dim(a$adjustment),c(1L,3L)))
cat('Passed exact-tail equality, legacy rank rounding, singleton matrices and precision-limited inversion regressions.\n')
