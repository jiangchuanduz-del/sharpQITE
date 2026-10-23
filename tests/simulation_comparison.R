library(sharpQITE)
# Preserve the calibration stream when comparing with known reference draws.
fixed_band <- function(z,y,...,tie_order=seq_along(z))
  quantile_band(z,y,...,tie_order=tie_order)
near <- function(a,b,tol=1e-8) stopifnot(isTRUE(all.equal(a,b,tolerance=tol,check.attributes=FALSE)))

# Balanced exact reflection probability agrees with direct integer path counts.
for (m in c(2,5,12)) for (h in seq_len(m)) {
  counts <- matrix(0,m+1,m+1)
  for (i in 0:m) for (j in 0:m) {
    if (j-i >= h) next
    if (i==0 && j==0) counts[1,1] <- 1 else
      counts[i+1,j+1] <- (if(i) counts[i,j+1] else 0)+(if(j) counts[i+1,j] else 0)
  }
  near(rks_tail(h/m,m,m),1-counts[m+1,m+1]/choose(2*m,m))
}
stopifnot(is.finite(rks_tail(.01,50000,50000)))

# Fast continuous-outcome lower limits match the original exact breakpoint
# inversion, including tied data and finite-sample Monte Carlo calibration.
for (case in 1:3) {
  set.seed(case)
  z <- rep(0:1,6); y <- rnorm(12)+2*z
  if (case==3) y <- round(y)
  cal <- simulation_comparison_calibration(12,alpha=.2,B=999,seed=71)
  for (st in c("rks","wilcoxon","stephenson")) {
    old <- fixed_band(z,y,target="treated",alpha=.2,statistic=st,s=6,
                          method=if(st=="rks")"exact" else "monte_carlo",B=999,seed=71)
    fast <- treated_lower_bounds(z,y,(1:6)/6,statistic=st,s=6,alpha=.2,calibration=cal,tol=1e-5)
    stopifnot(identical(is.infinite(old$lower),is.infinite(fast$lower)))
    ok <- is.finite(old$lower)
    stopifnot(all(fast$lower[ok]<=old$lower[ok]+1e-9),all(old$lower[ok]-fast$lower[ok]<1.1e-5))
  }
}

# Known Gaussian sharp bound: both information sets coincide if X is irrelevant.
model <- gaussian_covariate_model(mu0=function(x)0*x,mu1=function(x)2+0*x,
                                  sigma0=function(x)1+0*x,sigma1=function(x)1+0*x)
bound <- superpopulation_bounds(model,c(.1,.5,.9),outcome_grid_n=4001)
expected <- 2+2*qnorm(bound$beta/2)
near(bound$marginal,expected,tol=1e-5)
near(bound$covariate,expected,tol=1e-6)
near(bound$covariate_fixed,expected,tol=1e-5)
bound2 <- superpopulation_bounds(betas=c(.1,.5,.9))
stopifnot(all(bound2$covariate >= bound2$marginal-1e-6))
stopifnot(all(bound2$covariate_fixed <= bound2$covariate+1e-4),
          attr(bound2,'tail_proxy_nodes')==0)
near(bound2$covariate_fixed[2],bound2$covariate[2],tol=1e-5)
stopifnot(any(bound2$covariate_fixed[c(1,3)] < bound2$covariate[c(1,3)]-.01))

# Grid values equal the independently computed kernel empirical CDF at knots.
set.seed(64)
x <- runif(30,-1,1); z <- rep(0:1,15); y <- 3*x+rnorm(30)*(1+2*z)
gridfit <- fit_kernel_cdf_grid(x,z,y,x_grid=seq(-1,1,length.out=21),outcome_grid_n=25,bandwidth=.3)
fit <- fit_kernel_cdf(x,z,y,bandwidth=.3)
near(gridfit$F0,predict_kernel_cdf(fit,gridfit$x_grid,gridfit$outcome_grid,arm=0))
near(gridfit$F1,predict_kernel_cdf(fit,gridfit$x_grid,gridfit$outcome_grid,arm=1))
near(sum(gridfit$integration_weights),1)
stopifnot(all(gridfit$integration_weights>=0))
a <- estimate_grid_adjustment(gridfit,x,c(.2,.5,.8))
stopifnot(all(is.finite(a$adjustment)),all(diff(a$tau_lower)>=0))
# Common outcome translation changes s_bar by that amount and leaves tau_hat.
shiftfit <- fit_kernel_cdf_grid(x,z,y+7,x_grid=gridfit$x_grid,outcome_grid_n=25,bandwidth=.3)
aa <- estimate_grid_adjustment(shiftfit,x,c(.2,.5,.8))
near(a$tau_lower,aa$tau_lower)
# Argmins need not be unique. Compare attained minima, not a particular
# selection from a flat interval of equally good adjustments.
for (i in c(1,11,21)) for (j in 1:3) {
  objective <- function(f,e,b) stats::approx(f$outcome_grid,f$F1[i,],xout=b,rule=2)$y -
    stats::approx(f$outcome_grid,f$F0[i,],xout=b-e,rule=2)$y
  value <- objective(gridfit,a$tau_lower[j],a$adjustment_grid[i,j])
  shifted <- objective(shiftfit,aa$tau_lower[j],aa$adjustment_grid[i,j])
  near(value,shifted)
  candidates <- c(gridfit$outcome_grid,gridfit$outcome_grid+a$tau_lower[j])
  near(value,min(objective(gridfit,a$tau_lower[j],candidates)))
}

# Large integer-valued Wilcoxon rank sums must not overflow at n=100000.
cal <- simulation_comparison_calibration(100000,alpha=.1,B=39,seed=22)
stopifnot(all(is.finite(cal$critical)))
b <- treated_lower_bounds(rep(0:1,50000),seq_len(100000),.9,"wilcoxon",alpha=.1,calibration=cal)
stopifnot(!anyNA(b$lower))

# End-to-end study: independent replication seeds and external samples.
sim <- simulate_simulation_comparison(sample_sizes=20,reps=2,betas=c(.2,.5,.8),B=99,progress=FALSE)
s <- summarize_simulation_comparison(sim)
stopifnot(nrow(sim$results)==42,nrow(s)==21,all(sim$diagnostics$external_n==20),
          length(unique(sim$diagnostics$seed))==2,sim$settings$adjusted_alpha==.05/3,
          sim$settings$fixed_adjustment_beta==.5,sim$settings$fixed_adjustment_alpha==.05)
stopifnot(all(sim$results$lower[sim$results$method=="RKS adjusted"] <=
               Inf))
# Reconstruct one replication and compare the full fixed-adjustment band with
# the original Corollary 2 inversion using the same external beta=.5 fit.
set.seed(sim$diagnostics$seed[1])
draw <- getFromNamespace('.draw_gaussian_population','sharpQITE')
pop <- draw(sim$model,20)
z <- integer(20); z[sample.int(20,10)] <- 1L
y <- ifelse(z==1,pop$y1,pop$y0)
train <- draw(sim$model,20)
zt <- integer(20); zt[sample.int(20,10)] <- 1L
yt <- ifelse(zt==1,train$y1,train$y0)
fit <- fit_kernel_cdf_grid(train$x,zt,yt,x_grid=seq(-1,1,length.out=41))
adj <- estimate_grid_adjustment(fit,pop$x,.5)
original <- fixed_band(z,y,target='treated',alpha=.05,adjustment=adj$adjustment[,1])
fixed <- sim$results[sim$results$replication==1 &
                      sim$results$method=='RKS adjusted (beta=0.5)',]
near(fixed$lower,original$lower[fixed$k])
# Infinite replicates are retained in the empirical quantiles, never discarded.
fake <- list(results=data.frame(n=100,method="RKS",beta=.2,lower=c(-Inf,-Inf,-1,0)))
q <- summarize_simulation_comparison(fake)
stopifnot(q$median == -Inf,q$finite_fraction==.5)
cat("Simulation comparison scaling, reference-bound, adjustment, and simulation checks passed.\n")
