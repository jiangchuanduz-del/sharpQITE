library(sharpQITE)
near <- function(a,b,tol=1e-10)
  stopifnot(isTRUE(all.equal(a,b,tolerance=tol,check.attributes=FALSE)))

# Fractional quantiles at integer boundaries must select exactly their ranks,
# in both arms. A genuine increment above the boundary must still round up.
for(m in c(50L,100L)) {
  z <- rep(0:1,each=m); y <- seq_along(z)
  for(target in c("treated","control")) for(k in seq_len(m)) {
    a <- quantile_test(z,y,beta=k/m,target=target,seed=31)
    b <- quantile_test(z,y,k=k,target=target,tie_order=a$tie_order)
    stopifnot(a$k==k,identical(a$p.value,b$p.value))
  }
  stopifnot(quantile_test(z,y,beta=.07+1e-9,seed=31)$k==ceiling(m*(.07+1e-9)))
}

# The seed controls one stored random tie permutation and subsequent reference
# draws, and preserves the caller's RNG. Explicit priorities bypass tie draws.
z <- rep(c(1,0),each=6); y <- rep(c(0,1,2),4)
set.seed(2026); before <- .Random.seed
a <- quantile_test(z,y,k=6,target="control",seed=18)
b <- quantile_test(z,y,k=6,target="control",seed=18)
stopifnot(identical(before,.Random.seed),identical(a,b),
          identical(sort(a$tie_order),seq_along(z)),!identical(a$tie_order,seq_along(z)))
near(a$p.value,quantile_test(z,y,k=6,target="control",tie_order=a$tie_order)$p.value)
stat <- rank_statistic(z,y,seed=18)
stopifnot(identical(attr(stat,"tie_order"),a$tie_order))

for(target in c("treated","control","all")) {
  b <- quantile_band(z,y,alpha=.1,target=target,seed=18)
  stopifnot(identical(before,.Random.seed),identical(attr(b,"tie_order"),a$tie_order))
  if(target!="all") for(k in seq_len(sum(if(target=="treated")z else 1-z))) {
    p <- quantile_test(z,y,k=k,c=c(-1,0,1),target=target,seed=18)$p.value
    included <- c(-1,0,1)>b$lower[k] |
      (c(-1,0,1)==b$lower[k] & b$lower_closed[k])
    stopifnot(identical(included,p>.1))
  }
}
# A single randomization is retained for all ranks and both sides of a band.
mc <- quantile_band(z,y,alpha=.2,target="treated",statistic="stephenson",s=3,
                    method="monte_carlo",B=99,seed=18)
for(k in 1:6) {
  p <- quantile_test(z,y,k=k,c=c(-1,0,1),statistic="stephenson",s=3,
                     method="monte_carlo",B=99,seed=18)$p.value
  stopifnot(identical(p>.2,c(-1,0,1)>mc$lower[k] |
    (c(-1,0,1)==mc$lower[k] & mc$lower_closed[k])))
}
stopifnot(identical(before,.Random.seed))

# Explicit unit priorities travel with the units when rows are sorted by arm.
# This proves row-reordering invariance, including control and two-sided cases.
priority <- c(7,1,10,2,5,12,3,8,4,11,6,9)
perm <- order(-z,y)
for(target in c("treated","control","all")) for(alt in c("greater","less","two.sided")) {
  a <- quantile_band(z,y,alpha=.1,target=target,alternative=alt,tie_order=priority)
  b <- quantile_band(z[perm],y[perm],alpha=.1,target=target,alternative=alt,
                     tie_order=priority[perm])
  near(a,b)
}
# Fold inference uses the corresponding subset of the stored global priority,
# never raw row numbers. Include tied outcomes and treatment-sorted rows.
folds <- rep(1:2,6)
b <- split_adjusted_band(seq_along(z),z,y,folds=folds,seed=918)
stopifnot(identical(before,.Random.seed),
  identical(b,split_adjusted_band(seq_along(z),z,y,folds=folds,seed=918)))
for(f in 1:2) stopifnot(identical(attr(attr(b,"fold_bands")[[f]],"tie_order"),
                                attr(b,"tie_order")[folds==f]))

# Independent tiny Monte Carlo example: the smallest achievable p is 1/(B+1).
d <- randomization_distribution(3,3,"stephenson",3,"monte_carlo",B=9,seed=8)
stopifnot(randomization_tail(d,max(d$values)+1)==.1)
cat("Passed quantile-boundary, random-tie, reproducibility, row-order, split and Monte Carlo regressions.\n")
