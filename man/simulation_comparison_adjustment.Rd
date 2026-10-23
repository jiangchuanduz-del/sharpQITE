\name{fit_kernel_cdf_grid}
\alias{fit_kernel_cdf_grid}
\alias{estimate_grid_adjustment}
\title{Scalable Kernel-CDF Grid Approximation for External Adjustment}
\description{Compress externally trained kernel-weighted CDFs onto outcome and covariate grids, then estimate quantile-specific adjustment functions.}
\usage{
fit_kernel_cdf_grid(x, z, y, x_grid = NULL, outcome_grid_n = NULL,
  bandwidth = NULL)
estimate_grid_adjustment(fit, newx, betas = seq(.02,.98,by=.02),
  tol = .005)
}
\arguments{
\item{x}{One finite scalar training covariate per unit.}
\item{z}{Binary external training assignments, containing both arms.}
\item{y}{Finite external training outcomes.}
\item{x_grid}{Increasing covariate nodes; default grid grows approximately as the training size to power 0.3.}
\item{outcome_grid_n}{Number of pooled empirical-quantile outcome nodes; grows approximately as the square root of training size.}
\item{bandwidth}{Positive Gaussian kernel bandwidth; default uses the sample-SD rule.}
\item{fit}{An object returned by \code{fit_kernel_cdf_grid}.}
\item{newx}{Covariates of inference units.}
\item{betas}{Target quantile levels in (0,1).}
\item{tol}{Numerical inversion tolerance for the estimated sharp bound.}
}
\details{
At grid nodes, the fitted CDFs equal the Gaussian-kernel weighted empirical
CDFs. Between outcome nodes they are linearly interpolated. The minimum of the
difference of these piecewise-linear CDFs is evaluated at the union of the
treated knots and shifted control knots. Conditional minima are integrated
using the external empirical covariate distribution projected onto x_grid.
The estimated sharp bound is inverted by cached bisection; adjustment values
are then linearly interpolated across covariate nodes.
If floating-point spacing prevents attaining tol, bisection terminates with a
warning and returns the upper brackets. Rescale outcomes for finer precision.

This is an explicitly approximate, scalable version of the uncompressed
\code{fit_kernel_cdf} workflow, not a parametric Gaussian fit. Both grids grow
with training size. Since the interpolants have finite support, a zero minimum
can be attained at a finite threshold below both supports; this avoids creating
infinite adjusted outcomes. The outcome and covariate grids and effective sample
sizes are available for diagnostics. Grid refinement can change the selected
argmin, especially when there are several minima.

Training observations must be independent of inference observations and
assignments. Numerical or estimation error affects power but does not invalidate
randomization inference conditional on the fixed external adjustment. Different
columns of the adjustment matrix correspond to different beta-specific
functions and require multiplicity correction for joint coverage.
}
\value{The fit contains grids, arm-specific CDF matrices and integration
weights. The estimate contains an adjustment matrix (inference units by beta,
including a one-row matrix for a single inference unit),
estimated sharp bounds, adjustment values on x_grid, and diagnostics.}
\examples{
set.seed(4)
x <- runif(30,-1,1)
z <- rep(0:1,15)
y <- 3*x + rnorm(30)*(1+2*z)
fit <- fit_kernel_cdf_grid(x,z,y)
estimate_grid_adjustment(fit, newx = c(-.5,.5), betas = c(.5,.8))
}
