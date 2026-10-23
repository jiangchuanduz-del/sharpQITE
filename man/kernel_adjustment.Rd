\name{fit_kernel_cdf}
\alias{fit_kernel_cdf}
\alias{predict_kernel_cdf}
\alias{sharp_cdf_bound}
\alias{estimate_adjustment}
\title{Kernel Conditional CDFs and Quantile-Targeted Covariate Adjustment}
\description{Estimate equations (35)-(36) and the covariate adjustment in Theorem 11 using Gaussian-kernel weighted empirical distribution functions.}
\usage{
fit_kernel_cdf(x, z, y, bandwidth = NULL)
predict_kernel_cdf(fit, newx, thresholds, arm = 1L)
sharp_cdf_bound(fit, effect, newx = fit$x)
estimate_adjustment(fit, newx, beta, integration_x = fit$x,
  max_breaks = 1000000L)
}
\arguments{
\item{x}{Finite numeric covariate vector or matrix, one row per training unit. Encode categorical variables explicitly.}
\item{z}{Binary training assignment vector containing both arms.}
\item{y}{Finite numeric training outcomes.}
\item{bandwidth}{Positive scalar or vector with one bandwidth per column. Defaults to a standard-deviation rule with dimension-dependent shrinkage.}
\item{fit}{A \code{sharp_kernel_cdf} object.}
\item{newx}{Query covariates. For multiple covariates, use a matrix or a vector describing one query unit.}
\item{thresholds}{Numeric CDF evaluation thresholds; infinite values are allowed.}
\item{arm}{Potential-outcome arm, 0 or 1.}
\item{effect}{Finite effect thresholds for the plug-in sharp CDF bound.}
\item{beta}{Target quantile level in (0,1).}
\item{integration_x}{Covariate observations over which conditional bounds are averaged. Defaults to training covariates.}
\item{max_breaks}{Cap on training treated-control pairs used for discrete inversion.}
}
\details{
Within each arm, weights are proportional to
\eqn{\exp[-\frac12\sum_j\{(X_{ij}-x_j)/h_j\}^2]}.
The estimated CDF is their weighted mean of outcome threshold indicators.

The estimator computes
\eqn{\widehat{\overline F}(e)=1+\operatorname{average}_x
\inf_b\{\widehat F_1(b|x)-\widehat F_0(b-e|x)\}}.
Minimization searches all shifted control outcome jumps and positive infinity.
The generalized inverse at beta is found over training treated-control outcome
differences. The adjustment is an argmin evaluated at that estimated sharp
lower bound. Positive infinity is selected for a numerically zero minimum.
This is a plug-in bound, not a confidence bound.

Finite-sample inference requires training independent of inference assignments
and outcomes. Oracle asymptotic sharpness in Theorem 11 additionally requires
the theorem's distributional assumptions and appropriate consistency of the
estimated adjustment. The estimated function should be held fixed while
inverting thresholds for one quantile. Different adjustments at different
quantiles do not automatically produce a simultaneous band.
}
\value{\code{fit_kernel_cdf} returns training data and bandwidths in a
\code{sharp_kernel_cdf} object. \code{predict_kernel_cdf} returns a matrix,
with query units in rows and thresholds in columns. \code{sharp_cdf_bound}
returns one estimated upper CDF bound per effect threshold.
\code{estimate_adjustment} returns a list with \code{adjustment},
\code{tau_lower}, beta, estimated CDF at the bound, conditional minima,
bandwidth, and effective local sample sizes in each arm.}
\references{Du, J. and Li, X. Sharp Design-based Inference for Quantiles of Individual Treatment Effects. Supplied manuscript, equations (35)-(36) and Theorem 11.}
\examples{
x <- c(-1,0,1,-1,0,1)
z <- c(0,0,0,1,1,1)
y <- c(-2,0,2,-1,2,4)
fit <- fit_kernel_cdf(x, z, y, bandwidth = 1)
predict_kernel_cdf(fit, c(0,1), c(-Inf,0,2,Inf))
sharp_cdf_bound(fit, c(0,1,2))
estimate_adjustment(fit, newx = c(-.5,.5), beta = .8)
}
