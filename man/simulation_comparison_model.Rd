\name{gaussian_covariate_model}
\alias{gaussian_covariate_model}
\alias{superpopulation_bounds}
\title{Covariate-Dependent Gaussian Models and Superpopulation Sharp Bounds}
\description{Specify conditional Gaussian potential outcomes and numerically evaluate the marginal and covariate-informed sharp lower quantile bounds in equations (24) and (35).}
\usage{
gaussian_covariate_model(mu0 = function(x) 3 * x,
  mu1 = function(x) 4 * x, sigma0 = function(x) exp(-0.25 * x),
  sigma1 = function(x) 3 * exp(0.25 * x), xlim = c(-1, 1))
superpopulation_bounds(model = gaussian_covariate_model(),
  betas = seq(.005, .995, by = .005), quadrature_n = 161L,
  outcome_grid_n = 8001L, fixed_beta = .5)
}
\arguments{
\item{mu0,mu1}{Vectorized conditional-mean functions.}
\item{sigma0,sigma1}{Vectorized strictly positive conditional-SD functions.}
\item{xlim}{Support of the uniform scalar covariate X.}
\item{model}{An object returned by \code{gaussian_covariate_model}.}
\item{betas}{Quantile levels strictly between zero and one.}
\item{quadrature_n}{Number of Gauss-Legendre quadrature nodes for integrating over X.}
\item{outcome_grid_n}{Number of outcome nodes for numerical inversion of marginal CDFs.}
\item{fixed_beta}{Quantile used to select a single oracle adjustment, held fixed across all reported quantiles.}
}
\details{
The model draws X uniformly on xlim and draws the two normal errors independently
conditional on X. Both marginal outcome distributions are generally Gaussian
mixtures. They are not obtained by substituting average means and SDs into a
single normal distribution. The induced potential outcomes need not be
unconditionally independent.

For (24), integrate the conditional normal CDFs to obtain F0 and F1, invert
them numerically, and maximize \eqn{F_1^{-1}(u)-F_0^{-1}(1+u-\beta)}.
All grid-local maxima are refined. For (35), the conditional density-equality
equation is quadratic in the outcome threshold; evaluate both roots and the
zero tail limit, then integrate the conditional minimum over X and invert.
These are deterministic numerical superpopulation bounds, not estimates from
simulation samples. Increasing quadrature_n and outcome_grid_n checks numerical
accuracy. There is no replication ribbon for these reference curves.

The additional \code{covariate_fixed} reference applies (24) to the distributions
of Y(1)-s(X) and Y(0)-s(X), where s(X) is the conditional CDF minimizer at the
covariate-informed sharp bound for \code{fixed_beta}. Its adjustment is held
fixed as the reported quantile changes. It agrees with (35) at fixed_beta, but
need not be sharp for the full covariate information elsewhere. If a conditional
minimum is attained only in the tails, a finite 11-SD tail proxy is used; the
\code{tail_proxy_nodes} attribute counts such nodes (zero in the default study).
}
\value{A \code{sharp_gaussian_model} list, or a data frame with beta,
\code{marginal} (24), \code{covariate} (35), and \code{covariate_fixed}
(sharp bound for the fixed-adjusted marginal distributions).}
\examples{
model <- gaussian_covariate_model()
superpopulation_bounds(model, betas = c(.1,.5,.9))
}
