\name{treated_lower_bounds}
\alias{simulation_comparison_calibration}
\alias{treated_lower_bounds}
\title{Scalable Treated-Group Lower Prediction Bounds for Balanced Experiments}
\description{Calculate selected lower prediction limits without materializing all cross-arm outcome differences.}
\usage{
simulation_comparison_calibration(n, alpha = .05, s = c(6L,10L), B = 9999L,
  seed = NULL)
treated_lower_bounds(z, y, betas = seq(.02,.98,by=.02),
  statistic = "rks", s = 6L, alpha = .05, calibration = NULL,
  tol = .002)
}
\arguments{
\item{n}{Even total sample size.}
\item{alpha}{One-sided treated-group error level. No across-arm correction is needed for this target.}
\item{s}{Stephenson subset sizes for calibration, or one size for inference.}
\item{B}{Monte Carlo reference draws for Wilcoxon and Stephenson statistics.}
\item{seed}{Optional local random seed.}
\item{z}{Binary assignment vector with equal arm sizes.}
\item{y}{Finite observed or independently adjusted outcomes.}
\item{betas}{Quantile levels in (0,1].}
\item{statistic}{\code{"rks"}, \code{"wilcoxon"}, or \code{"stephenson"}.}
\item{calibration}{Matching \code{simulation_comparison_calibration} for rank-sum inference; not needed for RKS.}
\item{tol}{Absolute tolerance for conservative rank-sum interval inversion.}
}
\details{
For equal group sizes m, the exact RKS tail is
\eqn{P(T\ge h/m)=\binom{2m}{m+h}/\binom{2m}{m}}. This reflection formula
equals the lattice-path DP. With h the smallest integer whose tail is at most
alpha, the lower limit at rank k is negative infinity for k<h, and otherwise
\eqn{\max_{1\le i\le k-h+1}\{Y^t_{(i)}-Y^c_{(i+m-k+h-1)}\}}.
It is the interval infimum; including that endpoint is conservative.

Rank-sum inference evaluates ranks using sorted control searches and numerically
inverts the test. When the requested precision is attainable, the returned
rejected-side endpoint lies at most tol below the boundary. If no distinct
floating-point midpoint exists before meeting tol, inversion stops, warns,
and returns the conservative lower bracket; its gap can exceed tol. Rescale
large-magnitude outcomes when finer absolute precision is needed. Reference critical values
implement the inclusive-tail plus-one test. One fixed set of reference draws
is shared across quantiles, preserving simultaneous validity with respect to
both randomization and reference simulation. All Wilcoxon reference calculations
here use Monte Carlo to avoid a very large exact rank-sum distribution.
No asymptotic KS or normal calibration is used.
}
\value{A calibration list containing scores, reference draws and critical
values, or a data frame with k, beta and lower. Infinite limits are retained.}
\examples{
z <- rep(0:1, 6)
y <- seq_len(12) + z
treated_lower_bounds(z, y, betas = c(.5,.8,.9))
cal <- simulation_comparison_calibration(12, s = c(6,10), B = 199, seed = 8)
treated_lower_bounds(z, y, betas = .8, statistic = "stephenson", s = 6,
  calibration = cal)
}
