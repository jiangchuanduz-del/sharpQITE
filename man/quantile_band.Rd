\name{quantile_band}
\alias{quantile_band}
\title{Simultaneous Prediction and Confidence Bands for ITE Quantiles}
\description{Invert bounded-null tests and pool arm-specific intervals according to Corollary 2.}
\usage{
quantile_band(z, y, alpha = 0.05, statistic = "rks", s = 2L,
  target = c("all", "treated", "control"),
  alternative = c("greater", "less", "two.sided"),
  adjustment = NULL, tie_order = NULL,
  method = c("auto", "exact", "monte_carlo"), B = 9999L, seed = NULL)
}
\arguments{
\item{z}{Binary assignment vector with both arms.}
\item{y}{Finite observed outcomes.}
\item{alpha}{Overall familywise error level in (0,1).}
\item{statistic}{\code{"rks"}, \code{"wilcoxon"}, or \code{"stephenson"}.}
\item{s}{Stephenson subset size.}
\item{target}{All-unit confidence band, or treated/control-unit prediction band.}
\item{alternative}{\code{"greater"} returns lower bounds; \code{"less"} upper bounds; \code{"two.sided"} both.}
\item{adjustment}{One fixed prespecified or independently trained adjustment per unit, used for the entire band. Positive infinity is allowed.}
\item{tie_order}{Distinct priorities independent of assignment and outcomes. If NULL, draw one uniform random permutation, retain it for every rank, arm and direction, and store it as an attribute.}
\item{method}{Reference calibration method.}
\item{B}{Monte Carlo reference size, when needed.}
\item{seed}{Optional local seed controlling random tie priorities followed by reference draws. The caller's RNG state is restored.}
}
\details{
Default random priorities protect against ties in assignment-sorted input.
When row order depends on assignment, validity averages over the independent
tie randomization as well as treatment assignment. A reproducibility seed is
not permission to choose favorable tie priorities. Explicit row-index priorities
are appropriate only for an independently ordered or randomly shuffled file.

For an all-unit one-sided band, alpha/2 is spent on each arm. For a two-sided
all-unit band, alpha/4 is spent on each arm and direction. There is no additional
correction over quantiles when the adjustment is held fixed. Bounds are found
by discrete inversion over outcome-difference breakpoints, with exact endpoint
inclusion. Monte Carlo inversion reuses the same reference draws for all
thresholds and ranks in an arm. Uninformative infinite bounds are retained.

For the all-unit band, the two lists of arm-specific lower intervals are sorted
by containment. At equal bounds, a closed lower interval comes before an open
one. Upper limits follow by reversing and negating lower limits for negative
outcomes. Plotting omits infinite limits; the returned object retains them.
}
\value{A \code{sharp_band} data frame with k, beta, lower, upper,
lower_closed, and upper_closed. Attributes record alpha, target, statistic,
alternative, per-arm error allocation, simultaneous validity, the realized tie_order, and seed.}
\references{Du, J. and Li, X. Sharp Design-based Inference for Quantiles of Individual Treatment Effects. Supplied manuscript, Corollary 2.}
\examples{
z <- rep(c(0,1), 5)
y <- c(0,1,1,3,2,4,3,6,4,8)
b <- quantile_band(z, y, alpha = .1, alternative = "two.sided")
as.data.frame(b)
plot(b)
}
