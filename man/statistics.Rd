\name{rank_statistic}
\alias{rank_statistic}
\alias{rks_tail}
\alias{randomization_distribution}
\alias{randomization_tail}
\title{Rank Statistics and Their Randomization Distributions}
\description{Compute one-sided rank statistics and inclusive upper-tail probabilities under complete randomization.}
\usage{
rank_statistic(z, y, statistic = "rks", s = 2L, tie_order = NULL, seed = NULL)
rks_tail(q, n1, n0)
randomization_distribution(n1, n0, statistic = "rks", s = 2L,
  method = c("auto", "exact", "monte_carlo"), B = 9999L,
  seed = NULL, max_enumerate = 200000L)
randomization_tail(distribution, q)
}
\arguments{
\item{z}{Binary assignment vector containing both arms.}
\item{y}{Numeric outcomes of the same length; infinite values are allowed for rank calculations.}
\item{statistic}{One of \code{"rks"}, \code{"wilcoxon"}, or \code{"stephenson"}.}
\item{s}{Stephenson subset size, an integer from 2 to the total sample size.}
\item{tie_order}{Distinct numeric priorities, smaller first. If NULL, generate one uniform random permutation. Explicit priorities must be independent of assignment and outcomes and travel with units when rows are reordered.}
\item{q}{Numeric vector of statistic thresholds.}
\item{n1,n0}{Positive treated and control sample sizes.}
\item{method}{Automatic calibration, required exact calibration, or Monte Carlo calibration.}
\item{B}{Number of independent uniform random assignments for Monte Carlo calibration.}
\item{seed}{Optional seed. When supplied, the caller's RNG state is restored.}
\item{max_enumerate}{Maximum number of assignments to enumerate for Stephenson statistics.}
\item{distribution}{A reference object returned by \code{randomization_distribution}.}
}
\details{
The RKS statistic is the maximum of the control minus treated empirical CDFs
after unique ranking. Its distribution is evaluated by a backward lattice-path
boundary-hitting recursion, in \eqn{O(n_1 n_0)} work and \eqn{O(n_1)} storage per
threshold. Equal arm sizes use the equivalent reflection formula
\eqn{\binom{2m}{m+h}/\binom{2m}{m}}, where \eqn{h=\lceil mq\rceil}.
The returned probability is \eqn{P(T \ge q)}, including ties.
This uses the exact finite reference law with floating-point arithmetic.

Wilcoxon uses ranks 1 through n. Stephenson uses
\eqn{\binom{r-1}{s-1}/\binom{n-1}{s-1}}, a common positive rescaling of the
usual statistic that avoids overflow. For s=2 its law is computed using the
exact Wilcoxon distribution. Other Stephenson statistics use enumeration
when feasible, or Monte Carlo. Requesting infeasible exact enumeration raises
an error. Monte Carlo tails use the correction \eqn{(1+E)/(B+1)}, where E is
the number of reference statistics at least as large as the observed statistic.
Reference simulations are held fixed within an object. RKS tails are cached.
}
\value{\code{rank_statistic} returns a scalar with the realized \code{tie_order} and \code{seed} attributes; \code{rks_tail} and
\code{randomization_tail} return numeric probability vectors.
\code{randomization_distribution} returns a \code{sharp_null} object with
group sizes, statistic, calibration method, and reference draws or a DP cache.}
\examples{
rank_statistic(c(0,1,0,1), c(1,2,3,4))
rks_tail(c(0, 0.5, 1), n1 = 2, n0 = 2)
d <- randomization_distribution(4, 4, "stephenson", s = 3)
randomization_tail(d, c(0, 1))
}
