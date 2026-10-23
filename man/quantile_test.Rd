\name{quantile_test}
\alias{quantile_test}
\title{Randomization Tests for Quantile Bounded Null Hypotheses}
\description{Test a bound on an order statistic of individual treatment effects among the treated or control units, following Theorem 1 and its rank-sum analogues.}
\usage{
quantile_test(z, y, k = NULL, beta = NULL, c = 0,
  statistic = "rks", s = 2L, target = c("treated", "control"),
  alternative = c("greater", "less", "two.sided"),
  adjustment = NULL, tie_order = NULL, distribution = NULL,
  method = c("auto", "exact", "monte_carlo"), B = 9999L, seed = NULL)
}
\arguments{
\item{z}{Assignment vector with both 0 and 1.}
\item{y}{Finite observed outcomes, one per unit.}
\item{k}{Rank within the target arm, starting at 1. Supply exactly one of k and beta.}
\item{beta}{Quantile level in (0,1]; k is the ceiling of beta times target-arm size, after snapping floating-point values within machine tolerance of an integer.}
\item{c}{One or more finite null thresholds.}
\item{statistic}{\code{"rks"}, \code{"wilcoxon"}, or \code{"stephenson"}.}
\item{s}{Stephenson subset size.}
\item{target}{Arm whose ITE quantile is tested.}
\item{alternative}{Direction of the alternative; \code{"greater"} tests a null of ITE quantile at most c.}
\item{adjustment}{Optional prespecified or independently trained covariate adjustment, one value per unit. Finite values or positive infinity are allowed.}
\item{tie_order}{Distinct priorities for ties. If NULL, generate one uniform random permutation and save it in the result. Explicit priorities must be independent of assignment and outcomes and aligned with the units.}
\item{distribution}{Optional matching reference distribution from \code{randomization_distribution}. For control targets its treated size must equal the original control size.}
\item{method}{Reference calibration method; ignored if distribution is supplied.}
\item{B}{Monte Carlo reference size.}
\item{seed}{Optional seed for random tie priorities followed by reference calibration, preserving the caller's RNG state. Supply a saved tie_order to reuse priorities across calls.}
}
\details{
Default random priorities protect against ties in assignment-sorted input.
When row order depends on assignment, validity averages over the independent
tie randomization as well as treatment assignment. A reproducibility seed is
not permission to choose favorable tie priorities. Explicit row-index priorities
are appropriate only for an independently ordered or randomly shuffled file.

For the greater alternative, the largest m-k treated adjusted outcomes are
replaced by negative infinity, and c is subtracted from the remaining treated
outcomes. The returned p-value is the inclusive reference tail at this minimized
statistic. Control-arm inference switches labels and negates outcomes. The less
alternative negates outcomes and maps k to m+1-k. Two-sided testing uses twice
the smaller one-sided p-value, capped at one.

The target-arm null depends on the random treatment assignment. Validity means
\eqn{P(p \le \alpha, H_0\ \mathrm{true})\le\alpha}, not conditional type-I
error given that this random null is true. Theorem 1 is for treated-group
ranks, not ranks among all experimental units.

Training an adjustment on the same inference observations is not justified.
Use external training or \code{split_adjusted_band}. For estimated adjustments
that vary by quantile, separate tests have pointwise validity; their collection
does not automatically have simultaneous coverage.
}
\value{A \code{sharp_test} list containing \code{p.value}, the minimized
\code{statistic}, k, beta, threshold c, target, alternative, and calibration
metadata, the realized tie_order, and seed. Both p-value and statistic are vectors when c is a vector.}
\references{Du, J. and Li, X. Sharp Design-based Inference for Quantiles of Individual Treatment Effects. Supplied manuscript, Theorem 1 and Section 3.}
\examples{
z <- rep(c(0,1), 4)
y <- c(1,2,2,4,3,6,4,8)
quantile_test(z, y, beta = .75, c = c(0,1))
quantile_test(z, y, k = 3, statistic = "stephenson", s = 3)
}
