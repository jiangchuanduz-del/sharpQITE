\name{simulate_comparison}
\alias{simulation_scenarios}
\alias{simulate_comparison}
\alias{summarize_simulation}
\title{Compare Randomization Tests in Fixed Finite Populations}
\description{Generate populations once and repeatedly randomize treatment, comparing unadjusted and adjusted RKS, Wilcoxon, and several Stephenson tests.}
\usage{
simulation_scenarios()
simulate_comparison(n = 80L, n1 = floor(n / 2), nsim = 500L,
  betas = c(0.5, 0.8, 0.9), c = 0,
  stephenson_s = c(2L, 5L, 10L, 20L), training_n = 240L,
  alpha = 0.05, B = 19999L, scenarios = simulation_scenarios(),
  seed = 20261003L)
summarize_simulation(simulation, conf_level = 0.95)
}
\arguments{
\item{n}{Inference population size.}
\item{n1}{Number assigned to treatment in each trial.}
\item{nsim}{Number of assignments per scenario.}
\item{betas}{Treated-group quantile levels in (0,1).}
\item{c}{One finite null threshold.}
\item{stephenson_s}{Vector of subset sizes; must not exceed n.}
\item{training_n}{Independent external training sample size for adjusted RKS.}
\item{alpha}{Test rejection level.}
\item{B}{Monte Carlo reference size for non-enumerated Stephenson tests.}
\item{scenarios}{Named list of functions taking n and returning a list with x, y0, and y1.}
\item{seed}{Local seed.}
\item{simulation}{A simulation result or its replicate-level data frame.}
\item{conf_level}{Confidence level for Wilson intervals on rejection probabilities.}
}
\details{
The default scenarios have X standard normal and Y(0)=2X+epsilon with independent
standard normal epsilon. Effects are zero, constant 2.5, equal to 6 only when
X exceeds its population 75th percentile, or 1+1.5X. Populations and independent
training samples are drawn once per scenario and held fixed. Reference
calibrations are also shared across trials. Adjusted RKS uses additional training
data and a separate externally trained adjustment for each quantile. This is
an extra-data comparison, not an equal-total-sample comparison.

The null is about a realized treated-group order statistic and can vary with
assignment. The summary therefore includes the joint rate of rejection and
true null, the quantity controlled by the theorem. Conditional power given a
false null is descriptive. MCSE and Wilson intervals concern assignment
variability for the fixed finite population and fixed reference calibration;
they exclude population generation and reference-calibration uncertainty.
Runtime excludes model training and initial reference calibration.
}
\value{\code{simulation_scenarios} returns a named list of generators.
\code{simulate_comparison} returns a \code{sharp_simulation} list with replicate
results, fixed populations, training data, adjustments, and settings.
\code{summarize_simulation} returns a data frame grouped by scenario, method,
and beta, with rejection rates, MCSE, Wilson intervals, null probability,
joint false rejection, conditional power, and mean test runtime.}
\examples{
sim <- simulate_comparison(n = 12, nsim = 2, betas = .8,
  stephenson_s = c(2,3), training_n = 16, B = 99,
  scenarios = simulation_scenarios()[1], seed = 4)
summarize_simulation(sim)
}
