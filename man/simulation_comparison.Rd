\name{simulate_simulation_comparison}
\alias{simulate_simulation_comparison}
\alias{summarize_simulation_comparison}
\alias{plot_simulation_comparison}
\title{Simulation Comparisons with External Covariate Adjustment}
\description{Simulate lower prediction limits over increasing sample sizes and plot median curves, interquartile ribbons, and superpopulation reference bounds, including an adjustment held fixed at beta=0.5.}
\usage{
simulate_simulation_comparison(sample_sizes = c(100L,1000L,10000L,100000L),
  reps = 100L, betas = seq(.02,.98,by=.02), alpha = .05,
  s = c(6L,10L), external_ratio = 1, B = 9999L,
  model = gaussian_covariate_model(), inversion_tol = .002,
  adjustment_tol = .005, seed = 20261004L,
  checkpoint_dir = NULL, progress = interactive())
summarize_simulation_comparison(simulation)
plot_simulation_comparison(summary, sharp_bounds,
  adjusted = c("fixed","simultaneous","pointwise"), alpha = .05,
  reps = 100L, xlim = NULL)
}
\arguments{
\item{sample_sizes}{Distinct even inference sample sizes.}
\item{reps}{Independent population, assignment and external-training replications per size.}
\item{betas}{Increasing plotted quantiles strictly in (0,1).}
\item{alpha}{Nominal one-sided treated-group error level.}
\item{s}{Stephenson subset sizes.}
\item{external_ratio}{External sample size relative to n, rounded up to an even size.}
\item{B}{Shared Monte Carlo reference draws per sample size for rank-sum tests.}
\item{model}{Covariate-dependent Gaussian model.}
\item{inversion_tol}{Conservative rank-sum inversion tolerance.}
\item{adjustment_tol}{Tolerance for kernel-grid estimated sharp bounds.}
\item{seed}{Reproducibility seed; caller RNG state is preserved.}
\item{checkpoint_dir}{Optional directory for per-size checkpoints saved every ten replications.}
\item{progress}{Whether to print progress messages.}
\item{simulation}{Result of \code{simulate_simulation_comparison}.}
\item{summary}{Result of \code{summarize_simulation_comparison}.}
\item{sharp_bounds}{Data frame returned by \code{superpopulation_bounds}.}
\item{adjusted}{Default \code{fixed} plots only the fixed beta=0.5 adjusted statistic alongside unadjusted methods. Optional \code{simultaneous} and \code{pointwise} retain access to earlier varying-adjustment comparisons. Both displayed references are the marginal and full covariate sharp bounds; the fixed-adjusted reference is not plotted.}
\item{xlim}{Optional common horizontal plotting limits.}
}
\details{
Each replication draws a new population and an independent external sample
from the same conditional Gaussian model. Both experiments have half their
units treated. Only observed external outcomes train the kernel CDF estimator;
neither true potential outcomes nor known model parameters are used to train it.
True potential outcomes are retained only for coverage diagnostics.

The fixed beta=0.5 adjusted curve trains one adjustment on the external sample,
then holds that same function fixed across all quantiles. It spends alpha with
no quantile-wise correction and has native simultaneous validity.

The Bonferroni adjusted curve spends alpha divided by the number of plotted quantiles
at each quantile-specific adjustment. A cumulative maximum of these lower limits
is a valid monotone improvement. This gives joint coverage at the plotted ranks;
a constant step extension to subsequent ranks preserves coverage. Straight
segments in the plot are only visual guides. A separate pointwise adjusted
curve spends alpha per quantile and is explicitly not a simultaneous band.
Unadjusted methods share a statistic and reference calibration across quantiles
and retain their native simultaneous validity.

The summary uses empirical inverse-CDF quantiles (R type 1), retaining negative
infinity. Ribbons are 25--75 percent replication ranges, not confidence intervals
for the median curve. Infinite ribbon edges are clipped to the plotting window;
all limits and finite fractions remain in the returned data. Reference sharp
bounds are deterministic and have no replication ribbons.
The fixed-adjustment reference is sharp for the adjusted marginal distributions;
the pointwise covariate reference is sharp given the full covariate information.
Reference rank-sum
draws are held fixed across replications at each sample size, so ribbons condition
on that calibration. The RNG seeds for every replication and calibration are
saved for reconstruction.
}
\value{A simulation list with replicate results, training diagnostics,
calibrations, model and settings; a summary data frame; or invisibly the plotted
summary. Plotting uses the active graphics device and supports one to four
sample sizes.}
\examples{
sim <- simulate_simulation_comparison(sample_sizes = 20, reps = 2,
  betas = c(.2,.5,.8), B = 99, progress = FALSE)
summarize_simulation_comparison(sim)
}
