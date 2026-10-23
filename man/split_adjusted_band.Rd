\name{split_adjusted_band}
\alias{split_adjusted_band}
\title{Sample-Split Covariate-Adjusted Simultaneous Confidence Bands}
\description{Use opposite folds to train fixed kernel adjustments, infer each fold's ITEs, and pool the results into an all-unit simultaneous RKS band.}
\usage{
split_adjusted_band(x, z, y, beta_adjust = 0.8, alpha = 0.05,
  alternative = c("greater", "less", "two.sided"),
  folds = NULL, bandwidth = NULL, seed = NULL, tie_order = NULL)
}
\arguments{
\item{x}{Numeric covariates, one row per unit.}
\item{z}{Binary assignment vector with both arms.}
\item{y}{Finite observed outcomes.}
\item{beta_adjust}{Quantile level used to train one fixed adjustment for each fold's entire band.}
\item{alpha}{Overall error level.}
\item{alternative}{Lower, upper, or two-sided bands, specified as in \code{quantile_band}.}
\item{folds}{Prespecified assignment-independent fold labels 1 and 2. If NULL, randomly split into two roughly equal folds.}
\item{tie_order}{Optional independent distinct unit priorities. If NULL, generate one random permutation after constructing folds; each fold uses its corresponding subset.}
\item{bandwidth}{Kernel bandwidth, passed to \code{fit_kernel_cdf}.}
\item{seed}{Optional seed for the random split followed by random tie priorities, preserving caller RNG state.}
}
\details{
The held-out assignment is completely randomized conditional on the training
assignments and its realized arm counts. Each fold receives error level alpha/2;
arm and direction corrections are then performed within that fold. Sorted
unit bounds from both folds are pooled. Dependence between the fold errors
does not invalidate the Bonferroni guarantee. This constructs simultaneous
full-population bounds using one fixed adjustment per fold.

Folds must not be chosen by examining outcomes or treatment assignments. If a
fold lacks an arm, uninformative intervals are returned for the affected
analysis; the routine never repeatedly splits until a favorable split appears.
The same adjustment targets beta_adjust in both directions, so upper bounds
are valid but not separately optimized for upper-tail power.
}
\value{A \code{sharp_band} data frame. Additional attributes contain folds,
fold-specific bands, fitted adjustments, beta_adjust, the realized tie_order, and seed.}
\examples{
set.seed(3)
x <- rnorm(20)
z <- sample(rep(0:1, each = 10))
y <- x + rnorm(20) + 2*z
b <- split_adjusted_band(x, z, y, seed = 41)
as.data.frame(b)[c(10,16,18,20), ]
}
