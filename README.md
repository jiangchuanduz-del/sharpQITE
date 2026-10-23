# sharpQITE

sharpQITE provides design-based inference for quantiles of individual treatment effects in completely randomized experiments. It implements Rank-based Kolmogorov-Smirnov (RKS), Wilcoxon, and Stephenson randomization tests, simultaneous confidence bands, and covariate adjustments estimated using training samples. Examples reproduce simulation comparisons and a vocabulary-training application, including quantile inference tables and adjustment function plots by pretest score and sex.

## Installation

Requires **R >= 4.1.0**.

```r
devtools::install_github("jiangchuanduz-del/sharpQITE")
```

## Load the package

```r
library(sharpQITE)
```

## Main functions

| Function | Purpose |
|---|---|
| `quantile_test()` | Test a bounded null for an ITE quantile among treated or control units. |
| `quantile_band()` | Construct simultaneous lower, upper, or two-sided bands for ITE quantiles. |
| `rank_statistic()` | Compute the RKS, Wilcoxon, or Stephenson statistic. |
| `randomization_distribution()` | Obtain an exact or Monte Carlo null distribution. |
| `fit_kernel_cdf()` | Estimate conditional distributions for outcomes from training data. |
| `estimate_adjustment()` | Estimate an optimal adjustment for a specified quantile from fitted conditional distributions. |
| `split_adjusted_band()` | Construct a two-fold, covariate-adjusted RKS band with Bonferroni correction. |
| `simulate_simulation_comparison()` / `summarize_simulation_comparison()` | Run and summarize the Gaussian covariate simulation study. |
| `plot_simulation_comparison()` | Plot the simulation comparison with sharp bounds. |
| `simulation_comparison_calibration()` | Generate the reference calibration for the simulation rank statistics. |


Use the R help pages for arguments, targets, calibration details, and examples:

```r
?quantile_test
?quantile_band
?estimate_adjustment
?split_adjusted_band
```

## A simple example

Generate observed outcomes under complete randomization, test a quantile null, and construct a lower confidence band for all units:

```r
library(sharpQITE)
set.seed(42)
n <- 40L
z <- sample(rep(c(0L, 1L), each = n / 2L))
y0 <- rnorm(n)
tau <- 1 + rnorm(n, sd = 0.5)
y <- y0 + z * tau

# Test H0: the 90th ITE percentile among treated units is <= 0.
test <- quantile_test(z, y, beta = 0.9, c = 0,
  statistic = "rks", target = "treated", alternative = "greater", seed = 1)
print(test)

# 90% simultaneous lower confidence bounds for ITE quantiles.
band <- quantile_band(z, y, alpha = 0.1, statistic = "rks",
  target = "all", alternative = "greater", seed = 1)
print(as.data.frame(band)[c(10, 20, 30, 40), c("beta", "lower")])
plot(band)

```

## Reproduce the results


```sh
# Regenerate the four deliverables from saved simulations and training fits.
Rscript demo_code/run_all.R render

# Validate the saved results without regenerating figures.
Rscript demo_code/run_all.R validate

# Refit the vocabulary analysis from the original 2,200-person dataset.
Rscript demo_code/run_all.R real-data

# Rerun both simulations with 300 replications per sample size.
Rscript demo_code/run_all.R simulations 300

# Rerun both simulations and the vocabulary analysis.
Rscript demo_code/run_all.R all 300
```

## Repository structure

```text
sharpQITE/
├── DESCRIPTION, NAMESPACE, LICENSE
├── R/                 # statistical functions
├── man/               # R help pages
├── tests/             # package tests
├── data-raw/          # original WSCdata, documentation and provenance
├── demo_code/
│   ├── example.R      # small usage example
│   ├── run_all.R      # reproduction entry point
│   ├── simulation/   # two simulation studies and saved results
│   ├── real-data/    # vocabulary application
│   ├── shared/       # calibration and kernel helpers
│   ├── scripts/      # export and validation scripts
│   └── results/      # the four selected deliverables
├── sharpQITE.Rproj
└── README.md
```
