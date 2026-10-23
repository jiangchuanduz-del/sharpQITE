# Run after installing sharpQITE. No additional packages are required.
library(sharpQITE)
set.seed(2026)
n <- 40L
z <- sample(rep(c(0L, 1L), each = n / 2L))
y0 <- rnorm(n)
tau <- 1 + rnorm(n, sd = 0.5)
y <- y0 + z * tau

# Test H0: the 90th ITE percentile among treated units is <= 0.
test <- quantile_test(z, y, beta = 0.9, c = 0,
  statistic = "rks", target = "treated", alternative = "greater", seed = 1)
print(test)

# 90% simultaneous lower confidence bounds for all-unit ITE quantiles.
band <- quantile_band(z, y, alpha = 0.1, statistic = "rks",
  target = "all", alternative = "greater", seed = 1)
print(as.data.frame(band)[c(10, 20, 30, 40), c("beta", "lower")])
if (interactive()) plot(band)

stopifnot(length(test$p.value) == 1L, test$p.value >= 0, test$p.value <= 1,
  nrow(band) == n, all(diff(band$lower) >= 0 | is.nan(diff(band$lower))))
