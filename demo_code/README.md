# Reproduction examples

Start with [the main README](../README.md#reproduce-the-results).

- `example.R`: a small RKS quantile test and simultaneous band example.
- `run_all.R`: renders or recomputes the four selected results.
- `results/index.html`: local results viewer.

Run from the repository root with `Rscript demo_code/run_all.R render`.
Use `validate` for checks, `real-data` to refit the vocabulary application,
and `simulations 300` to repeat both simulation studies. Panel B retains
pooled Bonferroni fold bounds. The original data are in `../data-raw/`.

The vocabulary adjustment figure contains three panels: pretest only, and
pretest plus sex evaluated for males and females. All panels share axes and
a legend for beta = 0.1, 0.3, 0.5, 0.7, 0.9; the table uses beta = 0.5.
