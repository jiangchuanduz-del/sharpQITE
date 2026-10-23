# Keller et al. (2025) source data

`WSCdata.rda` is the public dataset used for the vocabulary-training application.
The original 2,200 rows and all source columns are retained. The analysis uses
randomized assignment (`mathGrp`), vocabulary posttest (`vocabPost`), vocabulary
pretest (`vocabPre`) and the female indicator (`female`).

Upstream: https://github.com/jzangela/WSCdata

The exact source commit, file hashes, documentation and original licenses are
in `sources/`. `demo_code/real-data/keller-2025/run.R` reads this file directly;
no download, imputation or exclusion of participants is needed. This raw-data
folder is excluded from the installable statistical package and retained in
the full repository/reproducibility ZIP.
