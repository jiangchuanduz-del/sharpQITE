a <- readRDS('real-data/keller-2025/results/vocabulary-table.rds')
b <- readRDS('real-data/keller-2025/full-data-results/vocabulary-table.rds')
stopifnot(identical(dim(a$lower),c(10L,6L)),identical(dim(b$lower),c(10L,6L)),
  b$settings$fold_combination=='pooled order statistics',b$settings$fold_arm_alpha==.025,
  identical(as.numeric(b$lower[,5]),c(-12,-6,-4,-3,-2,-1,0,1,2,2)),
  identical(as.numeric(b$lower[,6]),c(-12,-6,-4,-3,-2,-1,0,1,2,2)))
for(path in c('simulation/results/simulation_comparison.rds','simulation/adjustment-stress/results/simulation.rds')) {
  s <- readRDS(path)
  d <- s$results[s$results$n %in% c(100,1000,100000) & s$results$method %in% c('RKS','RKS adjusted (beta=0.5)','Wilcoxon','s=6'),]
  expected <- expand.grid(n=c(100,1000,100000),
    method=c('RKS','RKS adjusted (beta=0.5)','Wilcoxon','s=6'),
    beta=s$settings$betas,stringsAsFactors=FALSE)
  key <- function(x) paste(sprintf('%d',x$n),x$method,sprintf('%.12g',x$beta),sep='|')
  groups <- split(d$replication,key(d))
  stopifnot(setequal(names(groups),key(expected)),!anyNA(d$lower),
    all(vapply(groups,function(ids)
      identical(sort(as.integer(ids)),seq_len(s$settings$reps)),logical(1))))
}
g <- readRDS('real-data/keller-2025/results/adjustment-quantiles/vocabulary-adjustment-quantiles.rds')
stopifnot(all(g$curves$pretest==floor(g$curves$pretest)),setequal(unique(g$curves$panel),c('pretest','sex0','sex1')),
  all(table(g$curves$panel,g$curves$beta)==21L))
for(stem in c('01-simulation-comparison','02-vocabulary-table','03-vocabulary-adjustment','04-population-adjustments'))
  stopifnot(all(file.exists(paste0('results/',stem,c('.pdf','.png','.svg')))))
cat('PASS: screenshot table, simulation counts, integer-score curves and four deliverables verified.\n')

