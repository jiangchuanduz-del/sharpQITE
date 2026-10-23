# Run from the repository/bundle root after the main Keller analysis.
# Reuse its split, fitted pretest adjustment, and frozen calibrations.
local({
  library(sharpQITE)
  source('shared/ri-helpers.R')
  source('shared/kernel-helpers.R')
  base <- 'real-data/keller-2025'
  source(file.path(base,'helpers.R'))
  a <- readRDS(file.path(base,'results/analysis.rds'))
  model <- a$models$vocabulary
  te <- model$inference_idx
  n <- length(te)
  stopifnot(n==1100L, a$settings$alpha==.1,
    identical(model$feature_columns,'vocabPre'),
    a$settings$beta_adjust==.5,
    is.null(a$settings$outcome_scale) || a$settings$outcome_scale=='posttest')
  z <- 1-a$data$mathGrp[te]
  yy <- a$data$vocabPost[te]-a$data$vocabPre[te]
  parts <- list(
    treated=fast_arm_lower(z,yy,'RKS',model$calibration$treated,a$settings$tol),
    control=fast_arm_lower(1-z,-yy,'RKS',model$calibration$control,a$settings$tol))
  naive <- getFromNamespace('.pool_lower','sharpQITE')(parts)
  ks <- as.integer(seq_len(10L)*(n/10L))
  methods <- c('Wilcoxon','Stephenson (s=6)','RKS','RKS adjusted (beta=0.5)')
  bounds <- closed <- matrix(NA,10,6)
  for(j in seq_along(methods)) {
    b <- a$bands[a$bands$outcome=='vocabulary' & a$bands$method==methods[j],]
    stopifnot(nrow(b)==n,identical(as.integer(b$k),seq_len(n)),all(diff(b$lower)>=0 | is.nan(diff(b$lower))))
    bounds[,j] <- b$lower[ks]
    closed[,j] <- b$lower_closed[ks]
  }
  bounds[,5] <- naive$bound[ks]
  closed[,5] <- naive$closed[ks]
  extra_sets <- list(vocab_pretest_sex=c('vocabPre','female'))
  extra_labels <- 'RKS adjusted (vocab pretest + sex)'
  colnames(bounds) <- colnames(closed) <- c(methods,'RKS adjusted (naive)',extra_labels)
  tr <- model$training_idx
  stopifnot(length(intersect(tr,te))==0L,all(a$data$female %in% 0:1))
  extra_models <- list()
  for(j in seq_along(extra_sets)) {
    key <- names(extra_sets)[j]; cols <- extra_sets[[j]]
    x <- as.matrix(a$data[,cols]); storage.mode(x) <- 'double'
    zz <- 1-a$data$mathGrp; y <- a$data$vocabPost
    message('Fitting vocabulary adjustment: ',paste(cols,collapse=' + '))
    # Identical folds and candidate multipliers to the own-pretest analysis.
    # Each covariate set is fixed separately; no selection using inference data.
    tuning <- tune_keller_kernel(x[tr,,drop=FALSE],zz[tr],y[tr],
      setNames(list(cols),key),a$settings$kernel_config$multipliers,
      seed=7102L,folds=a$settings$kernel_config$folds)
    stopifnot(identical(tuning$fold,model$tuning$fold))
    bw <- kernel_bandwidth(x[tr,,drop=FALSE],zz[tr],tuning$selected$multiplier)
    fit <- fit_kernel_cdf(x[tr,,drop=FALSE],zz[tr],y[tr],bandwidth=bw)
    adj <- estimate_adjustment_fast(fit,x[te,,drop=FALSE],beta=.5)
    adjusted_y <- y[te]-adj$adjustment
    parts <- list(
      treated=fast_arm_lower(z,adjusted_y,'RKS',model$calibration$treated,a$settings$tol),
      control=fast_arm_lower(1-z,-adjusted_y,'RKS',model$calibration$control,a$settings$tol))
    pooled <- getFromNamespace('.pool_lower','sharpQITE')(parts)
    stopifnot(all(is.finite(adj$adjustment)),
      all(diff(pooled$bound)>=0 | is.nan(diff(pooled$bound))))
    bounds[,5+j] <- pooled$bound[ks]; closed[,5+j] <- pooled$closed[ks]
    extra_models[[key]] <- list(feature_columns=cols,training_idx=tr,inference_idx=te,
      tuning=tuning,fit=fit,adjustment=adj,local_bands=parts,all_ranks=pooled)
    message('  Selected multiplier ',tuning$selected$multiplier,
      '; plug-in lower median ',adj$tau_lower)
  }
  previous_path <- file.path(base,'results/vocabulary-table.rds')
  if(file.exists(previous_path)) {
    previous <- readRDS(previous_path)
    common <- intersect(colnames(bounds),colnames(previous$lower))
    stopifnot(identical(bounds[,common],previous$lower[,common]),
      identical(closed[,common],previous$lower_closed[,common]))
  }
  # Cross-check against the independently saved earlier naive comparison,
  # when available in the development workspace.
  prior <- file.path(base,'naive-comparison/results/bands.csv')
  if(file.exists(prior)) {
    b <- read.csv(prior)
    b <- b[b$outcome=='vocabulary' & b$method=='RKS adjusted (naive)',]
    stopifnot(identical(naive$bound,b$lower),identical(naive$closed,b$lower_closed))
  }
  fmt <- function(v) ifelse(is.infinite(v),'-Inf',formatC(v,format='f',digits=2))
  column_order <- c('Wilcoxon','Stephenson (s=6)','RKS','RKS adjusted (naive)',
    'RKS adjusted (beta=0.5)','RKS adjusted (vocab pretest + sex)')
  bounds <- bounds[,column_order]; closed <- closed[,column_order]
  headings <- c('Quantile','Wilcoxon','Stephenson (s=6)','RKS',
    'RKS + naive adjusted','RKS + optimal adjusted (vocab pretest)*',
    'RKS + optimal adjusted (vocab pretest + sex)*')
  display <- cbind(paste0(10*seq_len(10L),'%'),apply(bounds,2,fmt))
  note <- paste('90% simultaneous lower confidence bounds for quantiles of individual',
    'treatment effects on vocabulary posttest score (score points). Vocabulary training',
    'versus math training; 1,100 held-out inference participants (564 vocabulary, 536 math),',
    'and 1,100 training participants. Quantiles use k = ceiling(n * q); 100% is the maximum ITE.',
    'Coverage is simultaneous across quantiles within each method, not jointly across methods.',
    '*Optimal adjusted means the package\'s estimated optimal adjustment using the indicated covariates,',
    'trained at beta_0 = 0.5 and held fixed across quantiles; it is not an oracle or',
    'separately optimized at every quantile. Naive adjusted uses vocabPost - vocabPre.',
    'Sex is the released female indicator (0 = non-female, 1 = female). Both',
    'optimal adjustments use Gaussian product-kernel conditional CDFs, with each',
    'covariate set fixed in advance and bandwidth selected using the same five',
    'training-only folds and CRPS criterion. The binary female indicator enters',
    'the kernel as a numeric 0/1 covariate, scaled by its training standard deviation.',
    'Subtracting the same pretreatment quantity from both potential outcomes preserves',
    'the posttest treatment effect. -Inf denotes an uninformative lower bound.',
    'Entries are rounded to two decimals; exact endpoints and open/closed flags are saved',
    'in vocabulary-table.rds. Calibration: exact RKS DP; 100,000 Monte Carlo assignments',
    'for Wilcoxon and Stephenson; alpha = 0.05 per arm, pooled to alpha = 0.10.')
  out <- file.path(base,'results')
  line <- function(v) paste0('| ',paste(v,collapse=' | '),' |')
  writeLines(c('# Keller et al. (2025): vocabulary training', '',
    line(headings),line(rep('---:',length(headings))),apply(display,1,line),'',note),
    file.path(out,'vocabulary-table.md'))
  rows <- apply(display,1,function(v)paste0('<tr><td>',paste(v,collapse='</td><td>'),'</td></tr>'))
  writeLines(c('<!doctype html><html lang="en"><meta charset="utf-8">',
    '<title>Keller vocabulary quantile bounds</title><style>body{font:16px Georgia,serif;max-width:1500px;margin:40px auto;padding:20px;color:#172b3a}.table-wrap{overflow-x:auto}table{border-collapse:collapse;width:100%;min-width:1000px;font-variant-numeric:tabular-nums}th,td{padding:12px;text-align:right;border-bottom:1px solid #ddd}th{background:#edf2f6;font-size:14px}p{line-height:1.55;font-size:14px}</style>',
    '<h1>Vocabulary training: posttest outcomes</h1><p>Keller et al. (2025) &middot; 90% simultaneous lower confidence bounds</p>',
    paste0('<div class="table-wrap"><table><thead><tr><th>',paste(headings,collapse='</th><th>'),'</th></tr></thead><tbody>'),
    rows,'</tbody></table></div>',paste0('<p>',note,'</p></html>')),
    file.path(out,'vocabulary-table.html'))
  saveRDS(list(quantile=seq_len(10)/10,k=ks,lower=bounds,lower_closed=closed,
    note=note,naive_all_ranks=naive,extra_models=extra_models,settings=a$settings),file.path(out,'vocabulary-table.rds'))
  writeLines(c('PASS: 10 quantiles, six methods, posttest outcome.',
    'PASS: additional fits use the same disjoint training sample, CV folds and frozen calibration.',
    'PASS: original five columns unchanged; all additional bounds monotone across ranks.',
    'PASS: existing bounds reused exactly; naive RKS recomputed with the same split and calibration.'),
    file.path(out,'vocabulary-table-validation.txt'))
  cat(paste(readLines(file.path(out,'vocabulary-table.md')),collapse='\n'),'\n')
})
