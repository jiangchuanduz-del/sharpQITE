# Full-sample comparison; run from the repository or reproducibility-bundle root.
local({
  library(sharpQITE)
  source('shared/ri-helpers.R')
  source('shared/kernel-helpers.R')
  base <- 'real-data/keller-2025'
  source(file.path(base,'helpers.R'))
  a <- readRDS(file.path(base,'results/analysis.rds'))
  d <- a$data; z <- 1-d$mathGrp; y <- d$vocabPost; n <- nrow(d)
  halves <- list(a$models$vocabulary$training_idx,a$models$vocabulary$inference_idx)
  stopifnot(n==2200L,identical(sort(unlist(halves)),seq_len(n)),
    length(intersect(halves[[1]],halves[[2]]))==0L,
    all(lengths(halves)==1100L),sum(z)==1136L,
    all(is.finite(y)),all(d$female %in% 0:1))
  out <- file.path(base,'full-data-results')
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  alpha <- .1; B <- 100000L; tol <- 1e-6
  pool <- getFromNamespace('.pool_lower','sharpQITE')
  # Fixed transformations require no training, so use complete-sample inference.
  key <- list(n1=sum(z),n0=sum(1-z),alpha=alpha/2,B=B,seed=22002025L)
  cache <- file.path(out,'calibration.rds')
  saved <- if(file.exists(cache))readRDS(cache) else NULL
  if(!is.null(saved) && identical(saved$key,key)) cal <- saved$calibration else {
    message('Calibrating full-sample rank statistics (100,000 assignments).')
    cal <- calibrate_rank_statistics(key$n1,key$n0,key$alpha,key$B,key$seed)
    saveRDS(list(key=key,calibration=cal),cache)
  }
  cals <- list(treated=cal,control=reverse_calibration(cal))
  arm_bands <- function(zz,yy,label,cc) list(
    treated=fast_arm_lower(zz,yy,label,cc$treated,tol),
    control=fast_arm_lower(1-zz,-yy,label,cc$control,tol))
  labels <- c('Wilcoxon','Stephenson (s=6)','RKS','RKS + naive adjusted',
    'RKS + optimal adjusted (vocab pretest)*',
    'RKS + optimal adjusted (vocab pretest + sex)*')
  bands <- parts <- models <- list()
  for(label in labels[1:4]) {
    message('Full-sample inference: ',label)
    yy <- if(label==labels[4])y-d$vocabPre else y
    statistic <- if(label==labels[4])'RKS' else label
    parts[[label]] <- arm_bands(z,yy,statistic,cals)
    bands[[label]] <- pool(parts[[label]])
  }
  sets <- list(vocab_pretest='vocabPre',vocab_pretest_sex=c('vocabPre','female'))
  fold_calibrations <- lapply(halves,function(te) {
    m <- sum(z[te]); n0 <- length(te)-m
    h <- rks_critical_unbalanced(m,n0,alpha/4)
    list(treated=list(n1=m,n0=n0,alpha=alpha/4,rks_h=h),
      control=list(n1=n0,n0=m,alpha=alpha/4,rks_h=h))
  })
  tail_probability <- getFromNamespace('.rks_boundary_tail','sharpQITE')
  for(cc in c(list(cals$treated),lapply(fold_calibrations,`[[`,'treated'))) {
    stopifnot(tail_probability(cc$rks_h,cc$n1,cc$n0)<=cc$alpha,
      tail_probability(cc$rks_h-1,cc$n1,cc$n0)>cc$alpha)
  }
  for(j in seq_along(sets)) {
    cols <- sets[[j]]; x <- as.matrix(d[,cols,drop=FALSE])
    storage.mode(x) <- 'double'
    label <- labels[4+j]; fold_bands <- fits <- vector('list',2)
    # First train on the original first half, then swap training/inference.
    for(step in 1:2) {
      tr <- halves[[step]]; f <- 3L-step; te <- halves[[f]]
      message(label,': train half ',step,', infer half ',f)
      tuning <- tune_keller_kernel(x[tr,,drop=FALSE],z[tr],y[tr],
        setNames(list(cols),names(sets)[j]),a$settings$kernel_config$multipliers,
        seed=7101L+step,folds=a$settings$kernel_config$folds)
      bw <- kernel_bandwidth(x[tr,,drop=FALSE],z[tr],tuning$selected$multiplier)
      fit <- fit_kernel_cdf(x[tr,,drop=FALSE],z[tr],y[tr],bandwidth=bw)
      adj <- estimate_adjustment_fast(fit,x[te,,drop=FALSE],beta=.5)
      stopifnot(all(is.finite(adj$adjustment)))
      if(j==1L && step==1L) stopifnot(
        identical(adj$adjustment,a$models$vocabulary$adjustment$adjustment))
      local <- arm_bands(z[te],y[te]-adj$adjustment,'RKS',fold_calibrations[[f]])
      fold_bands[[f]] <- pool(local)
      if(j==1L && step==1L) {
        old <- a$bands[a$bands$outcome=='vocabulary' &
          a$bands$method=='RKS adjusted (beta=0.5)',]
        stopifnot(all(fold_bands[[f]]$bound<=old$lower))
      }
      fits[[f]] <- list(training_idx=tr,inference_idx=te,feature_columns=cols,
        tuning=tuning,fit=fit,adjustment=adj,arm_bands=local)
    }
    bands[[label]] <- pool(fold_bands)
    parts[[label]] <- fold_bands; models[[label]] <- fits
    direct <- pool(unlist(lapply(fits,function(m)m$arm_bands),recursive=FALSE))
    stopifnot(identical(bands[[label]]$bound,direct$bound))
  }
  ks <- seq_len(10L)*220L
  lower <- sapply(bands,function(b)b$bound[ks])
  closed <- sapply(bands,function(b)b$closed[ks])
  for(b in bands) stopifnot(nrow(b)==n,all(diff(b$bound)>=0 | is.nan(diff(b$bound))))
  # Independent order-statistic pooling check at every displayed finite bound.
  for(label in labels[5:6]) for(k in ks) {
    v <- bands[[label]]$bound[k]
    if(is.finite(v)) stopifnot(sum(vapply(parts[[label]],function(b)sum(b$bound<v),integer(1)))<k,
      sum(vapply(parts[[label]],function(b)sum(b$bound<=v),integer(1)))>=k)
  }
  note <- paste('90% simultaneous lower confidence bounds for individual treatment-effect',
    'quantiles on vocabulary posttest scores, vocabulary versus math training, all 2,200',
    'participants (1,136 vocabulary, 1,064 math). The original seed-1 shuffled halves',
    'are retained. For both learned adjustments, train on the first half and infer on',
    'the second, then swap. Allocate alpha=0.05 per inference half and alpha=0.025',
    'per treatment arm within each half. Pool and sort the two half-sample bound',
    'sequences to obtain full-sample quantile bounds, as in the selected screenshot.',
    'Bonferroni validity does not require independence between the two inference results.',
    'Wilcoxon, Stephenson, RKS, and naive subtraction use all 2,200 observations directly,',
    'with alpha=0.05 per treatment arm. Naive subtraction uses vocabPost-vocabPre and',
    'does not require fitting. *Estimated optimal adjustments use the specified covariates',
    'and beta_0=0.5 in each training half, fixed across reported quantiles. Bandwidths',
    'are selected by five-fold training-only CRPS. Sex is female (0=non-female, 1=female),',
    'treated as a numeric kernel covariate. All methods target the same posttest ITEs.',
    'Coverage is simultaneous across quantiles within each method, not across methods.',
    'RKS uses exact DP calibration; rank-sum methods use 100,000 Monte Carlo assignments',
    'with plus-one correction. -Inf means uninformative; 100% is the maximum ITE.',
    'Values are rounded to two decimals; RDS retains exact endpoints and inclusion flags.')
  settings <- list(alpha=alpha,fold_alpha=alpha/2,fold_arm_alpha=alpha/4,
    full_sample_arm_alpha=alpha/2,B=B,beta_adjust=.5,halves=halves,tol=tol,
    calibration_seed=key$seed,cv_seeds=c(7102L,7103L),
    fold_combination='pooled order statistics')
  saveRDS(list(lower=lower,lower_closed=closed,k=ks,quantile=ks/n,bands=bands,
    parts=parts,models=models,
    fold_calibrations=fold_calibrations,settings=settings,
    note=note,row_id=d$row_id),file.path(out,'vocabulary-table.rds'))
  fmt <- function(v)ifelse(is.infinite(v),'-Inf',formatC(v,format='f',digits=2))
  display <- cbind(paste0(seq_len(10)*10,'%'),apply(lower,2,fmt))
  headings <- c('Quantile',labels)
  line <- function(v)paste0('| ',paste(v,collapse=' | '),' |')
  writeLines(c('# Keller et al. (2025): all 2,200 participants','',line(headings),
    line(rep('---:',length(headings))),apply(display,1,line),'',note),file.path(out,'vocabulary-table.md'))
  rows <- apply(display,1,function(v)paste0('<tr><td>',paste(v,collapse='</td><td>'),'</td></tr>'))
  writeLines(c('<!doctype html><html lang="en"><meta charset="utf-8">',
    '<title>Keller vocabulary: all 2200 participants</title><style>body{font:16px Georgia,serif;max-width:1400px;margin:40px auto;padding:20px;color:#172b3a}.table-wrap{overflow-x:auto}table{border-collapse:collapse;width:100%;min-width:950px;font-variant-numeric:tabular-nums}th,td{padding:12px;text-align:right;border-bottom:1px solid #ddd}th{background:#edf2f6;font-size:14px}p{line-height:1.55;font-size:14px}</style>',
    '<h1>Vocabulary training: all 2,200 participants</h1><p>90% simultaneous lower confidence bounds; two-way sample splitting with Bonferroni for learned adjustments</p>',
    paste0('<div class="table-wrap"><table><thead><tr><th>',paste(headings,collapse='</th><th>'),'</th></tr></thead><tbody>'),
    rows,'</tbody></table></div>',paste0('<p>',note,'</p></html>')),file.path(out,'vocabulary-table.html'))
  writeLines(c('PASS: two disjoint exhaustive halves; every participant contributes inference once.',
    'PASS: original first-half vocab-pretest fit reproduced exactly.',
    'PASS: exact RKS tails verify all critical heights and Bonferroni arm levels.',
    'PASS: stricter first-direction Bonferroni bounds do not exceed original bounds.',
    'PASS: pooled half-sample bounds equal pooled four-arm bounds.',

    'PASS: all six 2200-element bound sequences are monotone; 10 requested ranks exported.'),
    file.path(out,'validation.txt'))
  cat(paste(readLines(file.path(out,'vocabulary-table.md')),collapse='\n'),'\n')

})
