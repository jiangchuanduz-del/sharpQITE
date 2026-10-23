library(sharpQITE)
base <- 'real-data/keller-2025'
out <- getOption('sharpQITE.keller.output',file.path(base,'results'))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
source('shared/ri-helpers.R')
source('shared/kernel-helpers.R')
source(file.path(base,'helpers.R'))
env <- new.env(); load('../data-raw/WSCdata.rda',env)
raw <- as.data.frame(env$WSCdata)
stopifnot(nrow(raw)==2200L,ncol(raw)==33L,all(vapply(raw,function(x)all(is.finite(x)),logical(1))),
  all(raw$mathGrp %in% 0:1),all(raw$mathSel %in% 0:1))
raw$row_id <- seq_len(nrow(raw))
set.seed(1); d <- raw[sample.int(nrow(raw)),]; rownames(d) <- NULL
train <- seq_len(1100L); test <- 1100L+seq_len(1100L)
d$split <- rep(c('training','inference'),each=1100L)
d$tie_priority <- seq_len(nrow(d))
# Default to the outcome's own pretest; alternative runs can supply candidate
# sets drawn only from the documented pretreatment variables below.
config <- getOption('sharpQITE.keller.kernel_config',list(
  sets_by_outcome=list(vocabulary=list(own_pretest='vocabPre')),
  multipliers=c(.25,.5,.75,1,1.5,2,3,4,6,8),folds=5L,
  selection=paste('Own-outcome pretest fixed in advance; bandwidth minimizes',
    'five-fold conditional-CDF CRPS, averaging the two arm means equally')))
allowed_baseline <- c('female','white','black','asian','hisp','married','logAge','income',
  'collegeS','collegeM','collegeD','calc','logBooks','mathLike','big5O','big5C',
  'big5E','big5A','big5N','AMAS','logBDI','MCS','GSES','vocabPre','mathPre','mathSel')
baseline_names <- unique(unlist(if(is.null(config$sets_by_outcome))config$sets else
  config$sets_by_outcome,use.names=FALSE))
stopifnot(length(baseline_names)>0L,all(baseline_names %in% allowed_baseline))
x <- as.matrix(d[,baseline_names,drop=FALSE]); storage.mode(x) <- 'double'
specifications <- data.frame(outcome=c('math','vocabulary'),variable=c('mathPost','vocabPost'),
  treatment=c('mathGrp','1 - mathGrp'),title=c('Math training: math posttest','Vocabulary training: vocabulary posttest'))
outcome_scale <- match.arg(getOption('sharpQITE.keller.outcome_scale','posttest'),c('posttest','gain'))
if(outcome_scale=='gain') {
  d$mathGain <- d$mathPost-d$mathPre
  d$vocabGain <- d$vocabPost-d$vocabPre
  specifications$variable <- c('mathGain','vocabGain')
  specifications$title <- c('Math training: posttest - pretest','Vocabulary training: posttest - pretest')
}
labels <- c('RKS','RKS adjusted (beta=0.5)','Wilcoxon','Stephenson (s=6)')
alpha <- .1; local_alpha <- .05; B <- 100000L; tol <- 1e-6
saveRDS(config,file.path(out,'configuration.rds'))
write.csv(d[c('row_id','split','mathGrp','mathSel','tie_priority')],file.path(out,'split-membership.csv'),row.names=FALSE)
write.csv(data.frame(variable=names(raw),missing=vapply(raw,function(v)sum(is.na(v)),integer(1))),
  file.path(out,'missingness-audit.csv'),row.names=FALSE)
write.csv(as.data.frame(table(split=d$split,preference_math=d$mathSel,assigned_math=d$mathGrp)),
  file.path(out,'design-counts.csv'),row.names=FALSE)
start <- proc.time()[[3]]
# Both outcomes use the same held-out assignment, with opposite treatment
# orientation for vocabulary. Reuse complementary calibrations exactly.
math_cal <- calibrate_rank_statistics(sum(d$mathGrp[test]),sum(1-d$mathGrp[test]),local_alpha,B,seed=2025L)
vocab_cal <- reverse_calibration(math_cal)
pool <- getFromNamespace('.pool_lower','sharpQITE')
results <- models <- diagnostics <- summaries <- cv_tables <- cv_scores <- cv_ids <- bandwidths <- list()
for(i in 2L) {
  spec <- specifications[i,]; key <- spec$outcome
  y <- d[[spec$variable]]; z <- if(key=='math')d$mathGrp else 1-d$mathGrp
  message(key,': training n=',length(train),'; inference n=',length(test),'; tuning kernel.')
  sets <- if(is.null(config$sets_by_outcome))config$sets else config$sets_by_outcome[[key]]
  tuning <- tune_keller_kernel(x[train,,drop=FALSE],z[train],y[train],sets,
    config$multipliers,seed=7100L+i,folds=config$folds)
  best <- tuning$selected; cols <- sets[[best$set]]
  bw <- kernel_bandwidth(x[train,cols,drop=FALSE],z[train],best$multiplier)
  fit <- fit_kernel_cdf(x[train,cols,drop=FALSE],z[train],y[train],bandwidth=bw)
  adj <- estimate_adjustment_fast(fit,x[test,cols,drop=FALSE],beta=.5)
  cals <- if(key=='math')list(treated=math_cal,control=vocab_cal)else list(treated=vocab_cal,control=math_cal)
  message('  Selected ',best$set,', bandwidth multiplier ',best$multiplier,
    '; plug-in lower median=',format(adj$tau_lower,digits=6))
  week_bands <- local_bands <- list()
  for(label in labels) {
    yy <- y[test]
    if(label=='RKS adjusted (beta=0.5)') yy <- yy-adj$adjustment
    parts <- list(treated=fast_arm_lower(z[test],yy,label,cals$treated,tol),
      control=fast_arm_lower(1-z[test],-yy,label,cals$control,tol))
    pooled <- pool(parts); n <- length(test)
    b <- data.frame(outcome=key,target='all',method=label,k=seq_len(n),beta=seq_len(n)/n,
      lower=pooled$bound,lower_closed=pooled$closed,alpha=alpha,local_alpha=local_alpha)
    finite <- which(is.finite(b$lower))
    summaries[[paste(key,label)]] <- data.frame(outcome=key,method=label,N=n,
      first_finite_k=min(finite),finite_limits=length(finite),median=b$lower[ceiling(.5*n)],
      q75=b$lower[ceiling(.75*n)],q90=b$lower[ceiling(.9*n)],q95=b$lower[ceiling(.95*n)])
    week_bands[[label]] <- b; local_bands[[label]] <- parts
  }
  results[[key]] <- do.call(rbind,week_bands)
  models[[key]] <- list(training_idx=train,inference_idx=test,fit=fit,adjustment=adj,
    tuning=tuning,feature_columns=cols,calibration=cals,local_bands=local_bands)
  diagnostics[[key]] <- data.frame(outcome=key,training_n=length(train),training_treated=sum(z[train]),
    inference_n=length(test),inference_treated=sum(z[test]),selected_set=best$set,features=length(cols),
    bandwidth_multiplier=best$multiplier,balanced_cv_crps=best$balanced_cv_crps,beta_adjust=.5,
    tau_plugin=adj$tau_lower,median_ess_treated=median(adj$effective_n_treated),
    median_ess_control=median(adj$effective_n_control),infinite_argmins=adj$infinite_argmins)
  cv_tables[[key]] <- cbind(outcome=key,tuning$table)
  cv_scores[[key]] <- cbind(outcome=key,tuning$scores)
  cv_ids[[key]] <- data.frame(outcome=key,row_id=d$row_id[train],treatment=z[train],fold=tuning$fold)
  bandwidths[[key]] <- data.frame(outcome=key,feature=cols,training_sd=apply(fit$x,2,sd),
    multiplier=best$multiplier,bandwidth=bw)
  saveRDS(models[[key]],file.path(out,paste0(key,'-model.rds')))
  write.csv(results[[key]],file.path(out,paste0(key,'-bands.csv')),row.names=FALSE)
}
combine <- function(xx) {v <- do.call(rbind,xx); rownames(v) <- NULL; v}
bands <- combine(results); summaries <- combine(summaries); diagnostics <- combine(diagnostics)
cv_tables <- combine(cv_tables); cv_scores <- combine(cv_scores); cv_ids <- combine(cv_ids)
bandwidths <- combine(bandwidths)
settings <- list(alpha=alpha,local_alpha=local_alpha,B=B,tol=tol,split_seed=1L,calibration_seed=2025L,
  raw_n=nrow(d),training_n=length(train),inference_n=length(test),beta_adjust=.5,
  baseline_names=baseline_names,kernel_config=config,outcomes=specifications,labels=labels,
  outcome_scale=outcome_scale,
  target='ITE quantiles among all 1100 inference participants, separately for each subject-training contrast',
  design='Complete randomization conditional on held-out arm counts; assignment is mathGrp, not mathSel',
  calibration='Exact RKS DP; 100,000 Monte Carlo permutations with plus-one correction for rank-sum statistics',
  adjustment='Training-only Gaussian kernel conditional CDFs, fixed beta=.5 over all quantiles')
saveRDS(list(data=d,features=x,bands=bands,models=models,settings=settings,summary=summaries,
  diagnostics=diagnostics,cv_table=cv_tables,cv_scores=cv_scores,cv_membership=cv_ids),file.path(out,'analysis.rds'))
for(name in c('bands','summaries','diagnostics','cv_tables','cv_scores','cv_ids','bandwidths')) {
  filename <- switch(name,summaries='method-summary',diagnostics='adjustment-diagnostics',
    cv_tables='training-cv-scores',cv_scores='training-cv-fold-scores',cv_ids='training-cv-membership',
    bandwidths='selected-bandwidths',name)
  write.csv(get(name),file.path(out,paste0(filename,'.csv')),row.names=FALSE)
}
write.csv(d,file.path(out,'data-with-split.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'session-info.txt'))
writeLines(sprintf('Elapsed seconds: %.3f',proc.time()[[3]]-start),file.path(out,'timing.txt'))
print(diagnostics,row.names=FALSE); print(summaries,row.names=FALSE)

