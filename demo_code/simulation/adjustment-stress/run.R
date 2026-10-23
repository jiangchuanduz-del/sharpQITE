# Run from the repository root: Rscript simulation/adjustment-stress/run.R [reps]
# Source this checkout to avoid silently using another installed package version.
source_dir <- "../R"
for (f in list.files(source_dir,pattern="[.]R$",full.names=TRUE)) source(f)
args <- commandArgs(trailingOnly=TRUE)
reps <- if(length(args)) as.integer(args[1]) else 300L
stopifnot(length(reps)==1L,!is.na(reps),reps>=2)
outdir <- "simulation/adjustment-stress/results"
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
start <- proc.time()[[3]]
model <- gaussian_covariate_model(mu0=function(x)0*x,mu1=function(x)6*x,
  sigma0=function(x)1+0*x,sigma1=function(x)exp(x))
settings <- list(sample_sizes=c(100L,1000L,10000L,100000L),reps=reps,
  betas=seq(.02,.98,.02),alpha=.05,adjustment_beta=.5,external_ratio=1,
  adjustment_tol=.005,seed=202610051L,power_beta=.9,power_threshold=1.5,
  description="X ~ U(-1,1), Y0 = epsilon0, Y1 = 6X + exp(X) epsilon1; independent standard normal errors")
# This DGP and the illustrative q=.9 test threshold were fixed after exploratory
# screening, before these fresh replication seeds were generated.
sharp <- superpopulation_bounds(model,settings$betas,quadrature_n=321,
                                 outcome_grid_n=16001,fixed_beta=.5)
stopifnot(attr(sharp,"tail_proxy_nodes")==0L)
oracle_effect <- attr(sharp,"fixed_effect")
write.csv(sharp,file.path(outdir,"superpopulation-bounds.csv"),row.names=FALSE)
set.seed(settings$seed)
seeds <- matrix(sample.int(.Machine$integer.max,length(settings$sample_sizes)*reps),
                nrow=reps,ncol=length(settings$sample_sizes))
allrows <- alldiag <- list()
for (nn in seq_along(settings$sample_sizes)) {
  n <- settings$sample_sizes[nn]; rows <- diagrows <- vector("list",reps)
  for(r in seq_len(reps)) {
    set.seed(seeds[r,nn]); t0 <- proc.time()[[3]]
    pop <- .draw_gaussian_population(model,n)
    z <- integer(n); z[sample.int(n,n/2)] <- 1L
    y <- ifelse(z==1,pop$y1,pop$y0)
    train <- .draw_gaussian_population(model,n)
    zt <- integer(n); zt[sample.int(n,n/2)] <- 1L
    yt <- ifelse(zt==1,train$y1,train$y0)
    fit <- fit_kernel_cdf_grid(train$x,zt,yt,
      x_grid=seq(-1,1,length.out=max(41L,2L*ceiling(n^.3)+1L)))
    adj <- estimate_grid_adjustment(fit,pop$x,betas=.5,tol=settings$adjustment_tol)
    oracle <- .normal_cdf_minimum(oracle_effect,.model_parameters(model,pop$x))$adjustment
    stopifnot(all(is.finite(oracle)))
    limits <- list(RKS=treated_lower_bounds(z,y,settings$betas,alpha=settings$alpha)$lower,
      `RKS adjusted (beta=0.5)`=treated_lower_bounds(z,y-adj$adjustment[,1],
           settings$betas,alpha=settings$alpha)$lower,
      `RKS oracle (beta=0.5)`=treated_lower_bounds(z,y-oracle,
           settings$betas,alpha=settings$alpha)$lower)
    k <- ceiling(.snap_integer(n/2*settings$betas))
    truth <- sort((pop$y1-pop$y0)[z==1])[k]
    rows[[r]] <- do.call(rbind,lapply(names(limits),function(method)
      data.frame(n=n,replication=r,method=method,beta=settings$betas,k=k,
                 lower=limits[[method]],true_quantile=truth,seed=seeds[r,nn])))
    diagrows[[r]] <- data.frame(n=n,replication=r,seed=seeds[r,nn],external_n=n,
      bandwidth=fit$bandwidth,x_knots=length(fit$x_grid),
      outcome_knots=length(fit$outcome_grid),tau_hat=adj$tau_lower,
      median_ess0=adj$median_ess0,median_ess1=adj$median_ess1,
      adjustment_min=min(adj$adjustment),adjustment_max=max(adj$adjustment),
      oracle_min=min(oracle),oracle_max=max(oracle),seconds=proc.time()[[3]]-t0)
    if(r%%25L==0L || r==reps) message("n=",n,": ",r,"/",reps)
  }
  allrows[[nn]] <- do.call(rbind,rows); alldiag[[nn]] <- do.call(rbind,diagrows)
  saveRDS(list(results=allrows[[nn]],diagnostics=alldiag[[nn]],model=model,
               settings=settings,oracle_effect=oracle_effect),
          file.path(outdir,paste0("checkpoint-n-",n,".rds")))
}
sim <- list(results=do.call(rbind,allrows),diagnostics=do.call(rbind,alldiag),
            model=model,settings=settings,oracle_effect=oracle_effect)
saveRDS(sim,file.path(outdir,"simulation.rds"))
write.csv(sim$results,file.path(outdir,"replicates.csv"),row.names=FALSE)
write.csv(sim$diagnostics,file.path(outdir,"diagnostics.csv"),row.names=FALSE)
write.csv(summarize_simulation_comparison(sim),file.path(outdir,"summary.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(outdir,"session-info.txt"))
writeLines(sprintf("Elapsed seconds: %.3f",proc.time()[[3]]-start),file.path(outdir,"timing.txt"))
source("simulation/adjustment-stress/add-rank-methods.R")
