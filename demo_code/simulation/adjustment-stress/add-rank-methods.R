# Add Wilcoxon and Stephenson to the saved experiments without redrawing them.
# Run from the repository root; also called by run.R for fresh simulations.
source_dir <- "../R"
for(f in list.files(source_dir,pattern="[.]R$",full.names=TRUE)) source(f)
outdir <- "simulation/adjustment-stress/results"
sim <- readRDS(file.path(outdir,"simulation.rds"))
rank_names <- c("Wilcoxon","s=6","s=10")
original_results <- sim$results[!sim$results$method %in% rank_names,]
start <- proc.time()[[3]]
cal_path <- file.path(outdir,"rank-calibrations.rds")
if(file.exists(cal_path)) {
  calibrations <- readRDS(cal_path)
} else {
  # Freeze this study's original calibration independently of the other study.
  cal_seeds <- c(`100`=1665141477L,`1000`=1061016257L,
                 `10000`=1966495329L,`100000`=2017171219L)
  calibrations <- lapply(names(cal_seeds),function(n)
    simulation_comparison_calibration(as.integer(n),alpha=.05,s=c(6L,10L),B=9999L,seed=cal_seeds[[n]]))
  names(calibrations) <- names(cal_seeds)
  saveRDS(calibrations,cal_path)
}
allrows <- timings <- list()
for(n in sim$settings$sample_sizes) {
  cal <- calibrations[[as.character(n)]]
  stopifnot(inherits(cal,"sharp_simulation_comparison_calibration"),cal$n==n,
    cal$alpha==sim$settings$alpha,cal$B==9999L,identical(as.integer(cal$s),c(6L,10L)),
    identical(names(cal$scores),rank_names))
  diag <- sim$diagnostics[sim$diagnostics$n==n,]
  diag <- diag[order(diag$replication),]
  check_path <- file.path(outdir,paste0("rank-checkpoint-n-",n,".rds"))
  key <- list(seeds=diag$seed,betas=sim$settings$betas,calibration=cal,
              inversion_tol=.002,model=sim$settings$description)
  cached <- if(file.exists(check_path))readRDS(check_path) else NULL
  if(!is.null(cached) && identical(cached$key,key)) {
    message("Reusing completed rank-method results for n=",n)
    allrows[[as.character(n)]] <- cached$results
    timings[[as.character(n)]] <- cached$timings
    next
  }
  rows <- elapsed <- vector("list",nrow(diag))
  for(i in seq_len(nrow(diag))) {
    t0 <- proc.time()[[3]]; r <- diag$replication[i]; seed <- diag$seed[i]
    set.seed(seed)
    pop <- .draw_gaussian_population(sim$model,n)
    z <- integer(n); z[sample.int(n,n/2)] <- 1L
    y <- ifelse(z==1,pop$y1,pop$y0)
    old <- original_results[original_results$n==n & original_results$replication==r &
                             original_results$method=="RKS",]
    old <- old[order(old$beta),]
    check <- treated_lower_bounds(z,y,sim$settings$betas,alpha=sim$settings$alpha)
    truth <- sort((pop$y1-pop$y0)[z==1])[check$k]
    stopifnot(identical(check$lower,old$lower),identical(truth,old$true_quantile),
              identical(check$k,old$k))
    limits <- list(Wilcoxon=treated_lower_bounds(z,y,sim$settings$betas,
                     statistic="wilcoxon",alpha=sim$settings$alpha,
                     calibration=cal,tol=.002)$lower)
    for(ss in c(6L,10L)) limits[[paste0("s=",ss)]] <-
      treated_lower_bounds(z,y,sim$settings$betas,statistic="stephenson",s=ss,
                           alpha=sim$settings$alpha,calibration=cal,tol=.002)$lower
    rows[[i]] <- do.call(rbind,lapply(names(limits),function(method)
      data.frame(n=n,replication=r,method=method,beta=sim$settings$betas,
                 k=check$k,lower=limits[[method]],true_quantile=truth,seed=seed)))
    elapsed[[i]] <- data.frame(n=n,replication=r,seconds=proc.time()[[3]]-t0)
    if(i==1L || i%%25L==0L || i==nrow(diag))
      message("Rank methods n=",n,": ",i,"/",nrow(diag),"; last ",
               round(elapsed[[i]]$seconds,3)," s")
  }
  allrows[[as.character(n)]] <- do.call(rbind,rows)
  timings[[as.character(n)]] <- do.call(rbind,elapsed)
  saveRDS(list(key=key,results=allrows[[as.character(n)]],
               timings=timings[[as.character(n)]]),check_path)
}
sim$results <- rbind(original_results,do.call(rbind,allrows))
# rbind may change automatic numeric row names to character row names.
# Compare every saved data column, including order, without that display attribute.
preserved <- sim$results[!sim$results$method %in% rank_names,]
stopifnot(identical(names(preserved),names(original_results)),
  all(vapply(names(original_results),function(name)
    identical(preserved[[name]],original_results[[name]]),logical(1))))
rownames(sim$results) <- NULL
sim$calibrations <- calibrations
sim$settings$s <- c(6L,10L)
sim$settings$B <- 9999L
sim$settings$inversion_tol <- .002
sim$settings$rank_calibration_source <- "Frozen reference assignments from the initial 100-replication simulation comparison study"
sim$rank_method_timings <- do.call(rbind,timings)
saveRDS(sim,file.path(outdir,"simulation.rds"))
write.csv(sim$results,file.path(outdir,"replicates.csv"),row.names=FALSE)
write.csv(summarize_simulation_comparison(sim),file.path(outdir,"summary.csv"),row.names=FALSE)
write.csv(sim$rank_method_timings,file.path(outdir,"rank-method-timings.csv"),row.names=FALSE)
for(n in sim$settings$sample_sizes) {
  path <- file.path(outdir,paste0("checkpoint-n-",n,".rds"))
  saved <- readRDS(path)
  saved$results <- sim$results[sim$results$n==n,]
  saved$settings <- sim$settings
  saved$calibration <- calibrations[[as.character(n)]]
  saveRDS(saved,path)
}
writeLines(c(sprintf("Rank-method computation seconds: %.3f",sum(sim$rank_method_timings$seconds)),
             sprintf("Extension elapsed seconds this invocation: %.3f",proc.time()[[3]]-start)),
           file.path(outdir,"rank-method-timing.txt"))


