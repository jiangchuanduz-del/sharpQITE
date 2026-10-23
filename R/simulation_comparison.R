simulate_simulation_comparison <- function(sample_sizes=c(100L,1000L,10000L,100000L),reps=100L,
                             betas=seq(.02,.98,by=.02),alpha=.05,s=c(6L,10L),
                             external_ratio=1,B=9999L,model=gaussian_covariate_model(),
                             inversion_tol=.002,adjustment_tol=.005,seed=20261004L,
                             checkpoint_dir=NULL,progress=interactive()) {
  sample_sizes <- vapply(sample_sizes,.integer_scalar,integer(1),name="sample_sizes",lower=4)
  if (any(sample_sizes%%2 != 0) || anyDuplicated(sample_sizes)) stop("sample_sizes must be distinct even integers.")
  reps <- .integer_scalar(reps,"reps")
  if (!is.numeric(betas) || !length(betas) || anyNA(betas) || any(betas<=0 | betas>=1) ||
      any(diff(betas)<=0)) stop("betas must be increasing in (0,1).")
  alpha <- .probability(alpha,"alpha")
  if (!is.numeric(external_ratio) || length(external_ratio)!=1 || !is.finite(external_ratio) || external_ratio<=0)
    stop("external_ratio must be positive.")
  if (!inherits(model,"sharp_gaussian_model")) stop("Use gaussian_covariate_model().")
  if (!is.null(checkpoint_dir)) dir.create(checkpoint_dir,recursive=TRUE,showWarnings=FALSE)
  .with_seed(seed, {
    rows <- list(); diagnostics <- list(); calibrations <- list()
    for (n in sample_sizes) {
      if (progress) message("Calibrating n=",n," using ",B," assignments ...")
      cal_seed <- sample.int(.Machine$integer.max,1)
      cal <- simulation_comparison_calibration(n,alpha,s,B,cal_seed)
      calibrations[[as.character(n)]] <- cal
      raw <- diag_rows <- vector("list",reps)
      for (r in seq_len(reps)) {
        rep_seed <- sample.int(.Machine$integer.max,1)
        output <- .with_seed(rep_seed, {
          start <- proc.time()[[3]]
          pop <- .draw_gaussian_population(model,n)
          z <- integer(n); z[sample.int(n,n/2)] <- 1L
          y <- ifelse(z==1,pop$y1,pop$y0)
          nt <- max(4L,2L*ceiling(n*external_ratio/2))
          train <- .draw_gaussian_population(model,nt)
          zt <- integer(nt); zt[sample.int(nt,nt/2)] <- 1L
          yt <- ifelse(zt==1,train$y1,train$y0)
          x_grid <- seq(model$xlim[1],model$xlim[2],length.out=max(41L,2L*ceiling(nt^.3)+1L))
          fit <- fit_kernel_cdf_grid(train$x,zt,yt,x_grid=x_grid)
          adjustment_betas <- sort(unique(c(betas,.5)))
          adj <- estimate_grid_adjustment(fit,pop$x,adjustment_betas,adjustment_tol)
          res <- list()
          res[["RKS"]] <- treated_lower_bounds(z,y,betas,alpha=alpha)$lower
          res[["Wilcoxon"]] <- treated_lower_bounds(z,y,betas,"wilcoxon",alpha=alpha,
                                                     calibration=cal,tol=inversion_tol)$lower
          for (ss in s) res[[paste0("s=",ss)]] <- treated_lower_bounds(z,y,betas,"stephenson",ss,
                                                   alpha,cal,inversion_tol)$lower
          # One externally trained adjustment is held fixed for the entire band.
          fixed_y <- y-adj$adjustment[,match(.5,adjustment_betas)]
          res[["RKS adjusted (beta=0.5)"]] <- treated_lower_bounds(z,fixed_y,betas,alpha=alpha)$lower
          res[["RKS adjusted"]] <- res[["RKS adjusted (pointwise)"]] <- numeric(length(betas))
          for (j in seq_along(betas)) {
            ya <- y-adj$adjustment[,match(betas[j],adjustment_betas)]
            res[["RKS adjusted"]][j] <- treated_lower_bounds(z,ya,betas[j],alpha=alpha/length(betas))$lower
            res[["RKS adjusted (pointwise)"]][j] <- treated_lower_bounds(z,ya,betas[j],alpha=alpha)$lower
          }
          k <- ceiling(.snap_integer(n/2*betas))
          # Valid monotone improvement on the finite grid after Bonferroni.
          res[["RKS adjusted"]] <- cummax(res[["RKS adjusted"]])
          truth <- sort((pop$y1-pop$y0)[z==1])[k]
          out <- do.call(rbind,lapply(names(res),function(name)
            data.frame(n=n,replication=r,method=name,beta=betas,k=k,
                       lower=res[[name]],true_quantile=truth,seed=rep_seed)))
          diag <- data.frame(n=n,replication=r,external_n=nt,seed=rep_seed,
                             seconds=proc.time()[[3]]-start,bandwidth=fit$bandwidth,
                             x_knots=length(fit$x_grid),outcome_knots=length(fit$outcome_grid),
                             median_ess0=adj$median_ess0,median_ess1=adj$median_ess1)
          list(results=out,diagnostics=diag)
        })
        raw[[r]] <- output$results; diag_rows[[r]] <- output$diagnostics
        if (!is.null(checkpoint_dir) && (r%%10L==0L || r==reps))
          saveRDS(list(results=do.call(rbind,raw[seq_len(r)]),
                       diagnostics=do.call(rbind,diag_rows[seq_len(r)]),calibration=cal),
                   file.path(checkpoint_dir,paste0("n-",n,".rds")))
        if (progress && (r==1 || r%%10L==0L || r==reps))
          message("n=",n,": ",r,"/",reps,"; last replication ",round(output$diagnostics$seconds,2)," s")
      }
      rows[[as.character(n)]] <- do.call(rbind,raw)
      diagnostics[[as.character(n)]] <- do.call(rbind,diag_rows)
    }
    list(results=do.call(rbind,rows),diagnostics=do.call(rbind,diagnostics),calibrations=calibrations,
         model=model,settings=list(sample_sizes=sample_sizes,reps=reps,betas=betas,alpha=alpha,s=s,
           external_ratio=external_ratio,B=B,inversion_tol=inversion_tol,adjustment_tol=adjustment_tol,
           seed=seed,adjusted_alpha=alpha/length(betas),fixed_adjustment_beta=.5,
           fixed_adjustment_alpha=alpha,target="treated",external_redrawn=TRUE))
  })
}

summarize_simulation_comparison <- function(simulation) {
  d <- simulation$results
  groups <- split(seq_len(nrow(d)),interaction(d$n,d$method,d$beta,drop=TRUE))
  ans <- lapply(groups,function(i) {
    a <- d[i,,drop=FALSE]
    # Type 1 includes -Inf as an extended-real observation and never averages
    # -Inf with finite values. Dropping infinite replicates would bias curves.
    q <- stats::quantile(a$lower,c(.25,.5,.75),type=1,names=FALSE)
    data.frame(n=a$n[1],method=a$method[1],beta=a$beta[1],
                q25=q[1],median=q[2],q75=q[3],finite_fraction=mean(is.finite(a$lower)))
  })
  out <- do.call(rbind,ans); rownames(out)<-NULL
  out[order(out$n,out$method,out$beta),]
}
