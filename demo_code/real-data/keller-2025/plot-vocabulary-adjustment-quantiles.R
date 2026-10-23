# Saved training-half models only; no refitting and no inference outcomes used.
local({
  library(sharpQITE)
  source('shared/kernel-helpers.R')
  base <- 'real-data/keller-2025'
  a <- readRDS(file.path(base,'results/analysis.rds'))
  tab <- readRDS(file.path(base,'results/vocabulary-table.rds'))
  models <- list(pretest=a$models$vocabulary,
    pretest_sex=tab$extra_models$vocab_pretest_sex)
  stopifnot(!is.null(models$pretest_sex))
  betas <- c(.1,.3,.5,.7,.9)
  te <- a$models$vocabulary$inference_idx
  observed <- a$data$vocabPre[te]
  grid <- seq.int(ceiling(min(observed)),floor(max(observed)))
  out <- file.path(base,'results/adjustment-quantiles')
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  curves <- inference <- diagnostics <- list()
  for(key in names(models)) {
    model <- models[[key]]; fit <- model$fit
    stopifnot(identical(model$training_idx,a$models$vocabulary$training_idx),
      identical(model$inference_idx,te),length(intersect(model$training_idx,te))==0L)
    panels <- if(key=='pretest')'pretest' else c('sex0','sex1')
    actual_x <- as.matrix(a$data[te,model$feature_columns,drop=FALSE])
    for(beta in betas) {
      adj <- estimate_adjustment_fast(fit,actual_x,beta)
      if(beta==.5) stopifnot(identical(adj$adjustment,model$adjustment$adjustment))
      inference[[paste(key,beta)]] <- data.frame(model=key,beta=beta,
        row_id=a$data$row_id[te],adjustment=adj$adjustment)
      diagnostics[[paste(key,beta)]] <- data.frame(model=key,beta=beta,
        tau=adj$tau_lower,estimated_cdf=adj$estimated_cdf,
        inference_fallbacks=adj$infinite_argmins,fallback=adj$finite_fallback)
      for(panel in panels) {
        xx <- matrix(grid,ncol=1,dimnames=list(NULL,'vocabPre'))
        if(key=='pretest_sex')xx <- cbind(xx,female=as.integer(panel=='sex1'))
        pred <- conditional_min_fast(kernel_pair_fast(fit,xx),adj$tau_lower)
        fallback <- !is.finite(pred$adjustment)
        pred$adjustment[fallback] <- adj$finite_fallback
        # Check every integer-score value against the canonical calculation.
        idx <- seq_len(nrow(xx))
        direct <- getFromNamespace('.conditional_min','sharpQITE')(
          getFromNamespace('.kernel_pair','sharpQITE')(fit,xx[idx,,drop=FALSE]),adj$tau_lower)
        direct$adjustment[!is.finite(direct$adjustment)] <- adj$finite_fallback
        stopifnot(identical(as.numeric(direct$adjustment),as.numeric(pred$adjustment[idx])))
        curves[[paste(panel,beta)]] <- data.frame(panel=panel,model=key,beta=beta,
          pretest=grid,adjustment=pred$adjustment,fallback=fallback,tau=adj$tau_lower)
      }
    }
  }
  curves <- do.call(rbind,curves); diagnostics <- do.call(rbind,diagnostics)
  inference <- do.call(rbind,inference)
  colors <- c('#1265DC','#009C3B','#9932CC','#F0A000','#B2182B')
  symbols <- rep(16,length(betas))
  draw <- function() {
    layout(matrix(c(1,1,1,2,3,4),nrow=2,byrow=TRUE),heights=c(.7,4))
    par(mar=rep(0,4),family='sans',cex=1); plot.new()
    legend('center',legend=as.expression(lapply(betas,function(b)bquote(beta==.(b)))),
      col=colors,lty=1,lwd=1.5,ncol=5,bty='n',cex=1)
    panels <- c('pretest','sex0','sex1')
    titles <- c('Pretest only','Pretest + sex: male','Pretest + sex: female')
    limits <- range(c(curves$adjustment,grid))+c(-.5,.5)
    for(k in seq_along(panels)) {
      g0 <- curves[curves$panel==panels[k],]
      par(mar=c(4.4,4.8,2.4,1),mgp=c(3,.8,0),cex.lab=1,cex.axis=1,cex.main=1.05)
      plot(NA,xlim=range(grid),ylim=limits,
        xlab='Vocabulary pretest score',ylab='Optimal adjustment function',
        main=titles[k],las=1)
      grid(col='#E8E8E8',lty=1)
      abline(a=0,b=1,col='#555555',lwd=1.2,lty=2)
      legend('topleft',legend=expression(y==x),col='#555555',lty=2,lwd=1.2,bty='n',cex=1)
      for(j in seq_along(betas)) {
        g <- g0[g0$beta==betas[j],]
        lines(g$pretest,g$adjustment,col=colors[j],lty=1,lwd=1.5,type='l')
      }
      box()
    }
  }
  stem <- file.path(out,'vocabulary-adjustment-quantiles')
  png(paste0(stem,'.png'),width=4500,height=1740,res=300,pointsize=16); draw(); dev.off()
  svg(paste0(stem,'.svg'),width=15,height=5.8,pointsize=16); draw(); dev.off()
  cairo_pdf(paste0(stem,'.pdf'),width=15,height=5.8,pointsize=16); draw(); dev.off()
  saveRDS(list(curves=curves,inference=inference,diagnostics=diagnostics,
    betas=betas,plot_interpolation='linear between integer pretest scores',
    training_idx=models$pretest$training_idx,inference_idx=te),paste0(stem,'.rds'))
  writeLines(c('PASS: saved models use only the original training half.',
    'PASS: beta=0.5 reproduces saved held-out adjustments exactly for both models.',
    'PASS: every integer pretest score per panel and quantile matches canonical package argmins.',
    'Lines interpolate integer-score evaluations linearly for display.',
    'Quantile-specific curves do not change the existing fixed-beta=0.5 inference tables.',
    'At zero conditional-CDF minima, the plotted finite fallback is a valid nonunique minimizer.'),
    file.path(out,'validation.txt'))
  print(diagnostics,row.names=FALSE)
})
