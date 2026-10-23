# All kernel tuning inputs are from the training half. Covariates in this
# released dataset are complete, so only fold-specific scales are estimated.
tune_keller_kernel <- function(x,z,y,sets,multipliers,seed,folds=5L) {
  stopifnot(nrow(x)==length(y),length(y)==length(z),all(is.finite(x)),all(is.finite(y)))
  set.seed(seed); fold <- integer(length(y))
  for(arm in 0:1) {
    idx <- which(z==arm)
    fold[idx] <- sample(rep(seq_len(folds),length.out=length(idx)))
  }
  scores <- list(); row <- 0L
  for(set_name in names(sets)) {
    xx <- x[,sets[[set_name]],drop=FALSE]
    for(f in seq_len(folds)) {
      tr <- which(fold!=f); va <- which(fold==f)
      bw <- kernel_bandwidth(xx[tr,,drop=FALSE],z[tr],1)
      for(arm in 0:1) {
        support <- tr[z[tr]==arm]; support <- support[order(y[support],support)]
        query <- va[z[va]==arm]
        dist <- kernel_distance(xx[support,,drop=FALSE],xx[query,,drop=FALSE],bw)
        for(multiplier in multipliers) {
          loss <- weighted_crps(kernel_normalize(dist,multiplier),y[support],y[query])
          row <- row+1L
          scores[[row]] <- data.frame(set=set_name,features=ncol(xx),multiplier=multiplier,
            fold=f,arm=arm,n=length(query),crps_sum=sum(loss),crps_mean=mean(loss))
        }
      }
    }
    message('  Validated ',set_name,' (',ncol(xx),' covariates).')
  }
  scores <- do.call(rbind,scores)
  by_arm <- aggregate(cbind(crps_sum,n)~set+features+multiplier+arm,data=scores,FUN=sum)
  by_arm$crps <- by_arm$crps_sum/by_arm$n
  tab <- aggregate(crps~set+features+multiplier,data=by_arm,FUN=mean)
  names(tab)[names(tab)=='crps'] <- 'balanced_cv_crps'
  tab <- tab[order(tab$balanced_cv_crps,tab$features,tab$multiplier),]
  rownames(tab) <- NULL; tab$selected <- seq_len(nrow(tab))==1L
  list(selected=tab[1,],table=tab,scores=scores,by_arm=by_arm,fold=fold,seed=seed)
}
