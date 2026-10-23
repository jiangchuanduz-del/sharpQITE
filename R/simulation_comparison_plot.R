plot_simulation_comparison <- function(summary,sharp_bounds,adjusted=c("fixed","simultaneous","pointwise"),
                         alpha=.05,reps=100L,xlim=NULL) {
  adjusted <- match.arg(adjusted)
  sizes <- sort(unique(summary$n))
  if (length(sizes)>4L || !length(sizes)) stop("Plot one to four sample sizes.")
  st_names <- sort(unique(summary$method[grepl("^s=",summary$method)]))
  st_names <- st_names[order(as.numeric(sub("s=","",st_names)))]
  adj_name <- switch(adjusted,fixed=NULL,simultaneous="RKS adjusted",pointwise="RKS adjusted (pointwise)")
  fixed_name <- "RKS adjusted (beta=0.5)"
  has_fixed <- fixed_name %in% summary$method
  methods <- c("RKS","Wilcoxon",st_names,adj_name,if (has_fixed) fixed_name)
  colors <- c("#BE1E2D","#F18158",rep(c("#218B4D","#2567B2"),length.out=length(st_names)),
               if (length(adj_name)) "#8B47AE",if (has_fixed) "#8B47AE")
  markers <- c(16,17,rep(c(15,18),length.out=length(st_names)),
                if (length(adj_name)) 16,if (has_fixed) 17)
  if (is.null(xlim)) {
    visible <- summary$median[summary$method %in% methods & is.finite(summary$median)]
    xlim <- range(c(stats::quantile(visible,c(.005,.995)),sharp_bounds$marginal,
                    sharp_bounds$covariate))
    xlim <- c(floor(xlim[1]),ceiling(xlim[2]))
  }
  oldpar <- graphics::par(no.readonly=TRUE)
  on.exit(graphics::par(oldpar))
  rows <- ceiling(length(sizes)/2)
  panels <- matrix(c(seq_len(length(sizes))+1L,rep(0L,2L*rows-length(sizes))),ncol=2,byrow=TRUE)
  graphics::layout(rbind(c(1,1),panels),heights=c(.22,rep(1,rows)))
  graphics::par(mar=rep(0,4))
  graphics::plot.new()
  adj_label <- switch(adjusted,fixed=NULL,simultaneous="RKS adjusted (Bonferroni)",
                       pointwise="RKS adjusted (beta pointwise)")
  labels <- c("RKS","Wilcoxon",paste0("Stephenson (",st_names,")"),adj_label,if (has_fixed) fixed_name,
               "sharp: marginals","sharp: covariates")
  legend_colors <- c(colors,"#555555","#111111")
  legend_lty <- c(rep(1,length(methods)),2,3)
  legend_pch <- c(markers,NA,NA)
  # Two rows: test statistics (adjusted RKS immediately after RKS), then references.
  fixed_position <- if (has_fixed) match(fixed_name,labels) else integer()
  method_order <- c(1L,fixed_position,setdiff(seq_along(methods),c(1L,fixed_position)))
  legend_rows <- list(method_order,length(methods)+seq_len(2L))
  legend_y <- c(.7,.27)
  for (r in seq_along(legend_rows)) {
    indices <- legend_rows[[r]]
    widths <- graphics::strwidth(labels[indices],cex=1.35)
    x <- (1-sum(widths+.038)-.045*(length(indices)-1))/2
    for (j in seq_along(indices)) {
      i <- indices[j]
      graphics::segments(x,legend_y[r],x+.023,legend_y[r],
                         col=legend_colors[i],lty=legend_lty[i],lwd=2.3)
      if (!is.na(legend_pch[i]))
        graphics::points(x+.0115,legend_y[r],pch=legend_pch[i],col=legend_colors[i],cex=1)
      graphics::text(x+.038,legend_y[r],labels[i],adj=c(0,.5),cex=1.35)
      x <- x+.038+widths[j]+.045
    }
  }
  # Reserve more label space when callers enlarge axis text.
  axis_extra <- max(0,graphics::par("cex.lab")-1,graphics::par("cex.axis")-1)
  graphics::par(mar=c(4.4,4.5,2.6,1.1)+c(2.5,2.5,0,0)*axis_extra)
  for (n in sizes) {
    title <- if (abs(log10(n)-round(log10(n)))<1e-10) bquote(n==10^.(round(log10(n)))) else bquote(n==.(n))
    graphics::plot(NA,xlim=xlim,ylim=c(0,1),xlab=expression("Lower prediction limit for "*tau[(k)]^t),
                    ylab=expression("Quantile "*k/n[1]),main=title,las=1,xaxs="i",yaxs="i")
    graphics::grid(col="#E8E8E8",lty=1)
    for (i in seq_along(methods)) {
      d <- summary[summary$n==n & summary$method==methods[i],]
      d <- d[order(d$beta),]
      keep <- is.finite(d$q75)
      if (any(keep)) {
        lo <- pmax(xlim[1],d$q25[keep]); hi <- pmin(xlim[2],d$q75[keep])
        graphics::polygon(c(lo,rev(hi)),c(d$beta[keep],rev(d$beta[keep])),
                           col=grDevices::adjustcolor(colors[i],.12),border=NA)
      }
      med <- replace(d$median,!is.finite(d$median),NA_real_)
      graphics::lines(med,d$beta,col=colors[i],lwd=1.9)
      dots <- which(seq_len(nrow(d)) %% 5L == 0L)
      if (length(dots)) graphics::points(med[dots],d$beta[dots],col=colors[i],pch=markers[i],cex=.65)
    }
    graphics::lines(sharp_bounds$marginal,sharp_bounds$beta,col="#555555",lty=2,lwd=2.2)
    graphics::lines(sharp_bounds$covariate,sharp_bounds$beta,col="#111111",lty=3,lwd=2.2)
    graphics::box()
  }
  invisible(summary)
}
