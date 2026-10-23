# Run from the project root: Rscript simulation/plot-combined-dgps.R
# Combine the saved Figure 1 and Figure 2 results; no simulations are rerun.
sizes <- c(100L, 1000L, 100000L)
if (.Platform$OS.type == "windows") windowsFonts(Arial=windowsFont("Arial"))
methods <- c("RKS", "RKS adjusted (beta=0.5)", "Wilcoxon", "s=6")
colors <- c("#BE1E2D", "#8B47AE", "#F18158", "#218B4D")
labels <- expression("RKS", "RKS adjusted ("*beta[0]==0.5*")",
                     "Wilcoxon", "Stephenson (s=6)")
read_dgp <- function(path, label, xlim = NULL) {
  s <- read.csv(file.path(path, "summary.csv"))
  b <- read.csv(file.path(path, "superpopulation-bounds.csv"))
  # Preserve the original figure's automatic range, including all saved sizes.
  if (is.null(xlim)) {
    visible <- s$median[s$method %in% methods & is.finite(s$median)]
    limits <- range(c(quantile(visible, c(.005, .995)), b$marginal, b$covariate))
    xlim <- c(floor(limits[1]), ceiling(limits[2]))
  }
  s <- s[s$n %in% sizes & s$method %in% methods, ]
  stopifnot(all(table(factor(s$n, levels=sizes), factor(s$method, levels=methods)) == 49L))
  list(summary=s, bounds=b, label=label, xlim=xlim)
}
dgps <- list(
  read_dgp("simulation/results", "DGP I"),
  read_dgp("simulation/adjustment-stress/results", "DGP II", c(-10, 6))
)

draw_combined <- function() {
  layout(rbind(c(1,1,1), c(2,2,2), c(3,4,5), c(6,6,6), c(7,8,9)),
         heights=c(.45,.3,3.6,.3,3.6))
  par(mar=rep(0,4), family="Arial")
  plot.new()
  legend_labels <- c(labels, expression("sharp: marginals", "sharp: covariates"))
  legend_colors <- c(colors, "#555555", "#111111")
  legend_lty <- c(rep(1, length(methods)), 2, 3)
  widths <- strwidth(legend_labels, cex=1.3)
  gap <- .018
  key_width <- .028
  total_width <- sum(widths+key_width)+gap*(length(legend_labels)-1)
  stopifnot(total_width < 1)
  x <- (1-total_width)/2
  for (i in seq_along(legend_labels)) {
    segments(x,.5,x+.019,.5,col=legend_colors[i],lty=legend_lty[i],lwd=3.6)
    text(x+key_width,.5,legend_labels[i],adj=c(0,.5),cex=1.3)
    x <- x+key_width+widths[i]+gap
  }
  for (dgp in dgps) {
    par(mar=rep(0,4)); plot.new()
    text(.5,.5,dgp$label,font=2,cex=1.2)
    for (n in sizes) {
      par(mar=c(4.4,4.5,2.5,1), mgp=c(2.9,.8,0), cex.lab=1.15, cex.axis=1.1, cex.main=1.25)
      plot(NA, xlim=dgp$xlim, ylim=c(0,1),
           xlab=expression("Lower prediction limit for "*tau[(k)]^t),
           ylab=expression("Quantiles "*beta==k/n[1]),
           main=bquote(n==10^.(round(log10(n)))), las=1, xaxs="i", yaxs="i")
      grid(col="#E8E8E8", lty=1)
      for (i in seq_along(methods)) {
        d <- dgp$summary[dgp$summary$n==n & dgp$summary$method==methods[i], ]
        d <- d[order(d$beta), ]
        keep <- is.finite(d$q75)
        if (any(keep)) {
          lo <- pmax(dgp$xlim[1],d$q25[keep])
          hi <- pmin(dgp$xlim[2],d$q75[keep])
          polygon(c(lo,rev(hi)),c(d$beta[keep],rev(d$beta[keep])),
                  col=adjustcolor(colors[i],.25),border=NA)
        }
        med <- replace(d$median,!is.finite(d$median),NA_real_)
        lines(med,d$beta,col=colors[i],lwd=3.2)
      }
      b <- dgp$bounds
      lines(b$marginal,b$beta,col="#555555",lty=2,lwd=3.2)
      lines(b$covariate,b$beta,col="#111111",lty=3,lwd=3.2)
      box()
    }
  }
}

dir.create("output/pdf", recursive=TRUE, showWarnings=FALSE)
stem <- "output/pdf/combined-simulation-dgps"
cairo_pdf(paste0(stem,".pdf"), width=16, height=10, family="Arial", pointsize=22)
draw_combined(); dev.off()
png(paste0(stem,".png"), width=4800, height=3000, res=300, pointsize=22)
draw_combined(); dev.off()
svg(paste0(stem,".svg"), width=16, height=10, family="Arial", pointsize=22)
draw_combined(); dev.off()
cat("Created", stem, "in PDF, PNG, and SVG formats.\n")



