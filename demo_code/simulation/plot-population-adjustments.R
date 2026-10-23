# True conditional Gaussian distributions; no estimated CDFs or simulated fits.
source_dir <- "../R"
source(file.path(source_dir,"simulation_comparison_model.R"))
if (.Platform$OS.type=="windows") windowsFonts(Arial=windowsFont("Arial"))
betas <- c(.1,.3,.5,.7,.9)
xx <- seq(-1,1,length.out=4001)
paths <- c("simulation/results/simulation_comparison.rds",
           "simulation/adjustment-stress/results/simulation.rds")

# Select the local minimum analytically, rather than comparing rounded tail
# probabilities. At a stationary point the second derivative has sign 2*a*s+b.
population_minimum <- function(effect,p) {
  a <- .5*(1/p$sigma0^2-1/p$sigma1^2)
  b <- p$mu1/p$sigma1^2-(effect+p$mu0)/p$sigma0^2
  cc <- .5*((effect+p$mu0)^2/p$sigma0^2-p$mu1^2/p$sigma1^2)+log(p$sigma0/p$sigma1)
  unequal <- abs(a)>1e-12*pmax(1/p$sigma0^2,1/p$sigma1^2)
  disc <- sqrt(pmax(0,b^2-4*a*cc))
  shift <- rep(NA_real_,length(a))
  direct <- unequal & b<=0
  rational <- unequal & b>0
  shift[direct] <- (-b[direct]+disc[direct])/(2*a[direct])
  shift[rational] <- -2*cc[rational]/(b[rational]+disc[rational])
  linear <- !unequal & b>0
  shift[linear] <- -cc[linear]/b[linear]
  value <- numeric(length(a))
  ok <- is.finite(shift)
  z1 <- (shift[ok]-p$mu1[ok])/p$sigma1[ok]
  z0 <- (shift[ok]-effect-p$mu0[ok])/p$sigma0[ok]
  value[ok] <- ifelse(z1>0 & z0>0,
    pnorm(z0,lower.tail=FALSE)-pnorm(z1,lower.tail=FALSE),pnorm(z1)-pnorm(z0))
  list(adjustment=shift,minimum=pmin(0,value))
}
curves <- validations <- list()
for(j in 1:2) {
  model <- readRDS(paths[j])$model
  p <- .model_parameters(model,xx)
  fbar <- function(e,tol=1e-9) 1+integrate(function(x)
    population_minimum(e,.model_parameters(model,x))$minimum,
    -1,1,rel.tol=tol,abs.tol=tol,subdivisions=2000L)$value/2
  for(beta in betas) {
    effect <- uniroot(function(e) fbar(e)-beta,c(-100,100),tol=1e-9)$root
    result <- population_minimum(effect,p)
    old <- .normal_cdf_minimum(effect,p)
    reliable <- old$minimum < -1e-8 & is.finite(old$adjustment)
    stopifnot(max(abs(old$adjustment[reliable]-result$adjustment[reliable]))<1e-7,
              abs(fbar(effect,1e-11)-beta)<1e-7)
    curves[[length(curves)+1]] <- data.frame(dgp=paste("DGP",c("I","II")[j]),
      beta=beta,x=xx,tau_lower=effect,adjustment=result$adjustment)
    validations[[length(validations)+1]] <- data.frame(dgp=j,beta=beta,tau_lower=effect,
      inversion_error=abs(fbar(effect,1e-11)-beta),
      no_finite_minimizer=sum(!is.finite(result$adjustment)))
  }
}
curves <- do.call(rbind,curves)
cols <- c("#2567B2","#218B4D","#8B47AE","#F18158","#BE1E2D")
draw <- function() {
  layout(rbind(c(1,1),c(2,3)),heights=c(.55,4))
  par(mar=rep(0,4),family="Arial",cex=.66); plot.new()
  legend("center",legend=as.expression(lapply(betas,function(b) bquote(beta==.(b)))),
         col=cols,lty=1,lwd=3.2,ncol=5,bty="n",cex=1.3)
  for(j in 1:2) {
    d <- curves[curves$dgp==paste("DGP",c("I","II")[j]),]
    yl <- if(j==1) range(c(-5,floor(min(d$adjustment)),ceiling(max(d$adjustment)))) else c(-15,15)
    par(mar=c(4.4,4.8,2.8,1),mgp=c(3,.8,0),cex.lab=1.2,cex.axis=1.1,cex.main=1.3)
    plot(NA,xlim=c(-1,1),ylim=yl,yaxs="i",xlab=expression(X),
         ylab="Optimal adjustment function",main=paste("DGP",c("I","II")[j]),las=1)
    grid(col="#E8E8E8",lty=1)
    for(i in seq_along(betas)) {
      b <- d[d$beta==betas[i],]
      lines(b$x,b$adjustment,col=cols[i],lwd=3.2)
    }
    box()
  }
}
dir.create("output/pdf",recursive=TRUE,showWarnings=FALSE)
stem <- "output/pdf/population-adjustment-functions-dgps"
write.csv(curves,paste0(stem,".csv"),row.names=FALSE,na="NA")
write.csv(do.call(rbind,validations),paste0(stem,"-validation.csv"),row.names=FALSE)
cairo_pdf(paste0(stem,".pdf"),width=16,height=7,pointsize=22,family="Arial")
draw(); dev.off()
png(paste0(stem,".png"),width=4800,height=2100,res=300,pointsize=22)
draw(); dev.off()
svg(paste0(stem,".svg"),width=16,height=7,pointsize=22,family="Arial")
draw(); dev.off()
