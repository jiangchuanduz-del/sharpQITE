local({
  base <- 'real-data/keller-2025'
  a <- readRDS(file.path(base,'results/vocabulary-table.rds'))
  b <- readRDS(file.path(base,'full-data-results/vocabulary-table.rds'))
  stopifnot(b$settings$fold_combination=='pooled order statistics')
  dir.create('results',showWarnings=FALSE)
  names <- c('Wilcoxon','Stephenson (s=6)','RKS','RKS + naive (pretest)',
    'RKS + optimal (pretest)','RKS + optimal (pretest + sex)')
  fmt <- function(x)ifelse(is.infinite(x),'-Inf',formatC(x,format='f',digits=2))
  panels <- list(t(a$lower),t(b$lower))
  colnames(panels[[1]]) <- colnames(panels[[2]]) <- paste0(10*1:10,'%')
  titles <- c('Panel A: Inference sample (n = 1,100)','Panel B: Full sample (n = 2,200)')
  draw <- function() {
    par(mar=rep(0,4),family='serif'); plot.new(); plot.window(c(0,1),c(0,1),xaxs='i',yaxs='i')
    text(.5,.965,'Table 1',cex=1.1)
    text(.5,.915,'90% simultaneous lower confidence bounds for individual treatment-effect quantiles',font=3,cex=1)
    xx <- c(.025,seq(.33,.975,length.out=10))
    text(xx[1],.853,'Method',adj=0,cex=.95)
    text(xx[-1],.853,paste0(10*1:10,'%'),cex=.95)
    segments(.025,.885,.985,.885);segments(.025,.827,.985,.827)
    for(f in 1:2) {
      yy <- if(f==1).79 else .415
      text(.025,yy,titles[f],adj=0,font=3,cex=1)
      segments(.025,yy-.022,.985,yy-.022)
      for(j in 1:6) {
        y <- yy-j*.05
        text(.025,y,names[j],adj=0,cex=.90)
        for(k in 1:10) {
          v <- panels[[f]][j,k]
          text(xx[k+1],y,if(is.infinite(v))expression(-infinity)else fmt(v),cex=.9)
        }
      }
      segments(.025,yy-.325,.985,yy-.325)
    }
  }
  cairo_pdf('results/02-vocabulary-table.pdf',width=10,height=4.7,pointsize=12);draw();dev.off()
  png('results/02-vocabulary-table.png',width=2400,height=1128,res=240,pointsize=12);draw();dev.off()
  svg('results/02-vocabulary-table.svg',width=10,height=4.7,pointsize=12);draw();dev.off()
  rows <- '<tr><th>Method</th><th>10%</th><th>20%</th><th>30%</th><th>40%</th><th>50%</th><th>60%</th><th>70%</th><th>80%</th><th>90%</th><th>100%</th></tr>'
  for(f in 1:2) {
    rows <- c(rows,paste0('<tr><th colspan="11">',titles[f],'</th></tr>'))
    for(j in 1:6)rows <- c(rows,paste0('<tr><td>',names[j],'</td><td>',paste(fmt(panels[[f]][j,]),collapse='</td><td>'),'</td></tr>'))
  }
  notes <- paste('Vocabulary versus math training; outcome: vocabulary posttest score.',
    'Panel A: all methods infer on the same 1,100 held-out participants.',
    'Panel B: untrained methods use all 2,200 directly; learned adjustments use both training/inference directions.',
    'Bonferroni allocates alpha=0.05 per half and 0.025 per arm within each half;',
    'the half-sample sequences are pooled and sorted. This is the screenshot-selected pooled version.',
    'Optimal adjustments use beta_0=0.5; sex is the female indicator.',
    'Coverage is simultaneous across quantiles within each method, not across methods. Values rounded to two decimals.')
  writeLines(c('<html><meta charset="utf-8"><title>Vocabulary table</title><style>body{font-family:Georgia;margin:35px}td,th{padding:7px;text-align:right;border-bottom:1px solid #ddd}td:first-child{text-align:left}table{border-collapse:collapse}p{max-width:1100px;line-height:1.5}</style>',
    '<h1>Vocabulary training</h1><h2>90% simultaneous lower confidence bounds</h2><table>',rows,'</table>',paste0('<p>',notes,'</p></html>')),'results/02-vocabulary-table.html')
  saveRDS(list(panel_a=a$lower,panel_b=b$lower,quantiles=a$quantile,notes=notes),'results/02-vocabulary-table.rds')
  stems <- c('output/pdf/combined-simulation-dgps',
    'real-data/keller-2025/results/adjustment-quantiles/vocabulary-adjustment-quantiles',
    'output/pdf/population-adjustment-functions-dgps')
  targets <- c('01-simulation-comparison','03-vocabulary-adjustment','04-population-adjustments')
  for(i in 1:3)for(ext in c('.pdf','.png','.svg'))stopifnot(file.copy(paste0(stems[i],ext),paste0('results/',targets[i],ext),overwrite=TRUE))
  writeLines(c('<html><meta charset="utf-8"><title>Four retained results</title><style>body{font-family:Georgia;max-width:1100px;margin:30px auto}img{max-width:100%}section{margin:35px 0}a{margin-right:15px}</style><h1>Four retained results</h1>',
    vapply(c('01-simulation-comparison','02-vocabulary-table','03-vocabulary-adjustment','04-population-adjustments'),function(s)paste0('<section><h2>',s,'</h2><a href="',s,'.pdf">PDF</a><a href="',s,'.svg">SVG</a><img src="',s,'.png"></section>'),character(1)),
    '<p>Table 1 uses the selected pooled Panel B bounds. See README.md for reproduction and inference details.</p></html>'),'results/index.html')
})
