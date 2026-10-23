# Rscript run_all.R [render|real-data|simulations|all|validate]
files <- Filter(function(x)is.character(x)&&length(x)==1L,lapply(sys.frames(),function(f)f$ofile))
arg <- grep('^--file=',commandArgs(),value=TRUE)
script <- if(length(files))tail(files,1)[[1]] else if(length(arg))sub('^--file=','',arg[1]) else 'run_all.R'
root <- dirname(normalizePath(script,winslash='/',mustWork=TRUE))
run_bundle <- function(mode='render',reps=300L) {
  stopifnot(mode %in% c('render','real-data','simulations','all','validate'))
  old <- setwd(root); on.exit(setwd(old))
  dir.create('.library',showWarnings=FALSE); dir.create('logs',showWarnings=FALSE)
  lib <- normalizePath('.library',winslash='/')
  oldenv <- Sys.getenv(c('R_LIBS_USER','LANG','LC_ALL'),unset=NA_character_)
  on.exit(for(nm in names(oldenv))if(is.na(oldenv[[nm]]))Sys.unsetenv(nm) else do.call(Sys.setenv,setNames(list(oldenv[[nm]]),nm)),add=TRUE)
  Sys.setenv(R_LIBS_USER=lib,LANG='C',LC_ALL='C')
  exe <- function(x) file.path(R.home('bin'),paste0(x,if(.Platform$OS.type=='windows')'.exe' else ''))
  if(!file.exists('.library/sharpQITE/DESCRIPTION')) {
    code <- system2(exe('R'),c('CMD','INSTALL','--no-multiarch',shQuote(paste0('--library=',lib)),'..'),stdout='logs/install.log',stderr='logs/install.log')
    if(code!=0)stop('Install failed: logs/install.log')
  }
  run <- function(path,args=character()) {
    log <- file.path('logs',paste0(gsub('[^A-Za-z0-9]','_',path),'.log'))
    message('Running ',path)
    status <- system2(exe('Rscript'),c('--vanilla',shQuote(path),args),stdout=log,stderr=log)
    if(status!=0){cat(readLines(log),sep='\n');stop('Failed: ',path)}
  }
  if(mode %in% c('real-data','all'))run('real-data/keller-2025/run.R')
  if(mode %in% c('simulations','all')) {
    run('simulation/run.R',as.character(reps))
    run('simulation/adjustment-stress/run.R',as.character(reps))
  }
  if(mode!='validate') {
    run('real-data/keller-2025/vocabulary-table.R')
    run('real-data/keller-2025/vocabulary-full-data.R')
    run('real-data/keller-2025/plot-vocabulary-adjustment-quantiles.R')
    run('simulation/plot-combined-dgps.R')
    run('simulation/plot-population-adjustments.R')
    run('scripts/export-results.R')
  }
  run('scripts/validate.R')
  message('Complete. Open results/index.html.')
}
if(sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  run_bundle(if(length(args))args[1] else 'render',if(length(args)>1)as.integer(args[2])else 300L)
} else run_bundle()
