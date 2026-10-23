# Run from the package root after installing sharpQITE >= 0.2.2.
# Rscript simulation/run.R [replications, default 300]
library(sharpQITE)
args <- commandArgs(trailingOnly=TRUE)
reps <- if (length(args)) as.integer(args[1]) else 300L
dir.create("simulation/results",recursive=TRUE,showWarnings=FALSE)
start <- proc.time()[[3]]
model <- gaussian_covariate_model()
sharp <- superpopulation_bounds(model)
write.csv(sharp,"simulation/results/superpopulation-bounds.csv",row.names=FALSE)
sim <- simulate_simulation_comparison(reps=reps,model=model,seed=20261004,B=9999,
                        checkpoint_dir="simulation/results/checkpoints",progress=TRUE)
summary <- summarize_simulation_comparison(sim)
saveRDS(sim,"simulation/results/simulation_comparison.rds")
write.csv(sim$results,"simulation/results/replicates.csv",row.names=FALSE)
write.csv(sim$diagnostics,"simulation/results/diagnostics.csv",row.names=FALSE)
write.csv(summary,"simulation/results/summary.csv",row.names=FALSE)
writeLines(capture.output(sessionInfo()),"simulation/results/session-info.txt")
writeLines(sprintf("Total elapsed seconds: %.3f",proc.time()[[3]]-start),"simulation/results/timing.txt")

