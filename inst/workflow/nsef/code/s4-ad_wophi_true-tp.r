### This script is an approximation starting from true values since this is
### an simulation.

rm(list = ls())

suppressMessages(library(pbdMPI, quietly = TRUE))
init(set.seed = FALSE)
suppressMessages(library(cubfits, quietly = TRUE))

### Set environment.
source("00-set_env.r")
set.seed(simulation$seed)
case.name <- "ad_wophi_true"
case.name <- paste(model, "_", case.name, sep = "")
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)
fn.in <- paste(prefix$data, "simu_true_", model, ".rda", sep = "")
if(file.exists(fn.in)){
  load(fn.in)
} else{
  comm.stop(paste(fn.in, " is not found.", sep = ""))
}
phi.init.true <- EPhi
bInit <- Eb

### Initial.
nIter <- run.info$nIter
burnin <- run.info$burnin 
phi.DrawScale <- run.info$phi.DrawScale

### For configuration.
.CF.DP$dump <- run.info$dump
.CF.DP$prefix.dump <- run.info$prefix.dump
.CF.CT$parallel <- run.info$parallel
.CF.CT$adaptive <- run.info$adaptive

### Run.
phi.init.true <- phi.init.true / mean(phi.init.true)
ret <- cubappr(reu13.df.obs, phi.init.true, y, n,
               nIter = nIter, burnin = burnin,
               bInit = bInit,
               phi.DrawScale = phi.DrawScale,
               model = model, verbose = TRUE, report = 10)

### Dump results.
if(comm.rank() == 0){
  ret.time <- proc.time()
  print(ret.time)

  fn.out <- paste(prefix$output, case.name, "/output_mcmc.rda", sep = "")
  save(list = c("nIter", "burnin", "phi.DrawScale", "ret", "ret.time"),
       file = fn.out)

  fn.out <- paste(prefix$output, case.name, "/output_env.rda", sep = "")
  save(list = ls(envir = .cubfitsEnv),
       file = fn.out, envir = .cubfitsEnv)
}

finalize()
