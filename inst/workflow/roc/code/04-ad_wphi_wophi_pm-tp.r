### This is similar to "04-ad_wophi_bInit-tp.r", but the bInit is provided by
### the poster means of "04-ad_wphi_pm-tp.r".

rm(list = ls())

suppressMessages(library(pbdMPI, quietly = TRUE))
init(set.seed = FALSE)
suppressMessages(library(cubfits, quietly = TRUE))

### Set environment.
source("00-set_env.r")
set.seed(simulation$seed)
case.name <- "ad_wphi_wophi_pm"
case.name <- paste(model, "_", case.name, sep = "")

### Check output directory.
fn.out <- paste(prefix$output, case.name, sep = "")
if(!file.exists(fn.out)){
  comm.stop(paste(fn.out, " is not found.", sep = ""))
}

### Load data.
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)

### Load b.PM and phi.PM from previous summarized MCMC outputs.
fn.in <- paste(prefix$subset, model, "_ad_wphi_pm_PM.rda", sep = "")
load(fn.in)
bInit <- convert.bVec.to.b(b.PM, names(reu13.df.obs), model = model)
phi.init.PM <- phi.PM

### Change initial for fitsappr.
nIter <- run.info$nIter
burnin <- run.info$burnin 
phi.DrawScale <- run.info$phi.DrawScale

### For configuration.
.CF.DP$dump <- run.info$dump
.CF.DP$prefix.dump <- run.info$prefix.dump
.CF.CT$parallel <- run.info$parallel
.CF.CT$adaptive <- run.info$adaptive

### Run.
phi.init.PM <- phi.init.PM / mean(phi.init.PM)
ret <- cubappr(reu13.df.obs, phi.init.PM, y, n,
               nIter = nIter, burnin = burnin,
               # bInit = bInit,
               p.nclass = p.nclass,
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
