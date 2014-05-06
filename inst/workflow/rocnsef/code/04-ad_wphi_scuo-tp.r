rm(list = ls())

suppressMessages(library(pbdMPI, quietly = TRUE))
init(set.seed = FALSE)
suppressMessages(library(cubfits, quietly = TRUE))

### Set environment.
source("00-set_env.r")
set.seed(simulation$seed)
case.name <- "ad_wphi_scuo"
case.name <- paste(model, "_", case.name, sep = "")

### Check output directory.
fn.out <- paste(prefix$output, case.name, sep = "")
if(!file.exists(fn.out)){
  comm.stop(paste(fn.out, " is not found.", sep = ""))
}

### Load data.
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)
fn.in <- paste(prefix$data, "init_", model, ".rda", sep = "")
load(fn.in)

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
phi.Obs <- phi.Obs / mean(phi.Obs)
phi.init.SCUO <- phi.init.SCUO / mean(phi.init.SCUO)
ret <- cubfits(reu13.df.obs, phi.Obs, y, n,
               nIter = nIter, burnin = burnin,
               p.nclass = p.nclass,
               phi.Init = phi.init.SCUO,
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

