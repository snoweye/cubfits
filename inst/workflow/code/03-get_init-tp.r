### This script provdes several initial values.

rm(list = ls())

suppressMessages(library(pbdMPI, quietly = TRUE))
init(set.seed = FALSE)
library(cubfits, quiet = TRUE)

### Set environment.
source("00-set_env.r")
set.seed(simulation$seed)
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)

### Initial.
init.function(model = model, parallel = run.info$parallel)

### Run PM initial.
fitlist <- fitMultinom(reu13.df.obs, phi.Obs, y, n)
phi.init.PM <- estimatePhi(fitlist, reu13.df.obs.list, y.list, n.list,
                        E.Phi = median(phi.Obs),
                        lower.optim = min(phi.Obs) * 0.01,
                        upper.optim = max(phi.Obs) * 2.0)
phi.init.PM <- as.double(phi.init.PM)
names(fitlist) <- names(reu13.df.obs)
names(phi.init.PM) <- names(phi.Obs)
### Check NA of phi.init.PM.
id.na <- which(is.na(phi.init.PM))
phi.init.PM[id.na] <- phi.Obs[id.na]

### Run scuo initial.
phi.init.SCUO <- scuo.random(SCUO, meanlog = -simulation$sdlog^2 / 2,
                                 sdlog = simulation$sdlog)
names(phi.init.SCUO) <- names(phi.Obs)

### Scale to Mean 1.
phi.init.PM <- phi.init.PM / mean(phi.init.PM)
phi.init.SCUO <- phi.init.SCUO / mean(phi.init.SCUO)

### Save.
if(comm.rank() == 0){
  comm.size <- comm.size()
  ret.time <- proc.time()
  print(ret.time)

  fn.out <- paste(prefix$data, "init_", model, ".rda", sep = "")
  list.save <- c("fitlist", "phi.init.PM", "phi.init.SCUO",
                 "ret.time", "comm.size")
  save(list = list.save, file = fn.out)
} 

finalize()
