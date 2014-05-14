### cubfits control variables.

### For method controls.
.CF.CT <- list(
  model = c("roc", "nsef", "rocnsef"),          # main models
  type.p = c("lognormal_fix",
             "lognormal_MH",
             "lognormal_RW",
             "logmixture"),                     # proposal for hyperparameters
  type.Phi = c("RW_Norm"),                      # proposal for Phi
  model.Phi = c("lognormal", "logmixture"),     # prior of Phi
  init.Phi = c("PM"),                           # initial methods for Phi
  init.fit = c("current", "random", "RW_Norm"), # how is beta proposed
  parallel = c("lapply", "mclapply",
               "task.pull", "pbdLapply"),       # parallel functions
  adaptive = c("simple", "none"),               # method for adaptive mcmc
  scale.Phi = c("mean_one", "none")             # scaling for Phi
)

### For optimization.
.CF.OP <- list(
  optim.method = c("Brent"),                       # for optim()
  stable.min.exp = .Machine$double.max.exp * 0.1,  # minimum exponent
  stable.max.exp = .Machine$double.max.exp * 0.5,  # maximum exponent
  lower.optim = 1e-4,                              # lower of d logL(x)
  upper.optim = 1e2,                               # upper of d logL(x)
  lower.integrate = 0.0,                           # lower of \int L(x)
  upper.integrate = Inf                            # upper of \int L(x)
)

### For dumpping data.
.CF.DP <- list(
  dump = FALSE,                    # if dumping within MCMC
  iter = 1000,                     # iterations per dumping
  prefix.dump = "dump_",           # path and file names of dumping
  trace.acceptance = TRUE,         # if trace acceptance rate
  verbose = FALSE,                 # if verbose
  iterThin = 1,                    # iterations to thin chain
  report = 10,                     # iterations to report
  report.proc = 100                # iterations to report proc.time()
)

### For addaptive control.
.CF.AC <- list(
  renew.iter = 100,                # per renewing iterations
  target.accept.lower = 0.25,      # target acceptance lower bound
  target.accept.upper = 0.5,       # target acceptance upper bound
  scale.increase = 1.2,            # 20% more
  scale.decrease = 0.8,            # 20% less
  sigma2.lower = 0.01,             # lower bound of sigma^2
  sigma2.upper = 100               # upper bound of sigma^2
)

### For parameters as reestimated for Yeast according to Yassour's data.
.CF.PARAM <- list(
  phi.meanlog = -1.125,            # mean of phi in log scale
  phi.sdlog = 1.5,                 # standard deviation of phi in log scale
  # hp.gamma.mean = 0.966,           # mean of gamma distribution for sdlog
  # hp.gamma.sd = 0.01375,           # sd of gamma distribution for sdlog
  # hp.gamma.var = 0.0001890625,     # var of gamma distribution for sdlog
  hp.gamma.shape = 4935.701,       # 0.966^2 / 0.0001890625
  hp.gamma.scale = 0.0001957169,  # 0.0001890625 / 0.966
  hp.overwrite = TRUE              # if allow my.pInit() to overwrite
)
### Gamma has mean = alpha * beta and var = alpha * beta^2
### where alpha and beta are the shape and scale parameters in R,
### resepectively.

### For configuration of initial and draw scaling.
.CF.CONF <- list(
  init.b.Scale = 1,                # initial b scale
  init.phi.Scale = 1,              # initial phi scale
  p.nclass = 2,                    # # of classes if mixture phi
  b.DrawScale = 1,                 # drawing scale for b if random walk
  p.DrawScale = 0.1,               # drawing scale for p if random walk
  phi.DrawScale = 1,               # random walk scale for phi
  phi.DrawScale.pred = 1           # random walk scale for phi.pred
)
