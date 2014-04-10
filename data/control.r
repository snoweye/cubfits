# cubfits control variables.

# For method controls.
.CF.CT <- list(
  model = c("roc"),                             # main models
  type.B = c("Norm"),                           # proposal for M and S
  type.p = c("lognormal_fix", "lognormal",
             "lognormal_MG", "lognormal_MG0",
             "fixed_SM"),                       # proposal for hyperparameters
  type.Phi = c("RW_Norm"),                      # proposal for Phi
  model.Phi = c("lognormal"),                   # prior of Phi
  init.Phi = c("PM"),                           # initial methods for Phi
  init.fit = c("current", "random"),            # how is coef initialed in VGAM
  parallel = c("lapply", "mclapply",
               "task.pull", "pbdLapply"),       # parallel functions
  adaptive = c("simple", "none"),               # method for adaptive mcmc
  scale.Phi = c("mean_one", "none")             # scaling for Phi
)

# For optimization.
.CF.OP <- list(
  optim.method = c("Brent"),                       # for optim()
  stable.min.exp = .Machine$double.max.exp * 0.1,  # minimum exponent
  stable.max.exp = .Machine$double.max.exp * 0.5,  # maximum exponent
  E.Xg = 1.0,                                      # expected Phi
  lower.optim = 1e-4,                              # lower of d logL(x)
  upper.optim = 1e2,                               # upper of d logL(x)
  lower.integrate = 0.0,                           # lower of \int L(x)
  upper.integrate = Inf                            # upper of \int L(x)
)

# For dumpping data.
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

# For addaptive control.
.CF.AC <- list(
  renew.iter = 100,                # per renewing iterations
  target.accept.lower = 0.25,      # target acceptance lower bound
  target.accept.upper = 0.5,       # target acceptance upper bound
  scale.increase = 1.1,            # 10% more
  scale.decrease = 0.9,            # 10% less
  sigma2.lower = 0.01,             # lower bound of sigma^2
  sigma2.upper = 100               # upper bound of sigma^2
)
