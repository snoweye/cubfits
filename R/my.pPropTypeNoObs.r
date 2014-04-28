# Gibbs sampler for step 2, mainly for hyperparameters.

# Get the specific function according to the options.
get.my.pPropTypeNoObs <- function(type){
  if(!any(type[1] %in% .CF.CT$type.p)){
    stop("type is not found.")
  }
  ret <- eval(parse(text = paste("my.pPropTypeNoObs.", type[1], sep = "")))
  assign("my.pPropTypeNoObs", ret, envir = .cubfitsEnv)
  ret
} # End of get.my.pPropTypeNoObs().


# Drew Gibbs Sampler given current status for lognormal prior around fixed
# mean of log expression.
my.pPropTypeNoObs.lognormal_fix <- function(n.G, log.Phi.Obs.mean, log.Phi.Curr,
    nu.Phi.Curr, sigma.Phi.Curr,
    p.DrawScale = 1, p.DrawScale.prev = 1, Phi.Curr = NULL){
  # nu.Phi.Curr is unused in this case.
  # sigma.Phi.Curr is unused in this case.

  ### The next is assuming non-informative priors for mu.Phi and sigma.Phi.sq.
  ### i.e. p(mu.Phi, sigma.Phi.sq) \propto 1/sigma.Phi.sq.
  ### See Gelman et al. (2003), p.75 for details.

  ## Draw \sigma^{2*}_phi from IG((n_G - 1) / 2, (n_G - 1) S^{2(t)}_{phi} / 2)
  # sigma.Phi.Curr <- sqrt(1 / rgamma(1, shape = (n_G - 1) / 2,
  #                                  rate = sum((log.Phi.Curr - nu.Phi.Curr)^2) / 2))
  ## Draw \mu^*_phi from N(sum_g log(phi^{(t)}_g) / n_G, \sigma^{2*}_phi / n_G)
  # nu.Phi.Curr <- mean(log.Phi.Curr) + rnorm(1) * sigma.Phi.Curr / sqrt(n.G)

  ### Since there is no phi.Obs, log.Phi.Curr can draft potentially to anywhere.
  ### This is to avoid draft of Phi where log.Phi.Obs.mean is fixed for all
  ### iterations, but we still allow some uncertainty around the mean
  ### value for all Phi's.
  sigma.Phi.Curr <-
    sqrt(1 / rgamma(1, shape = (n.G - 1) / 2,
                       rate = sum((log.Phi.Curr - log.Phi.Obs.mean)^2) / 2))
  nu.Phi.Curr <- log.Phi.Obs.mean + rnorm(1) * sigma.Phi.Curr / sqrt(n.G)

  # Update prior's acceptance and adaptive.
  accept <- 1
  my.update.acceptance("p", accept)
  my.update.adaptive("p", accept)

  # Only nu.Phi and sigma.Phi are used.
  ret <- list(nu.Phi = nu.Phi.Curr, sigma.Phi = sigma.Phi.Curr)
  ret
} # End of my.pPropTypeNoObs.lognormal_fix().

# Drew Gibbs Sampler given current status for lognormal prior.
my.pPropTypeNoObs.lognormal <- function(n.G, log.Phi.Obs, log.Phi.Curr,
    nu.Phi.Curr, sigma.Phi.Curr,
    p.DrawScale = 1, p.DrawScale.prev = 1, Phi.Curr = NULL){
  # nu.Phi.Curr is unused in this case.
  # sigma.Phi.Curr is unused in this case.

  ### The next is assuming non-informative priors for mu.Phi and sigma.Phi.sq.
  ### i.e. p(mu.Phi, sigma.Phi.sq) \propto 1/sigma.Phi.sq.
  ### See Gelman et al. (2003), p.75 for details.

  ## Draw \sigma^{2*}_phi from IG((n_G - 1) / 2, (n_G - 1) S^{2(t)}_{phi} / 2)
  sigma.Phi.Curr <-
    sqrt(1 / rgamma(1, shape = (n.G - 1) / 2,
                       rate = sum((log.Phi.Curr - log.Phi.Obs)^2) / 2))
  ## Draw \mu^*_phi from N(sum_g log(phi^{(t)}_g) / n_G, \sigma^{2*}_phi / n_G)
  nu.Phi.Curr <- mean(log.Phi.Curr) + rnorm(1) * sigma.Phi.Curr / sqrt(n.G)

  # Update prior's acceptance and adaptive.
  accept <- 1
  my.update.acceptance("p", accept)
  my.update.adaptive("p", accept)

  # Only nu.Phi and sigma.Phi are used.
  ret <- list(nu.Phi = nu.Phi.Curr, sigma.Phi = sigma.Phi.Curr)
  ret
} # End of my.pPropTypeNoObs.lognormal().

# This method violates MCMC fundamental assumption.
my.pPropTypeNoObs.lognormal_MG <- function(n.G, log.Phi.Obs.mean, log.Phi.Curr,
    nu.Phi.Curr, sigma.Phi.Curr,
    p.DrawScale = 1, p.DrawScale.prev = 1, Phi.Curr = NULL){
  # nu.Phi.Curr is unused in this case.
  # sigma.Phi.Curr is unused in this case.

  ### The next is assuming non-informative priors for mu.Phi and sigma.Phi.sq.
  ### i.e. p(mu.Phi, sigma.Phi.sq) \propto 1/sigma.Phi.sq.
  ### See Gelman et al. (2003), p.75 for details.

  ## Draw \sigma^{2*}_phi from IG((n_G - 1) / 2, (n_G - 1) S^{2(t)}_{phi} / 2)
  # sigma.Phi.Curr <- sqrt(1 / rgamma(1, shape = (n_G - 1) / 2,
  #                                  rate = sum((log.Phi.Curr - nu.Phi.Curr)^2) / 2))
  ## Draw \mu^*_phi from N(sum_g log(phi^{(t)}_g) / n_G, \sigma^{2*}_phi / n_G)
  # nu.Phi.Curr <- mean(log.Phi.Curr) + rnorm(1) * sigma.Phi.Curr / sqrt(n_G)

  ### Since there is no phi.Obs, log.Phi.Curr can draft potentially to anywhere.
  ### This is to avoid draft of Phi where log.Phi.Obs.mean is fixed for all
  ### iterations, but we still allow some uncertainty around the mean
  ### value for all Phi's.
  sigma.Phi.Curr <-
    sqrt(1 / rgamma(1, shape = (n.G - 1) / 2,
                       rate = sum((log.Phi.Curr - log.Phi.Obs.mean)^2) / 2))
  # nu.Phi.Curr <- log.Phi.Obs.mean + rnorm(1) * sigma.Phi.Curr / sqrt(n.G)
  nu.Phi.Curr <- -sigma.Phi.Curr^2 / 2
  .cubfitsEnv$my.print(sigma.Phi.Curr)

  # Update prior's acceptance and adaptive.
  accept <- 1
  my.update.acceptance("p", accept)
  my.update.adaptive("p", accept)

  # Only nu.Phi and sigma.Phi are used.
  ret <- list(nu.Phi = nu.Phi.Curr, sigma.Phi = sigma.Phi.Curr)
  ret
} # End of my.pPropTypeNoObs.lognormal_MG().

# This method violates MCMC fundamental assumption.
my.pPropTypeNoObs.lognormal_MG0 <- function(n.G, log.Phi.Obs.mean, log.Phi.Curr,
    nu.Phi.Curr, sigma.Phi.Curr,
    p.DrawScale = 1, p.DrawScale.prev = 1, Phi.Curr = NULL){
  # nu.Phi.Curr is unused in this case.
  # sigma.Phi.Curr is unused in this case.

  ### The next is assuming non-informative priors for mu.Phi and sigma.Phi.sq.
  ### i.e. p(mu.Phi, sigma.Phi.sq) \propto 1/sigma.Phi.sq.
  ### See Gelman et al. (2003), p.75 for details.

  ## Draw \sigma^{2*}_phi from IG((n_G - 1) / 2, (n_G - 1) S^{2(t)}_{phi} / 2)
  # sigma.Phi.Curr <- sqrt(1 / rgamma(1, shape = (n_G - 1) / 2,
  #                                  rate = sum((log.Phi.Curr - nu.Phi.Curr)^2) / 2))
  ## Draw \mu^*_phi from N(sum_g log(phi^{(t)}_g) / n_G, \sigma^{2*}_phi / n_G)
  # nu.Phi.Curr <- mean(log.Phi.Curr) + rnorm(1) * sigma.Phi.Curr / sqrt(n.G)

  ### Since there is no phi.Obs, log.Phi.Curr can draft potentially to anywhere.
  ### This is to avoid draft of Phi where log.Phi.Obs.mean is fixed for all
  ### iterations, but we still allow some uncertainty around the mean
  ### value for all Phi's.
  sigma.Phi.Curr <-
    sqrt(1 / rgamma(1, shape = (n.G - 1) / 2,
                       rate = sum((log.Phi.Curr - log.Phi.Obs.mean)^2) / 2))
  # nu.Phi.Curr <- log.Phi.Obs.mean + rnorm(1) * sigma.Phi.Curr / sqrt(n.G)
  nu.Phi.Curr <- 0.0
  .cubfitsEnv$my.print(sigma.Phi.Curr)

  # Update prior's acceptance and adaptive.
  accept <- 1
  my.update.acceptance("p", accept)
  my.update.adaptive("p", accept)

  # Only nu.Phi and sigma.Phi are used.
  ret <- list(nu.Phi = nu.Phi.Curr, sigma.Phi = sigma.Phi.Curr)
  ret
} # End of my.pPropTypeNoObs.lognormal_MG0().

my.pPropTypeNoObs.fixed_SM <- function(n.G, log.Phi.Obs.mean, log.Phi.Curr,
    nu.Phi.Curr, sigma.Phi.Curr,
    p.DrawScale = 1, p.DrawScale.prev = 1, Phi.Curr = NULL){
  # Do nothing but skip this step.

  # Update prior's acceptance and adaptive.
  accept <- 1
  my.update.acceptance("p", accept)
  my.update.adaptive("p", accept)

  ret <- list(nu.Phi = nu.Phi.Curr, sigma.Phi = sigma.Phi.Curr)
  ret
} # End of my.pPropTypeNoObs.fixed_SM().
