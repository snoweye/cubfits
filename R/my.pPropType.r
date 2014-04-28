# Gibbs sampler for step 2, mainly for hyperparameters.

# Get the specific function according to the options.
get.my.pPropType <- function(type){
  if(!any(type[1] %in% .CF.CT$type.p)){
    stop("type is not found.")
  }
  ret <- eval(parse(text = paste("my.pPropType.", type[1], sep = "")))
  assign("my.pPropType", ret, envir = .cubfitsEnv)
  ret
} # End of get.my.pPropType().


# Drew Gibbs Sampler given current status for lognormal prior around fixed
# mean of log expression.
my.pPropType.lognormal_fix <- function(n.G, log.Phi.Obs, log.Phi.Curr,
    nu.Phi.Curr, sigma.Phi.Curr, log.Phi.Obs.mean = 0,
    p.DrawScale = 1, p.DrawScale.prev = 1, Phi.Curr = NULL){
  # sigma.Phi.Curr is unused in this case.

  # Draw \sigma^{2*}_W from IG((n_G - 1) / 2,
  #                            (n_G - 1) S^{2(t)}_{phi_{obs}} / 2)
  sigmaW.Curr <-
    sqrt(1 / rgamma(1, shape = (n.G - 1) / 2,
                       rate = sum((log.Phi.Obs - log.Phi.Curr)^2) / 2))

  ### The next is assuming non-informative priors for mu.Phi and sigma.Phi.sq.
  ### i.e. p(mu.Phi, sigma.Phi.sq) \propto 1/sigma.Phi.sq.
  ### See Gelman et al. (2003), p.75 for details.

  ### As my.pPropTypeNoObs.lognormal_fix(), we can set a fixed mean to log.Phi.Obs
  ### if phi.Obs has been normalized to mean=1 prior.
  # Draw \sigma^{2*}_phi from IG((n_G - 1) / 2, (n_G - 1) S^{2(t)}_{phi} / 2)
  sigma.Phi.Curr <-
    sqrt(1 / rgamma(1, shape = (n.G - 1) / 2,
                       rate = sum((log.Phi.Curr - log.Phi.Obs.mean)^2) / 2))
  # Draw \mu^*_phi from N(sum_g log(phi^{(t)}_g) / n_G, \sigma^{2*}_phi / n_G)
  nu.Phi.Curr <- log.Phi.Obs.mean + rnorm(1) * sigma.Phi.Curr / sqrt(n.G)

  # Update prior's acceptance and adaptive.
  accept <- 1
  my.update.acceptance("p", accept)
  my.update.adaptive("p", accept)

  # Return.
  ret <- list(sigmaW = sigmaW.Curr,
              nu.Phi = nu.Phi.Curr, sigma.Phi = sigma.Phi.Curr)
  ret
} # my.pPropType.lognormal_fix().

# Drew Gibbs Sampler given current status for lognormal prior.
my.pPropType.lognormal <- function(n.G, log.Phi.Obs, log.Phi.Curr, nu.Phi.Curr,
    sigma.Phi.Curr, log.Phi.Obs.mean = 0,
    p.DrawScale = 1, p.DrawScale.prev = 1, Phi.Curr = NULL){
  # sigma.Phi.Curr is unused in this case.
  # log.Phi.Obs.mean is unused in this case.

  # Draw \sigma^{2*}_W from IG((n_G - 1) / 2, (n_G - 1) S^{2(t)}_{phi_{obs}} / 2)
  sigmaW.Curr <- sqrt(1 / rgamma(1, shape = (n.G - 1) / 2,
                                   rate = sum((log.Phi.Obs - log.Phi.Curr)^2) / 2))

  ### The next is assuming non-informative priors for mu.Phi and sigma.Phi.sq.
  ### i.e. p(mu.Phi, sigma.Phi.sq) \propto 1/sigma.Phi.sq.
  ### See Gelman et al. (2003), p.75 for details.

  # Draw \sigma^{2*}_phi from IG((n_G - 1) / 2, (n_G - 1) S^{2(t)}_{phi} / 2)
  sigma.Phi.Curr <-
    sqrt(1 / rgamma(1, shape = (n.G - 1) / 2,
                       rate = sum((log.Phi.Curr - nu.Phi.Curr)^2) / 2))
  # Draw \mu^*_phi from N(sum_g log(phi^{(t)}_g) / n_G, \sigma^{2*}_phi / n_G)
  nu.Phi.Curr <- mean(log.Phi.Curr) + rnorm(1) * sigma.Phi.Curr / sqrt(n.G)

  # Update prior's acceptance and adaptive.
  accept <- 1
  my.update.acceptance("p", accept)
  my.update.adaptive("p", accept)

  # Return.
  ret <- list(sigmaW = sigmaW.Curr,
              nu.Phi = nu.Phi.Curr, sigma.Phi = sigma.Phi.Curr)
  ret
} # End of my.pPropType.lognormal().

### This method violates MCMC fundamental assumption.
my.pPropType.lognormal_MG <- function(n.G, log.Phi.Obs, log.Phi.Curr,
    nu.Phi.Curr, sigma.Phi.Curr, log.Phi.Obs.mean = 0,
    p.DrawScale = 1, p.DrawScale.prev = 1, Phi.Curr = NULL){
  # sigma.Phi.Curr is unused in this case.

  # Draw \sigma^{2*}_W from IG((n_G - 1) / 2,
  #                            (n_G - 1) S^{2(t)}_{phi_{obs}} / 2)
  sigmaW.Curr <-
    sqrt(1 / rgamma(1, shape = (n.G - 1) / 2,
                       rate = sum((log.Phi.Obs - log.Phi.Curr)^2) / 2))

  ### The next is assuming non-informative priors for mu.Phi and sigma.Phi.sq.
  ### i.e. p(mu.Phi, sigma.Phi.sq) \propto 1/sigma.Phi.sq.
  ### See Gelman et al. (2003), p.75 for details.

  ### As my.pPropTypeNoObs.lognormal_fix(), we can set a fixed mean to log.Phi.Obs
  ### if phi.Obs has been normalized to mean=1 prior.
  # Draw \sigma^{2*}_phi from IG((n_G - 1) / 2, (n_G - 1) S^{2(t)}_{phi} / 2)
  sigma.Phi.Curr <-
    sqrt(1 / rgamma(1, shape = (n.G - 1) / 2,
                       rate = sum((log.Phi.Curr - log.Phi.Obs.mean)^2) / 2))
  # Draw \mu^*_phi from N(sum_g log(phi^{(t)}_g) / n_G, \sigma^{2*}_phi / n_G)
  # nu.Phi.Curr <- log.Phi.Obs.mean + rnorm(1) * sigma.Phi.Curr / sqrt(n_G)
  nu.Phi.Curr <- -sigma.Phi.Curr^2 / 2

  # Update prior's acceptance and adaptive.
  accept <- 1
  my.update.acceptance("p", accept)
  my.update.adaptive("p", accept)

  # Return.
  ret <- list(sigmaW = sigmaW.Curr,
              nu.Phi = nu.Phi.Curr, sigma.Phi = sigma.Phi.Curr)
  ret
} # my.pPropType.lognormal_MG().

### This method violates MCMC fundamental assumption.
my.pPropType.lognormal_MG0 <- function(n.G, log.Phi.Obs, log.Phi.Curr,
    nu.Phi.Curr, sigma.Phi.Curr, log.Phi.Obs.mean = 0,
    p.DrawScale = 1, p.DrawScale.prev = 1, Phi.Curr = NULL){
  # sigma.Phi.Curr is unused in this case.

  # Draw \sigma^{2*}_W from IG((n_G - 1) / 2,
  #                            (n_G - 1) S^{2(t)}_{phi_{obs}} / 2)
  sigmaW.Curr <-
    sqrt(1 / rgamma(1, shape = (n.G - 1) / 2,
                       rate = sum((log.Phi.Obs - log.Phi.Curr)^2) / 2))

  ### The next is assuming non-informative priors for mu.Phi and sigma.Phi.sq.
  ### i.e. p(mu.Phi, sigma.Phi.sq) \propto 1/sigma.Phi.sq.
  ### See Gelman et al. (2003), p.75 for details.

  ### As my.pPropTypeNoObs.lognormal_fix(), we can set a fixed mean to log.Phi.Obs
  ### if phi.Obs has been normalized to mean=1 prior.
  # Draw \sigma^{2*}_phi from IG((n_G - 1) / 2, (n_G - 1) S^{2(t)}_{phi} / 2)
  sigma.Phi.Curr <-
    sqrt(1 / rgamma(1, shape = (n.G - 1) / 2,
                       rate = sum((log.Phi.Curr - log.Phi.Obs.mean)^2) / 2))
  # Draw \mu^*_phi from N(sum_g log(phi^{(t)}_g) / n_G, \sigma^{2*}_phi / n_G)
  # nu.Phi.Curr <- log.Phi.Obs.mean + rnorm(1) * sigma.Phi.Curr / sqrt(n_G)
  nu.Phi.Curr <- 0.0

  # Update prior's acceptance and adaptive.
  accept <- 1
  my.update.acceptance("p", accept)
  my.update.adaptive("p", accept)

  # Return.
  ret <- list(sigmaW = sigmaW.Curr,
              nu.Phi = nu.Phi.Curr, sigma.Phi = sigma.Phi.Curr)
  ret
} # my.pPropType.lognormal_MG0().


# Do nothing but skipt the step. 
my.pPropType.fixed_SM <- function(n.G, log.Phi.Obs, log.Phi.Curr, nu.Phi.Curr,
    sigma.Phi.Curr, log.Phi.Obs.mean = 0,
    p.DrawScale = 1, p.DrawScale.prev = 1, Phi.Curr = NULL){
  # Draw \sigma^{2*}_W from IG((n_G - 1) / 2,
  #                            (n_G - 1) S^{2(t)}_{phi_{obs}} / 2)
  sigmaW.Curr <-
    sqrt(1 / rgamma(1, shape = (n.G - 1) / 2,
                       rate = sum((log.Phi.Obs - log.Phi.Curr)^2) / 2))

  # Update prior's acceptance and adaptive.
  accept <- 1
  my.update.acceptance("p", accept)
  my.update.adaptive("p", accept)

  # Return.
  ret <- list(sigmaW = sigmaW.Curr,
              nu.Phi = nu.Phi.Curr, sigma.Phi = sigma.Phi.Curr)
  ret
} # my.pPropType.fixed_SM().

