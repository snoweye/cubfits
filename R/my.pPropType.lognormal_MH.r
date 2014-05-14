### Gibbs sampler and M-H for step 2, mainly for hyperparameters.

### Draw Gibbs Sampler given current status for measure error (sigmaW) and
### drew inv-/gamma for lognormal priors (nu.Phi, sigma.Phi) around current
### mean of log expression.
my.pPropType.lognormal_MH <- function(n.G, log.phi.Obs, phi.Curr,
    p.Curr, hp.param, p.DrawScale = 0.1, p.DrawScale.prev = 0.1){
  ### Dispatch.
  nu.Phi.Curr <- p.Curr[2]
  sigma.Phi.Curr <- p.Curr[3]
  log.phi.Curr <- log(phi.Curr)
  # hp.sigma.Phi <- hp.param$hp.sigma.Phi

  ### Draw \sigma^{2*}_W from IG((n_G - 1) / 2,
  ###                            (n_G - 1) S^{2(t)}_{phi_{obs}} / 2)
  sigmaW.Curr <-
    sqrt(1 / rgamma(1, shape = (n.G - 1) / 2,
                       rate = sum((log.phi.Obs - log.phi.Curr)^2) / 2))

  ### Propose sigma.Phi.Curr.
  proplist <- my.propose.sigma.Phi.Gamma(sigma.Phi.Curr) # , hp.sigma.Phi)

  ### M-H step.
  list.Curr <- list(nu.Phi = nu.Phi.Curr, sigma.Phi = sigma.Phi.Curr)
  ret <- my.draw.lognormal.hp.MH(proplist, list.Curr, phi.Curr)

  ### Update prior's acceptance and adaptive.
  my.update.acceptance("p", ret$accept)
  my.update.adaptive("p", ret$accept)

  ### Only nu.Phi and sigma.Phi are used.
  ret <- c(sigmaW.Curr, ret$nu.Phi, ret$sigma.Phi)
  ret
} # my.pPropType.lognormal_MH().


### Propose sigma.Phi in log scale via gamma distribution.
my.propose.sigma.Phi.Gamma <- function(sigma.Phi.Curr, # hp.sigma.Phi,
    hp.gamma.shape = .CF.PARAM$hp.gamma.shape,
    hp.gamma.scale = .CF.PARAM$hp.gamma.scale){
  ### Draw from proposal.
  ### Gamma distribution with a flat prior.
  # sigma2.Phi.New <- rgamma(1, 1 / hp.sigma.Phi, scale = hp.sigma.Phi)
  ### Gamma distribution with an informative  prior.
  sigma2.Phi.New <- rgamma(1, hp.gamma.shape, scale = hp.gamma.scale)
  ### Inverse Gamma distribution with an informative prior.
  # sigma2.Phi.New <- 1 / rgamma(1, hp.invgamma.shape,
  #                                 scale = hp.invgamma.scale)

  ### Adjust accroding to the restriction.
  sigma.Phi.New <- sqrt(sigma2.Phi.New)
  nu.Phi.New <- -sigma2.Phi.New / 2

  ### Compute log ratio of prior since lognormal is not symmetric.
  ### Gamma distribution with a flat prior.
  # lir <- dgamma(sigma.Phi.New, 1 / hp.sigma.Phi,
  #               scale = hp.sigma.Phi, log = TRUE) -
  #        dgamma(sigma.Phi.Curr, 1 / hp.sigma.Phi,
  #               scale = hp.sigma.Phi, log = TRUE)
  ### Gamma distribution.
  lir <- dgamma(sigma.Phi.New, hp.gamma.shape,
                scale = hp.gamma.scale, log = TRUE) -
         dgamma(sigma.Phi.Curr, hp.gamma.shape,
                scale = hp.gamma.scale, log = TRUE)
  ### Inverse Gamma distribution with an informative prior.
  ### Y = 1 / X, |J| = 1 / X^2, f_Y(y) = f_X(x = 1/y) |J|
  # lir <- dgamma(1 / sigma.Phi.New, hp.gamma.shape,
  #               scale = hp.invgamma.scale, log = TRUE) -
  #        2 * log(sigma.Phi.New) -
  #        dgamma(1 / sigma.Phi.Curr, hp.gamma.shape,
  #               scale = hp.invgamma.scale, log = TRUE) +
  #        2 * log(sigma.Phi.Curr)

  ### Return.
  ret <- list(nu.Phi = as.numeric(nu.Phi.New),
              sigma.Phi = as.numeric(sigma.Phi.New),
              lir = lir)
  ret
} # End of my.propose.sigma.Phi.Gamma().

### Do the M-H.
my.draw.lognormal.hp.MH <- function(proplist, list.Curr, phi.Curr){
  ### Compute probability ratio.
  lpr <- sum(dlnorm(phi.Curr, meanlog = proplist$nu.Phi,
                    sdlog = proplist$sigma.Phi, log = TRUE)) -
         sum(dlnorm(phi.Curr, meanlog = list.Curr$nu.Phi,
                    sdlog = list.Curr$sigma.Phi, log = TRUE))

  ### log Acceptance probability.
  logAcceptProb <- lpr - proplist$lir

  ### Error handling -- interpreting NaN etc. as ~= 0.
  if(!is.finite(logAcceptProb)){
    warning("log acceptance probability not finite in hyperparam draw")
    logAcceptProb <- -Inf
  }

  ### Run MH acceptance rule.
  if(-rexp(1) < logAcceptProb){
    ret <- proplist
    ret$accept <- 1
  } else{
    ret <- list.Curr
    ret$accept <- 0
  }

  ### Return.
  ret
} # End of my.draw.lognormal.hp.MH().

