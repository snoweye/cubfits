# Gibbs sampler and M-H for step 2, mainly for hyperparameters.

# Drew Gibbs Sampler given current status for measure error (sigmaW) and
# drew random walk for lognormal priors (nu.Phi, sigma.Phi) around current
# mean of log expression.
my.pPropType.lognormal_RW <- function(n.G, log.Phi.Obs, log.Phi.Curr,
    nu.Phi.Curr, sigma.Phi.Curr, log.Phi.Obs.mean = 0,
    p.DrawScale = 1, p.DrawScale.prev = 1, Phi.Curr = NULL){
  # Draw \sigma^{2*}_W from IG((n_G - 1) / 2,
  #                            (n_G - 1) S^{2(t)}_{phi_{obs}} / 2)
  sigmaW.Curr <-
    sqrt(1 / rgamma(1, shape = (n.G - 1) / 2,
                       rate = sum((log.Phi.Obs - log.Phi.Curr)^2) / 2))

  # Propose sigma.Phi.Curr.
  proplist <- my.propose.sigma.Phi.RW(sigma.Phi.Curr, p.DrawScale[1],
                                      p.DrawScale.prev[1])

  # M-H step.
  p.Curr <- list(nu.Phi = nu.Phi.Curr, sigma.Phi = sigma.Phi.Curr)
  ret <- my.drawHyperparam.MH(proplist, p.Curr, Phi.Curr)

  # Update prior's acceptance and adaptive.
  accept <- ret$accept
  my.update.acceptance("p", accept)
  my.update.adaptive("p", accept)

  # Only nu.Phi and sigma.Phi are used.
  ret <- list(sigmaW = sigmaW.Curr,
              nu.Phi = ret$nu.Phi, sigma.Phi = ret$sigma.Phi)
  ret
} # my.pPropType.lognormal_RW().


# Propose sigma.Phi in log scale via random walk.
my.propose.sigma.Phi.RW <- function(sigma.Phi,
    sigma.Phi.DrawScale = 1, sigma.Phi.DrawScale.prev = 1){
  # Draw from proposal.
  sigma.Phi.New <- exp(rnorm(1, mean = log(sigma.Phi),
                                sd = sigma.Phi.DrawScale))
  nu.Phi.New <- -sigma.Phi.New^2 / 2

  # Compute log ratio of prior since lognormal is not symmetric.
  lir <- dlnorm(sigma.Phi.New, meanlog = log(sigma.Phi),
                sdlog = sigma.Phi.DrawScale, log = TRUE) -
         dlnorm(sigma.Phi, meanlog = log(sigma.Phi.New),
                sdlog = sigma.Phi.DrawScale.prev, log = TRUE)

  # Return.
  ret <- list(nu.Phi = as.numeric(nu.Phi.New),
              sigma.Phi = as.numeric(sigma.Phi.New),
              lir = lir)
  ret
} # End of my.propose.sigma.Phi.RW().


# Do the M-H.
my.drawHyperparam.MH <- function(proplist, p.Curr, Phi.Curr){
  # Compute probability ratio
  lpr <- sum(dlnorm(Phi.Curr, meanlog = proplist$nu.Phi,
                    sdlog = proplist$sigma.Phi, log = TRUE)) -
         sum(dlnorm(Phi.Curr, meanlog = p.Curr$nu.Phi,
                    sdlog = p.Curr$sigma.Phi, log = TRUE))

  # log Acceptance probability.
  logAcceptProb <- lpr - proplist$lir

  # Error handling -- interpreting NaN etc. as ~= 0.
  if(!is.finite(logAcceptProb)){
    warning("log acceptance probability not finite in hyperparam draw")
    logAcceptProb <- -Inf
  }

  # Run MH acceptance rule.
  if(-rexp(1) < logAcceptProb){
    ret <- proplist
    ret$accept <- 1
  } else{
    ret <- p.Curr
    ret$accept <- 0
  }

  # Return.
  ret
} # End of my.drawHyperparam.MH().

