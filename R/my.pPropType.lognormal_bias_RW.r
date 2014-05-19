### Gibbs sampler and random walk for step 2, mainly for hyperparameters.

### Draw Gibbs Sampler given current status for measure error (sigmaW) and
### drew inv-/gamma for lognormal priors (nu.Phi, sigma.Phi) around current
### mean of log expression. Also, draw bias.Phi.
my.pPropType.lognormal_bias_RW <- function(n.G, log.phi.Obs, phi.Curr,
    p.Curr, hp.param){
  ### Dispatch.
  nu.Phi.Curr <- p.Curr[2]
  sigma.Phi.Curr <- p.Curr[3]
  bias.Phi.Curr <- p.Curr[4]
  log.phi.Curr <- log(phi.Curr)
  p.DrawScale <- .cubfitsEnv$all.DrawScale$p
  p.DrawScale.prev <- .cubfitsEnv$all.DrawScale$p.prev

  ### Draw \sigma^{2*}_W from IG((n_G - 1) / 2,
  ###                            (n_G - 1) S^{2(t)}_{phi_{obs}} / 2)
  sigmaW.Curr <-
    sqrt(1 / rgamma(1, shape = (n.G - 1) / 2,
                       rate = sum((log.phi.Obs - log.phi.Curr)^2) / 2))

  ### Propose sigma.Phi.Curr.
  proplist <- my.propose.sigma.Phi.RW(sigma.Phi.Curr, p.DrawScale[1],
                                      p.DrawScale.prev[1])

  ### M-H step for hyperparameters.
  list.Curr <- list(nu.Phi = nu.Phi.Curr, sigma.Phi = sigma.Phi.Curr)
  ret <- my.draw.lognormal.hp.MH(proplist, list.Curr, phi.Curr)

  ### Propose bias.Phi.
  proplist <- my.propose.bias.Phi.RW(bias.Phi.Curr, p.DrawScale[2],
                                     p.DrawScale.prev[2])

  ### M-H step for hyperparameters of bias.
  list.Curr <- list(bias.Phi = bias.Phi.Curr,
                    nu.Phi = ret$nu.Phi, sigma.Phi = ret$sigma.Phi)
  ret.bias <- my.draw.lognormal_bias.hp.MH(proplist, list.Curr, log.phi.Obs,
                                           phi.Curr)

  ### Update prior's acceptance and adaptive.
  accept <- c(ret$accept, ret.bias$accept)
  my.update.acceptance("p", accept)
  my.update.adaptive("p", accept)

  ### Only nu.Phi and sigma.Phi are used.
  ret <- c(sigmaW.Curr, ret$nu.Phi, ret$sigma.Phi, ret.bias$bias.Phi)
  ret
} # my.pPropType.lognormal_bias_RW().


### Propose bias.Phi in log scale via random walk.
my.propose.bias.Phi.RW <- function(bias.Phi.Curr,
    bias.Phi.DrawScale = 0.1, bias.Phi.DrawScale.prev = 0.1){
  ### Draw from proposcal.
  bias.Phi.New <- rnorm(1, mean = log(bias.Phi.Curr),
                           sd = bias.Phi.DrawScale)

  ### Compute log ratio of prior since lognormal is not symmetric.
  lir <- dlnorm(exp(bias.Phi.New), meanlog = bias.Phi.Curr,
                sdlog = bias.Phi.DrawScale, log = TRUE) -
         dlnorm(exp(bias.Phi.Curr), meanlog = bias.Phi.New,
                sdlog = bias.Phi.DrawScale.prev, log = TRUE)

  ### Return.
  ret <- list(bias.Phi = as.numeric(bias.Phi.New),
              lir = lir)
  ret
} # End of my.propose.bias.Phi.RW().

### Do the M-H.
my.draw.lognormal_bias.hp.MH <- function(proplist, list.Curr, log.phi.Obs,
    phi.Curr){
  ### Compute probability ratio.
  lpr <- sum(dlnorm(log.phi.Obs,
                    meanlog = list.Curr$nu.Phi + proplist$bias.Phi,
                    sdlog = list.Curr$sigma.Phi, log = TRUE)) -
         sum(dlnorm(log.phi.Obs,
                    meanlog = list.Curr$nu.Phi + list.Curr$bias.Phi,
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
} # End of my.draw.lognormal_bias.hp.MH().
