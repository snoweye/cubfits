### Gibbs sampler and M-H for step 2, mainly for hyperparameters.

### Draw inv-/gamma for lognormal priors (nu.Phi, sigma.Phi) around current
### mean of log expression.
my.pPropTypeNoObs.lognormal_MH <- function(n.G, phi.Curr,
    p.Curr, hp.param, p.DrawScale = 0.1, p.DrawScale.prev = 0.1){
  ### Dispatch.
  nu.Phi.Curr <- p.Curr[1]
  sigma.Phi.Curr <- p.Curr[2]
  # hp.sigma.Phi <- hp.param$hp.sigma.Phi

  ### Propose sigma.Phi.Curr.
  proplist <- my.propose.sigma.Phi.Gamma(sigma.Phi.Curr) # , hp.sigma.Phi)

  ### M-H step.
  list.Curr <- list(nu.Phi = nu.Phi.Curr, sigma.Phi = sigma.Phi.Curr)
  ret <- my.draw.lognormal.hp.MH(proplist, list.Curr, phi.Curr)

  ### Update prior's acceptance and adaptive.
  my.update.acceptance("p", ret$accept)
  my.update.adaptive("p", ret$accept)

  ### Only nu.Phi and sigma.Phi are used.
  ret <- c(ret$nu.Phi, ret$sigma.Phi)
  ret
} # my.pPropTypeNoObs.lognormal_MH().

