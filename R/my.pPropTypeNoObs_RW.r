# Gibbs sampler and M-H for step 2, mainly for hyperparameters.

# Drew Gibbs Sampler given current status for measure error (sigmaW) and
# drew random walk for lognormal priors (nu.Phi, sigma.Phi) around current
# mean of log expression.
my.pPropTypeNoObs.lognormal_RW <- function(n.G, log.Phi.Obs.mean, log.Phi.Curr,
    nu.Phi.Curr, sigma.Phi.Curr,
    p.DrawScale = 1, p.DrawScale.prev = 1,
    Phi.Curr = NULL){
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
  ret <- list(nu.Phi = ret$nu.Phi, sigma.Phi = ret$sigma.Phi)
  ret
} # my.pPropTypeNoObs.lognormal_RW().
