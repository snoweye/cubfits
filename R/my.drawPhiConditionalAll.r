### Function to run MH draw for phi given observed (with log-Normal error) phi.Obs,
### count y, number of codons n, coefficients b, and measurement error
### log-variance sigmaWsq.
###
### Performs vectorized draw.
### Assumes phi.Curr1, phi.Obs, y, & n are vectors of length # of genes;
### b is 2-dim vector; sigmaWsq is scalar.
###
### Returns list with elements:
###   * phi.New  :   new value for phi resulting from MH draw
###   * accept:   boolean vector indicating if each proposal was accepted

my.drawPhiConditionalAll <- function(phi.Curr1, phi.Obs, y, n, b,
    p.Curr, phi.DrawScale = 1, phi.DrawScale.prev = 1, reu13.df = NULL){
  ### Propose new phi.
  prop <- .cubfitsEnv$my.proposePhiAll(phi.Curr1, phi.DrawScale = phi.DrawScale,
                                       phi.DrawScale.prev = phi.DrawScale.prev)

  ### Calculate acceptance prob.
  lpCurr <- .cubfitsEnv$my.logPosteriorAll(
              phi.Curr1, phi.Obs, y, n, b, p.Curr, reu13.df = reu13.df)
  lpProp <- .cubfitsEnv$my.logPosteriorAll(
              prop$phi.Prop, phi.Obs, y, n, b, p.Curr, reu13.df = reu13.df)
  logAcceptProb <- lpProp - lpCurr - prop$lir

  ### Run MH acceptance rule.
  u <- runif(length(phi.Curr1))
  accept <- u < exp(logAcceptProb)
  phi.New <- phi.Curr1
  phi.New[accept] <- prop$phi.Prop[accept]

  ### Extra update and trace.
  my.update.acceptance("phi", accept)
  my.update.adaptive("phi", accept)
    
  ### Return.
  ret <- phi.New
  ret
} # End of my.drawPhiConditionalAll().

