### Since there is no phi, it is no sense to implement the function.
### However, need a fake function to avoid polymorphism error.

my.pPropTypeNoObs.lognormal_bias <- function(n.G, phi.Curr,
    p.Curr, hp.param){
  .cubfitsEnv$my.stop("It is no sense to estimate bias without phi.")
} # my.pPropTypeNoObs.lognormal_bias().
