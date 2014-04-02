# This either returns a log posterior vector or sum of the vector.
#
# These functions are for all genes.

# Get the specific function according to the options.
get.my.logPosteriorAll <- function(model.Phi){
  if(!any(model.Phi[1] %in% .CF.CT$model.Phi)){
    stop("model.Phi is not found.")
  }
  ret <- eval(parse(text = paste("my.logPosteriorAll.",
                                 model.Phi[1], sep = "")))
  assign("my.logPosteriorAll", ret, envir = .cubfitsEnv)
  ret
} # End of get.my.logPosteriorAll().


# Function to calculate complete log-posterior for
# (phi, b, sigmaWsq, mu.Phi, sigma.Phi.sq) given y, n, and phi.Obs
my.logPosteriorAll.lognormal <- function(phi, phi.Obs, y, n, b, sigmaWsq,
    mu.Phi = 0, sigma.Phi.sq = 1, reu13.df = NULL){
  ret <- dlnorm(phi.Obs, log(phi), sqrt(sigmaWsq), log = TRUE) +
         my.logdmultinomCodAllR(b, phi, y, n, reu13.df = reu13.df) +
         dlnorm(phi, mu.Phi, sqrt(sigma.Phi.sq), log = TRUE)

  ### Note that prior for all Phi is lognormal(mu.Phi, sigma.Phi), and assume
  ### proper priors for mu.Phi (normal) and sigma.Phi (inverse gamma.)
  ### It is necessary to add this for proportional posterior probability.
  ### No need to add this for acceptance ratio in our random walk, since
  ### the acceptance ratio cancels out this term, and we assume non-informative
  ### prior for mu.Phi and sigma.Phi.sq.
  ### i.e. p(mu.Phi, sigma.Phi.sq) \propto 1/sigma.Phi.sq.
  ### See Gelman et al. (2003), p.75 for details.

  # ret <- ret +
  #        dnorm(mu.Phi, some.mean, some.sd, log = TRUE) +
  #        dgamma(1/sigma.Phi.sq, some.alpha, some.beta, log = TRUE) + sigma.Phi.sq^2

  ret
} # End of my.logPosteriorAll.lognormal().

