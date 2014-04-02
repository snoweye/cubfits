# This either returns a logL vector or sum of the vector.
#
# These functions are for all genes.

# Get the specific function according to the options.
get.my.logLAll <- function(model.Phi){
  if(!any(model.Phi[1] %in% .CF.CT$model.Phi)){
    stop("model.Phi is not found.")
  }
  ret <- eval(parse(text = paste("my.logLAll.",
                                 model.Phi[1], sep = "")))
  assign("my.logLAll", ret, envir = .cubfitsEnv)
  ret
} # End of get.my.logLAll().


# Function to calculate complete logL for
# (phi, b, sigmaWsq) given y, n, and phi.Obs
my.logLAll.lognormal <- function(phi, phi.Obs, y, n, b, sigmaWsq, reu13.df = NULL){
  dlnorm(phi.Obs, log(phi), sqrt(sigmaWsq), log = TRUE) +
  my.logdmultinomCodAllR(b, phi, y, n, reu13.df = reu13.df)
} # End of my.logLAll.lognormal().

