### Provide number of coefficients for the given model.

### Get the constant according to the options.
get.my.ncoef <- function(model){
  if(model[1] == "rocnsef"){
    ret <- 3
  } else if(model[1] == "roc"){
    ret <- 2
  } else if(model[1] == "nsef"){
    ret <- 2
  } else{
    stop("model is not found.")
  }

  assign("my.ncoef", ret, envir = .cubfitsEnv)
  ret
} # End of get.my.ncoef().

get.my.coefnames <- function(model){
  if(model[1] == "rocnsef"){
    ret <- c("log.mu", "Delta.t", "omega")
  } else if(model[1] == "roc"){
    ret <- c("log.mu", "Delta.t")
  } else if(model[1] == "nsef"){
    ret <- c("log.mu", "omega")
  } else{
    stop("model is not found.")
  }

  assign("my.coefnames", ret, envir = .cubfitsEnv)
  ret
} # End of get.my.coefnames().
