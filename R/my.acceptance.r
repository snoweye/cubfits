# Inital global storages for acceptance rate.
my.set.acceptance <- function(nSave, n.aa, G = NULL, G.pred = NULL){
  if(.CF.DP$trace.acceptance){
    .cubfitsEnv$acceptance <- list(B = list(), phi = list(), phi.pred = list())

    # For acceptance rate in S/M.
    for(i.aa in 1:n.aa){
      .cubfitsEnv$acceptance$B[[i.aa]] <- 0L 
    }

    # For acceptance rate in training.
    if(!is.null(G)){
      .cubfitsEnv$acceptance$phi <- list(rep(0L, G))
    }

    # For acceptance rate in prediction.
    if(!is.null(G.pred)){
      .cubfitsEnv$acceptance$phi.pred <- list(rep(0L, G.pred))
    }
  }

  invisible()
} # End of my.set.acceptance().


# Updating function to the global storages based on the variable name.
my.update.acceptance <- function(var.name, accept, i = 1){
  if(.CF.DP$trace.acceptance){
    .cubfitsEnv$acceptance[[var.name]][[i]] <-
      .cubfitsEnv$acceptance[[var.name]][[i]] + accept
  }

  invisible()
} # End of my.update.acceptance().

