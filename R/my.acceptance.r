### Inital global storages for acceptance rate.
my.set.acceptance <- function(nSave, n.aa,
    n.p = NULL,
    n.G = NULL, n.G.pred = NULL){
  if(.CF.DP$trace.acceptance){
    .cubfitsEnv$acceptance <- list()

    ### For acceptance rate in S/M.
    .cubfitsEnv$acceptance$b <- rep(0L, n.aa)

    ### For acceptance rate in prior.
    .cubfitsEnv$acceptance$p <- rep(0L, n.p)

    ### For acceptance rate in training.
    if(!is.null(n.G)){
      .cubfitsEnv$acceptance$phi <- rep(0L, n.G)
    }

    ### For acceptance rate in prediction.
    if(!is.null(n.G.pred)){
      .cubfitsEnv$acceptance$phi.pred <- rep(0L, n.G.pred)
    }
  }

  invisible()
} # End of my.set.acceptance().


### Updating function to the global storages based on the variable name.
my.update.acceptance <- function(var.name, accept){
  if(.CF.DP$trace.acceptance){
    .cubfitsEnv$acceptance[[var.name]] <-
      .cubfitsEnv$acceptance[[var.name]] + accept
  }

  invisible()
} # End of my.update.acceptance().
