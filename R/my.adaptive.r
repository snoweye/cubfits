# These functions are for adaptive MCMC, and only for E[Phi | Phi^{obs}] and
# E[Phi]. See and use the acceptance function for B (S/M).

# Initial global storages.
my.set.adaptive <- function(nSave, G = NULL, phi.DrawScale = NULL,
    G.pred = NULL, phi.DrawScale.pred = NULL,
    renew.iter = .CF.AC$renew.iter, adaptive = .CF.CT$adaptive[1]){
  .cubfitsEnv$curr.renew <- 1
  .cubfitsEnv$adaptive <- list(phi = list(), phi.pred = list())
  .cubfitsEnv$DrawScale <- list(phi = list(), phi.pred = list())

  if(adaptive == "none"){
    total.renew <- 1
  } else{
    total.renew <- ceiling(nSave / renew.iter)
  }

  # Initials for updating acceptance rate.
  # Initial adaptive rate for phi.DrawScale, and phi.DrawScale.pred.
  if(!is.null(G)){
    for(i.renew in 1:total.renew){
      .cubfitsEnv$adaptive$phi[[i.renew]] <- rep(0L, G)

      if(length(phi.DrawScale) == 1){
        .cubfitsEnv$DrawScale$phi[[i.renew]] <- rep(phi.DrawScale, G)
      } else if(length(phi.DrawScale) == G){
        .cubfitsEnv$DrawScale$phi[[i.renew]] <- phi.DrawScale
      } else{
        stop("length of phi.DrawScale is incorrect.")
      }
    }
  }

  # For adaptive rate in prediction.
  if(!is.null(G.pred)){
    for(i.renew in 1:total.renew){
      .cubfitsEnv$adaptive$phi.pred[[i.renew]] <- rep(0L, G.pred)

      if(length(phi.DrawScale.pred) == 1){
        .cubfitsEnv$DrawScale$phi.pred[[i.renew]] <- rep(phi.DrawScale.pred, G.pred)
      } else if(length(phi.DrawScale.pred) == G.pred){
        .cubfitsEnv$DrawScale$phi.pred[[i.renew]] <- phi.DrawScale.pred
      } else{
        stop("length of phi.DrawScale.pred is incorrect.")
      }
    }
  }

  invisible()
} # End of my.set.adaptive().


# Updating function based on variable name and current iteration.
my.update.adaptive <- function(var.name, accept){
  .cubfitsEnv$adaptive[[var.name]][[.cubfitsEnv$curr.renew]] <-
    .cubfitsEnv$adaptive[[var.name]][[.cubfitsEnv$curr.renew]] + accept
  invisible()
} # End of my.update.adaptive().

