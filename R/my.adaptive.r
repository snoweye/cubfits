### These functions are for adaptive MCMC, and only for E[b],
### E[Phi | Phi^{obs}], and E[Phi].

### Initial global storages.
my.set.adaptive <- function(nSave,
    n.aa = NULL, b.DrawScale = NULL,
    n.p = NULL, p.DrawScale = NULL,
    n.G = NULL, phi.DrawScale = NULL,
    n.G.pred = NULL, phi.DrawScale.pred = NULL,
    renew.iter = .CF.AC$renew.iter, adaptive = .CF.CT$adaptive[1]){
  ### Check.
  if(adaptive == "none"){
    total.renew <- 1
  } else{
    total.renew <- ceiling(nSave / renew.iter)
  }

  ### Initials for updating acceptance rate.
  .cubfitsEnv$curr.renew <- 1
  .cubfitsEnv$adaptive <- list(b = list(),
                               p = list(),
                               phi = list(), phi.pred = list())
  .cubfitsEnv$DrawScale <- list(b = list(),
                                p = list(),
                                phi = list(), phi.pred = list())

  ### For adaptive rate in parameters.
  if(!is.null(n.aa)){
    for(i.renew in 1:total.renew){
      .cubfitsEnv$adaptive$b[[i.renew]] <- rep(0L, n.aa)
      if(length(b.DrawScale) == 1){
        .cubfitsEnv$DrawScale$b[[i.renew]] <- rep(b.DrawScale, n.aa)
      } else if(length(b.DrawScale) == n.aa){
        .cubfitsEnv$DrawScale$b[[i.renew]] <- b.DrawScale
      } else{
        stop("length of b.DrawScale is incorrect.")
      }
    }
  }

  ### For adaptive rate in prior
  if(!is.null(n.p)){
    for(i.renew in 1:total.renew){
      .cubfitsEnv$adaptive$p[[i.renew]] <- rep(0L, n.p)
      if(length(p.DrawScale) == 1){
        .cubfitsEnv$DrawScale$p[[i.renew]] <- rep(p.DrawScale, n.p)
      } else if(length(p.DrawScale) == n.p){
        .cubfitsEnv$DrawScale$p[[i.renew]] <- p.DrawScale
      } else{
        stop("length of p.DrawScale is incorrect.")
      }
    }
  }

  ### For adaptive rate in expectations.
  if(!is.null(n.G)){
    for(i.renew in 1:total.renew){
      .cubfitsEnv$adaptive$phi[[i.renew]] <- rep(0L, n.G)
      if(length(phi.DrawScale) == 1){
        .cubfitsEnv$DrawScale$phi[[i.renew]] <- rep(phi.DrawScale, n.G)
      } else if(length(phi.DrawScale) == n.G){
        .cubfitsEnv$DrawScale$phi[[i.renew]] <- phi.DrawScale
      } else{
        stop("length of phi.DrawScale is incorrect.")
      }
    }
  }

  ### For adaptive rate in predictions.
  if(!is.null(n.G.pred)){
    for(i.renew in 1:total.renew){
      .cubfitsEnv$adaptive$phi.pred[[i.renew]] <- rep(0L, n.G.pred)

      if(length(phi.DrawScale.pred) == 1){
        .cubfitsEnv$DrawScale$phi.pred[[i.renew]] <- rep(phi.DrawScale.pred,
                                                         n.G.pred)
      } else if(length(phi.DrawScale.pred) == n.G.pred){
        .cubfitsEnv$DrawScale$phi.pred[[i.renew]] <- phi.DrawScale.pred
      } else{
        stop("length of phi.DrawScale.pred is incorrect.")
      }
    }
  }

  invisible()
} # End of my.set.adaptive().


### Updating function based on variable name and current iteration.
my.update.adaptive <- function(var.name, accept){
  .cubfitsEnv$adaptive[[var.name]][[.cubfitsEnv$curr.renew]] <-
    .cubfitsEnv$adaptive[[var.name]][[.cubfitsEnv$curr.renew]] + accept
  invisible()
} # End of my.update.adaptive().
