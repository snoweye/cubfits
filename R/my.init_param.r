my.pInit <- function(p.Init, phi.Obs, model.Phi, p.nclass = 2,
    cub.method = c("fits", "appr", "pred")){
  if(is.null(p.Init)){
    sigmasqObs <- var(log(phi.Obs))
    sigmaW.Init <- sqrt(.2 * sigmasqObs)        # \sigma_W

    if(model.Phi == "lognormal"){
      nu.Phi.Init <- mean(log(phi.Obs))
      sigma.Phi.Init <- sqrt(.8 * sigmasqObs)     # \sigma_phi

      ret <- c(sigmaW.Init, nu.Phi.Init, sigma.Phi.Init)
    } else if(model.Phi[1] == "logmixture"){
      if(p.nclass <= 1){
        stop("p.nclass > 1")
      }
      log.phi.Obs <- matrix(log(phi.Obs), ncol = 1)
      tmp <- EMCluster::rand.EM(log.phi.Obs, nclass = as.integer(p.nclass),
                                min.n = length(phi.Obs) * 0.025,
                                EMC = EMCluster::.EMControl(short.eps = Inf))
      tmp.em <- EMCluster::rand.EM(log.phi.Obs, nclass = as.integer(p.nclass),
                                   min.n = length(phi.Obs) * 0.025,
                                   EMC = EMCluster::.EMControl())
      if(tmp.em$llhdval > tmp$llhdval){
        tmp <- tmp.em
      } 
      id <- order(tmp$Mu)
      paramlog <- c(tmp$pi[id], tmp$Mu[id], sqrt(tmp$LTSigma[id]))

      ret <- c(sigmaW.Init, paramlog)
    } else{
      stop("model.Phi is not found.")
    }

    if(cub.method[1] == "appr"){
      ret <- ret[-1]
    }
  } else{
    ret <- p.Init
  }
  
  ret
} # End of my.pInit().
