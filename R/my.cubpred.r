# Function to run MCMC inference for following model:
#   Given phi, n, b; across n.G genes:
#     phi      ~ rlnorm(n.G, mu.Phi, sigma.Phi)
#     phi.Obs   ~ rlnorm(n.G, log(phi), sigmaW)
#     y      ~ for(aa) rmultinom(n.G, invmlogit( phi * b[[aa]] ), n[[aa]] )
#     phi.pred ~ rlnorm(n.G.pred, mu.Phi, sigma.Phi)
#
# Expects phi.Obs as vector of length n.G.
# 
# Subsequent Gibbs sampler:
# (1) Sample of b | phi, y, using VGAM fit for Gaussian proposal
# (2) Sample sigmaW | phi, phi.Obs
#     Sample mu.Phi, sigma.Phi | x 
# (3) Sample of phi | b, phi.Obs, nu.Phi, bsig.Phi, m.Phi, ww, sigmaW, y,
#     using logNormal proposal 
# (4) Sample of phi.pred | b, nu.Phi, bsig.Phi, m.Phi, ww, y.pred,
#     using logNormal proposal 

### All genes have observations in the training set and
### no observation in the testing set.
my.cubpred <- function(reu13.df.obs, phi.Obs, y, n,
    reu13.df.pred, y.pred, n.pred,
    nIter = 1000, burnin = 100,
    pInit = NULL,
    bInit = NULL, init.b.Scale = 1, b.DrawScale = 1,
    phi.Init = NULL, init.Phi.Scale = 1, phi.DrawScale = 1,
    phi.Init.pred = NULL, phi.DrawScale.pred = 1,
    model = .CF.CT$model[1], model.Phi = .CF.CT$model.Phi[1],
    adaptive = .CF.CT$adaptive[1],
    scale.Phi = .CF.CT$scale.Phi[1],
    verbose = .CF.DP$verbose,
    iterThin = .CF.DP$iterThin, report = .CF.DP$report){

### Setup functions ###
  # Setup function pointers by type or model.
  my.function <- my.init.function(model = model[1], model.Phi = model.Phi[1],
                                  adaptive = adaptive[1])
  my.ncoef <- my.function$my.ncoef

### Check Data ###
  # Check phi.Obs is well-behaved.
  if(!(all(is.finite(phi.Obs)) && all(phi.Obs > 0))){
    .cubfitsEnv$my.stop("phi.Obs is invalid.")
  }
  if(abs(mean(phi.Obs) - 1) > 1e-8 && scale.Phi == "mean_one"){
    .cubfitsEnv$my.stop(paste("mean(phi.Obs) =", mean(phi.Obs)))
  }
  if(!is.null(phi.Init)){
    if(!(all(is.finite(phi.Init)) && all(phi.Init > 0))){
      .cubfitsEnv$my.stop("phi.Init is invalid.")
    }
    if(abs(mean(phi.Init) - 1) > 1e-8){
      .cubfitsEnv$my.stop(paste("mean(phi.Init) =", mean(phi.Init)))
    }
  }
  if(!is.null(phi.Init.pred)){
    if(!(all(is.finite(phi.Init.pred)) && all(phi.Init.pred > 0))){
      .cubfitsEnv$my.stop("phi.Init.pred is invalid.")
    }
    if(abs(mean(phi.Init.pred) - 1) > 1e-8 && scale.Phi == "mean_one"){
      .cubfitsEnv$my.stop(paste("mean(phi.Init.pred) =", mean(phi.Init.pred)))
    }
  }

  # Check if sort by ORF.
  my.check.rearrange(reu13.df.obs, y, n, phi.Obs = phi.Obs,
                     reu13.df.pred = reu13.df.pred, y.pred = y.pred,
                     n.pred = n.pred, phi.Init = phi.Init,
                     phi.Init.pred = phi.Init.pred)

### Initial Storages ###
  # Setup data structures for results.
  n.G <- length(phi.Obs)                    # # of genes
  n.aa <- length(y)                         # # of amino acids
  nsyns <- sapply(y, function(ybit){ dim(ybit)[2] })
                                            # # of synomous codons
  nBparams <- my.ncoef * sum(nsyns - 1)     # total # of regression parameters
  nSave <- (nIter + burnin) / iterThin + 1  # # of space for iterations

  # Storages for saving posterior samples.
  b.Mat <- my.generate.list(NA, nBparams, nSave)    # S/M
  p.Mat <- my.generate.list(NA, 3, nSave)           # hyperparameters
  phi.Mat <- my.generate.list(NA, n.G, nSave)       # E[Phi | Phi^{obs}]

  # Setup data structures for prediction.
  n.G.pred <- nrow(y.pred[[1]])                     # # of genes for prediction

  # Storages for saving posterior samples in the prediction set.
  phi.Mat.pred <- my.generate.list(NA, n.G.pred, nSave) # E[Phi]

### Initial Parameters ###
  # Initial values for b.
  bInitList <- .cubfitsEnv$my.fitMultinomAll(reu13.df.obs, phi.Obs, y, n)
  bRInitList <- lapply(bInitList, function(B){ B$R })
  if(is.null(bInit)){
    bInit <- lapply(bInitList,
               function(B){
                 B$coefficients +
                 init.b.Scale * backsolve(B$R, rnorm(nrow(B$R)))
               })
  } else{
    bInit <- lapply(bInit, function(B){ B$coefficients })
  }
  bInitVec <- unlist(bInit)

  # Initial values for p.
  if(is.null(pInit)){
    nu.Phi.Init <- mean(log(phi.Obs))
    sigmasqObs <- var(log(phi.Obs))
    bsig.Phi.Init <- sqrt(.8 * sigmasqObs)     # \sigma_phi
    sigmaWInit <- sqrt(.2 * sigmasqObs)        # \sigma_W
    pInit <- c(sigmaW = sigmaWInit, nu.Phi = nu.Phi.Init,
               bsig.Phi = bsig.Phi.Init)
  } else{
    sigmaWInit <- pInit["sigmaW"]
    nu.Phi.Init <- pInit["nu.Phi"]
    bsig.Phi.Init <- pInit["bsig.Phi"]
  }

  # Initial values for training phi.
  if(is.null(phi.Init)){
    phi.Init <- phi.Obs * exp(init.Phi.Scale * sigmaWInit * rnorm(n.G))
  }

  # Initial values for predicting phi.
  if(is.null(phi.Init.pred)){
    phi.Init.pred <- median(phi.Obs) *
                     exp(init.Phi.Scale * sigmaWInit * rnorm(n.G.pred))
    names(phi.pred.pred) <- rownames(y.pred[[1]])
  }

  # Set current step for b.
  b.Mat[[1]] <- bInitVec
  bCurr <- bInit

  # Set current step for p.
  p.Mat[[1]] <- pInit 
  pCurr <- pInit
  nu.Phi.Curr <- nu.Phi.Init
  bsig.Phi.Curr <- bsig.Phi.Init
  sigmaWCurr <- sigmaWInit

  # Set current step for phi.
  phi.Mat[[1]] <- phi.Init
  phi.Curr <- phi.Init
  log.Phi.Obs <- log(phi.Obs)

  # Set hyper-prior parameters
  log.Phi.Obs.mean <- mean(log.Phi.Obs)

  # Set current step for phi.pred.
  phi.Mat.pred[[1]] <- phi.Init.pred
  phi.Curr.pred <- phi.Init.pred

### MCMC here ###
  # Set acceptance rate storage.
  my.set.acceptance(nSave, n.aa, n.G = n.G, n.G.pred = n.G.pred)

  # Set adaptive storage.
  my.set.adaptive(nSave,
                  n.aa = n.aa, b.DrawScale = b.DrawScale,
                  n.G = n.G, phi.DrawScale = phi.DrawScale,
                  n.G.pred = n.G.pred, phi.DrawScale.pred = phi.DrawScale.pred,
                  adaptive = adaptive[1])
  b.DrawScale <- .cubfitsEnv$DrawScale$b[[1]]
  b.DrawScale.prev <- b.DrawScale
  phi.DrawScale <- .cubfitsEnv$DrawScale$phi[[1]]
  phi.DrawScale.pred <- .cubfitsEnv$DrawScale$phi.pred[[1]]
  phi.DrawScale.prev <- phi.DrawScale
  phi.DrawScale.pred.prev <- phi.DrawScale.pred

  # Run MCMC iterations.
  my.verbose(verbose, 0, report)
  .cubfitsEnv$my.dump(0, list = c("b.Mat", "p.Mat", "phi.Mat", "phi.Mat.pred"))

  # MCMC start.
  for(iter in 1:(nIter + burnin)){
    # Step 1: Update b using M-H step.
    bUpdates <- .cubfitsEnv$my.drawBConditionalAll(
                  bCurr, phi.Curr, y, n, reu13.df.obs,
                  bRInitList = bRInitList,
                  b.DrawScale = b.DrawScale,
                  b.DrawScale.prev = b.DrawScale.prev)
    bCurr  <- lapply(bUpdates, function(U){ U$bNew })

    # Step 2: Draw other parameters.
    log.Phi.Curr <- log(phi.Curr)
    pUpdate <- .cubfitsEnv$my.pPropType(n.G, log.Phi.Obs, log.Phi.Curr,
                                        nu.Phi.Curr, bsig.Phi.Curr,
                                        log.Phi.Obs.mean)
    sigmaWCurr <- pUpdate$sigmaWCurr 
    nu.Phi.Curr <- pUpdate$ppCurr$y0
    bsig.Phi.Curr <- pUpdate$ppCurr$bb
    pCurr <- c(sigmaWCurr, nu.Phi.Curr, bsig.Phi.Curr)
    mu.Phi.Curr <- nu.Phi.Curr
    sigma.Phi.sqCurr <- bsig.Phi.Curr^2

    # Step 3: Update phi using M-H step.
    phi.Update <- my.drawPhiConditionalAll(phi.Curr, phi.Obs, y, n,
                    bCurr, sigmaWCurr^2,
                    mu.Phi = mu.Phi.Curr,
                    sigma.Phi.sq = sigma.Phi.sqCurr, 
                    phi.DrawScale = phi.DrawScale,
                    phi.DrawScale.prev = phi.DrawScale.prev,
                    reu13.df = reu13.df.obs)
    phi.Curr <- phi.Update$phi.New

    # Step 4: Predict phi using M-H step. This is different to the Step 3.
    phi.pred <- my.drawPhiConditionalAllPred(phi.Curr.pred, y.pred, n.pred,
                  bCurr,
                  mu.Phi = mu.Phi.Curr,
                  sigma.Phi.sq = sigma.Phi.sqCurr, 
                  phi.DrawScale.pred = phi.DrawScale.pred,
                  phi.DrawScale.pred.prev = phi.DrawScale.pred.prev,
                  reu13.df = reu13.df.pred)
    phi.Curr.pred <- phi.pred$phi.New

    # Step A: Update scaling factor.
    b.DrawScale.prev <- b.DrawScale
    phi.DrawScale.prev <- phi.DrawScale
    phi.DrawScale.pred.prev <- phi.DrawScale.pred
    if(iter %/% .CF.AC$renew.iter + 1 != .cubfitsEnv$curr.renew){
      # For each E[b].
      b.DrawScale <- .cubfitsEnv$my.update.DrawScale(
                         "b", b.DrawScale,
                         update.curr.renew = FALSE)
      # For each E[Phi | Phi^{obs}].
      phi.DrawScale <- .cubfitsEnv$my.update.DrawScale(
                         "phi", phi.DrawScale,
                         update.curr.renew = FALSE)
      # For each E[Phi].
      phi.DrawScale.pred <- .cubfitsEnv$my.update.DrawScale(
                              "phi.pred", phi.DrawScale.pred,
                              update.curr.renew = TRUE)
    }

    # Dump parameters out.
    if((iter %% iterThin) == 0){
      thinnedIter <- iter / iterThin + 1
      b.Mat[[thinnedIter]] <- do.call("c", bCurr)
      p.Mat[[thinnedIter]] <- pCurr
      phi.Mat[[thinnedIter]] <- phi.Curr
      phi.Mat.pred[[thinnedIter]] <- phi.Curr.pred
    }
    my.verbose(verbose, iter, report)
    .cubfitsEnv$my.dump(iter, list = c("b.Mat", "p.Mat", "phi.Mat",
                                       "phi.Mat.pred"))
  } # MCMC end

### Return ###
  ret <- list(b.Mat = b.Mat, p.Mat = p.Mat, phi.Mat = phi.Mat,
              phi.Mat.pred = phi.Mat.pred)
  ret
} # End of my.cubpred().

