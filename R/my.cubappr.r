### Special note for cubappr().
# G and phi.DrawScale are actually the same role of
# G.pred and phi.DrawScale.pred in cubpred().
###

# Function to run MCMC inference for following model:
#   Given phi, n, b; across G genes:
#     phi      ~ rlnorm(G, mu.Phi, sigma.Phi)
#     y      ~ for(aa) rmultinom(G, invmlogit( phi * b[[aa]] ), n[[aa]] )
#
# Expects phi.Init as vector of length G.
# 
# Subsequent Gibbs sampler:
# (1) Sample of b | phi, y, using VGAM fit for Gaussian proposal
# (2) Sample mu.Phi, sigma.Phi | x 
# (3) Sample of phi | b, nu.Phi, bsig.Phi, m.Phi, ww, y,
#     using logNormal proposal 

### No observation (phi) is required.
my.cubappr <- function(reu13.df.obs, phi.Init, y, n,
    nIter = 1000, burnin = 100,
    pInit = NULL, bInit = NULL, initBScale = 1,
    phi.DrawScale = 1,
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
  # check phi.Init is well-behaved.
  if(!(all(is.finite(phi.Init)) && all(phi.Init > 0))){
    .cubfitsEnv$my.stop("phi.Init is invalid.")
  }
  if(abs(mean(phi.Init) - 1) > 1e-8 && scale.Phi == "mean_one"){
    .cubfitsEnv$my.stop(paste("mean(phi.Init) =", mean(phi.Init)))
  }
  # Check if sort by ORF and length.
  my.check.rearrange(reu13.df.obs, y, n, phi.Obs = phi.Init)

### Initial Storages ###
  # Setup data structures for results.
  G <- nrow(y[[1]])                         # # of genes
  A <- length(y)                            # # of amino acids
  nsyns <- sapply(y, function(ybit){ dim(ybit)[2] })
                                            # # of synomous codons
  nBparams <- my.ncoef * sum(nsyns - 1)     # total # of regression parameters
  nSave <- (nIter + burnin) / iterThin + 1  # # of space for iterations

  # Storages for saving posterior samples.
  b.Mat <- my.generate.list(NA, nBparams, nSave)     # S/M
  p.Mat <- my.generate.list(NA, 2, nSave)            # hyperparameters
  phi.Mat.pred <- my.generate.list(NA, G, nSave)       # E[Phi]

### Initial Parameters ###
  # Initial values for b.
  if(is.null(bInit)){
    bInitList <- .cubfitsEnv$my.fitMultinomAll(reu13.df.obs, phi.Init, y, n)
    bInit <- lapply(bInitList,
               function(B){
                 B$coefficients + initBScale * backsolve(B$R, rnorm(nrow(B$R)))
               })
  } else{
    bInit <- lapply(bInit,
               function(B){
                 B$coefficients
               })
  }
  bInitVec <- unlist(bInit)

  # Initial values for p.
  if(is.null(pInit)){
    nu.Phi.Init <- mean(log(phi.Init))
    sigmasqObs <- var(log(phi.Init))
    bsig.Phi.Init <- sqrt(.8 * sigmasqObs)     # \sigma_phi
    pInit <- c(nu.Phi = nu.Phi.Init, bsig.Phi = bsig.Phi.Init)
  } else{
    nu.Phi.Init <- pInit["nu.Phi"]
    bsig.Phi.Init <- pInit["bsig.Phi"]
  }

  # Set current step for b.
  b.Mat[[1]] <- bInitVec
  bCurr <- bInit

  # Set current step for p.
  p.Mat[[1]] <- pInit 
  pCurr <- pInit
  nu.Phi.Curr <- nu.Phi.Init
  bsig.Phi.Curr <- bsig.Phi.Init

  # Set current step for phi.
  phi.Mat.pred[[1]] <- phi.Init
  phi.Curr <- phi.Init

  # For hyper-prior parameters.
  log.Phi.Obs.mean <- mean(log(phi.Curr))

### MCMC here ###
  # Set acceptance rate storage.
  my.set.acceptance(nSave, A, G.pred = G)

  # Set adaptive storage.
  my.set.adaptive(nSave, G.pred = G, phi.DrawScale.pred = phi.DrawScale,
                  adaptive = adaptive[1])
  phi.DrawScale.pred <- .cubfitsEnv$DrawScale$phi.pred[[1]]
  phi.DrawScale.pred.prev <- phi.DrawScale.pred

  # Run MCMC iterations.
  my.verbose(verbose, 0, report)
  .cubfitsEnv$my.dump(0, list = c("b.Mat", "p.Mat", "phi.Mat.pred"))

  # MCMC start.
  for(iter in 1:(nIter + burnin)){
    # Step 1: Update b using M-H step.
    bUpdates <- .cubfitsEnv$my.drawBConditionalAll(bCurr, phi.Curr, y, n,
                                                   reu13.df.obs)
    bCurr  <- lapply(bUpdates, function(U){ U$bNew })

    # Step 2: Draw other parameters.
    log.Phi.Curr <- log(phi.Curr)
    pUpdate <- .cubfitsEnv$my.pPropTypeNoObs(G, log.Phi.Obs.mean, log.Phi.Curr,
                 nu.Phi.Curr, bsig.Phi.Curr)
    nu.Phi.Curr <- pUpdate$ppCurr$y0
    bsig.Phi.Curr <- pUpdate$ppCurr$bb
    pCurr <- c(nu.Phi.Curr, bsig.Phi.Curr)
    mu.Phi.Curr <- nu.Phi.Curr
    sigma.Phi.sqCurr <- bsig.Phi.Curr^2

    # Step 3: Predict phi using M-H step.
    #         This is different to cubfits() and cubpred().
    phi.pred <- my.drawPhiConditionalAllPred(phi.Curr, y, n,
               bCurr,
               mu.Phi = mu.Phi.Curr,
               sigma.Phi.sq = sigma.Phi.sqCurr, 
               phi.DrawScale.pred = phi.DrawScale.pred,
               phi.DrawScale.pred.prev = phi.DrawScale.pred.prev,
               reu13.df = reu13.df.obs)
    phi.Curr <- phi.pred$phi.New

    # Step A: Update scaling factor for each E[Phi].
    phi.DrawScale.pred.prev <- phi.DrawScale.pred
    if(iter %/% .CF.AC$renew.iter + 1 != .cubfitsEnv$curr.renew){
      phi.DrawScale.pred <- .cubfitsEnv$my.update.DrawScale("phi.pred", phi.DrawScale)
    }

    # Dump parameters out.
    if((iter %% iterThin) == 0){
      thinnedIter <- iter / iterThin + 1
      b.Mat[[thinnedIter]] <- do.call("c", bCurr)
      p.Mat[[thinnedIter]] <- pCurr
      phi.Mat.pred[[thinnedIter]] <- phi.Curr
    }
    my.verbose(verbose, iter, report)
    .cubfitsEnv$my.dump(iter, list = c("b.Mat", "p.Mat", "phi.Mat.pred"))
  } # MCMC end.

### Return ###
  ret <- list(b.Mat = b.Mat, p.Mat = p.Mat, phi.Mat.pred = phi.Mat.pred)
  ret
} # End of my.cubpred().

