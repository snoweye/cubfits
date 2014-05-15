### Special note for cubappr().
###   n.G and phi.DrawScale are actually the same role of
###   n.G.pred and phi.DrawScale.pred in cubpred().
###

### Function to run MCMC inference for following model:
###   Given phi, n, b; across n.G genes:
###     phi      ~ rlnorm(n.G, mu.Phi, sigma.Phi)
###     y      ~ for(aa) rmultinom(n.G, invmlogit( phi * b[[aa]] ), n[[aa]] )
###
### Expects phi.Init as vector of length n.G.
### 
### Subsequent Gibbs sampler:
### (1) Sample of b | phi, y, using VGAM fit for Gaussian proposal
### (2) Sample mu.Phi, sigma.Phi | x 
### (3) Sample of phi | b, nu.Phi, sigma.Phi, m.Phi, ww, y,
###     using logNormal proposal 

### No observation (phi) is required.
my.cubappr <- function(reu13.df.obs, phi.Init, y, n,
    nIter = 1000, burnin = 100,
    bInit = NULL, init.b.Scale = .CF.CONF$init.b.Scale,
        b.DrawScale = .CF.CONF$b.DrawScale,
    p.Init = NULL, p.nclass = .CF.CONF$p.nclass,
        p.DrawScale = .CF.CONF$p.DrawScale,
    phi.DrawScale.pred = .CF.CONF$phi.DrawScale.pred,
    model = .CF.CT$model[1],
    model.Phi = .CF.CT$model.Phi[1],
    adaptive = .CF.CT$adaptive[1],
    verbose = .CF.DP$verbose,
    iterThin = .CF.DP$iterThin,
    report = .CF.DP$report){

### Setup functions ###
  ### Setup function pointers by type or model.
  my.function <- my.init.function(model = model[1], model.Phi = model.Phi[1],
                                  adaptive = adaptive[1])
  my.ncoef <- my.function$my.ncoef

### Check Data ###
  ### check phi.Init is well-behaved.
  if(!(all(is.finite(phi.Init)) && all(phi.Init > 0))){
    .cubfitsEnv$my.stop("phi.Init is invalid.")
  }
  if(abs(mean(phi.Init) - 1) > 1e-8 && .CF.CONF$scale.phi){
    .cubfitsEnv$my.stop(paste("mean(phi.Init) =", mean(phi.Init)))
  }
  ### Check if sort by ORF and length.
  my.check.rearrange(reu13.df.obs, y, n, phi.Obs = phi.Init)

### Initial Storages ###
  ### Setup data structures for results.
  n.G <- nrow(y[[1]])                       # # of genes
  n.aa <- length(y)                         # # of amino acids
  nsyns <- sapply(y, function(ybit){ dim(ybit)[2] })
                                            # # of synomous codons
  nBparams <- my.ncoef * sum(nsyns - 1)     # total # of regression parameters
  nSave <- (nIter + burnin) / iterThin + 1  # # of space for iterations
  nPrior <- 2                               # # of prior parameters
  if(model.Phi == "logmixture"){
    nPrior <- 3 * p.nclass 
  }

  ### Storages for saving posterior samples.
  b.Mat <- my.generate.list(NA, nBparams, nSave)     # S/M
  p.Mat <- my.generate.list(NA, nPrior, nSave)       # prior parameters
  phi.Mat.pred <- my.generate.list(NA, n.G, nSave)   # E[Phi]

### Initial Parameters ###
  ### Initial values for b.
  bInitList <- .cubfitsEnv$my.fitMultinomAll(reu13.df.obs, phi.Init, y, n)
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

  ### Initial values for p.
  p.Init <- my.pInit(p.Init, phi.Init, model.Phi[1],
                     p.nclass = p.nclass, cub.method = "appr")

  ### Set current step for b.
  b.Mat[[1]] <- bInitVec
  b.Curr <- bInit

  ### Set current step for p.
  p.Mat[[1]] <- p.Init 
  p.Curr <- p.Init

  ### Set current step for phi.
  phi.Mat.pred[[1]] <- phi.Init
  phi.Curr <- phi.Init

  ### For hyper-prior parameters.
  hp.param <- list(log.phi.Obs.mean = mean(log(phi.Init)),
                   # hp.sigma.Phi = 1 / sqrt(var(log(phi.Init))),
                   hp.Init = p.Init)

### MCMC here ###
  ### Set acceptance rate storage.
  my.set.acceptance(nSave, n.aa, n.p = 1,
                    n.G.pred = n.G)

  ### Set adaptive storage.
  my.set.adaptive(nSave,
                  n.aa = n.aa, b.DrawScale = b.DrawScale,
                  n.p = 1, p.DrawScale = p.DrawScale,
                  n.G.pred = n.G, phi.DrawScale.pred = phi.DrawScale.pred,
                  adaptive = adaptive[1])
  b.DrawScale <- .cubfitsEnv$DrawScale$b[[1]]
  b.DrawScale.prev <- b.DrawScale
  p.DrawScale <- .cubfitsEnv$DrawScale$p[[1]]
  p.DrawScale.prev <- p.DrawScale
  phi.DrawScale.pred <- .cubfitsEnv$DrawScale$phi.pred[[1]]
  phi.DrawScale.pred.prev <- phi.DrawScale.pred

  ### Run MCMC iterations.
  my.verbose(verbose, 0, report)
  .cubfitsEnv$my.dump(0, list = c("b.Mat", "p.Mat", "phi.Mat.pred"))

  ### MCMC start.
  for(iter in 1:(nIter + burnin)){
    ### Step 1: Update b using M-H step.
    bUpdate <- .cubfitsEnv$my.drawBConditionalAll(
                 b.Curr, phi.Curr, y, n, reu13.df.obs,
                 bRInitList = bRInitList,
                 b.DrawScale = b.DrawScale,
                 b.DrawScale.prev = b.DrawScale.prev)
    b.Curr <- lapply(bUpdate, function(U){ U$bNew })

    ### Step 2: Draw other parameters.
    p.Curr <- .cubfitsEnv$my.pPropTypeNoObs(
                n.G, phi.Curr, p.Curr, hp.param,
                p.DrawScale = p.DrawScale,
                p.DrawScale.prev = p.DrawScale.prev)

    ### Step 3: Predict phi using M-H step.
    ###         This is different to cubfits() and cubpred().
    phi.Curr <- my.drawPhiConditionalAllPred(
                  phi.Curr, y, n, b.Curr, p.Curr,
                  phi.DrawScale.pred = phi.DrawScale.pred,
                  phi.DrawScale.pred.prev = phi.DrawScale.pred.prev,
                  reu13.df = reu13.df.obs)

    ### Step A: Update scaling factor.
    b.DrawScale.prev <- b.DrawScale
    p.DrawScale.prev <- p.DrawScale
    phi.DrawScale.pred.prev <- phi.DrawScale.pred
    if(iter %/% .CF.AC$renew.iter + 1 != .cubfitsEnv$curr.renew){
      ### For each E[b].
      b.DrawScale <- .cubfitsEnv$my.update.DrawScale(
                       "b", b.DrawScale,
                       update.curr.renew = FALSE)
      ### For prior.
      p.DrawScale <- .cubfitsEnv$my.update.DrawScale(
                       "p", p.DrawScale,
                       update.curr.renew = FALSE)
      ### For each E[Phi].
      phi.DrawScale.pred <- .cubfitsEnv$my.update.DrawScale(
                              "phi.pred", phi.DrawScale.pred,
                              update.curr.renew = TRUE)
    }

    ### Dump parameters out.
    if((iter %% iterThin) == 0){
      thinnedIter <- iter / iterThin + 1
      b.Mat[[thinnedIter]] <- do.call("c", b.Curr)
      p.Mat[[thinnedIter]] <- p.Curr
      phi.Mat.pred[[thinnedIter]] <- phi.Curr
    }
    my.verbose(verbose, iter, report)
    .cubfitsEnv$my.dump(iter, list = c("b.Mat", "p.Mat", "phi.Mat.pred"))
  } ### MCMC end.

### Return ###
  ret <- list(b.Mat = b.Mat, p.Mat = p.Mat, phi.Mat.pred = phi.Mat.pred)
  ret
} # End of my.cubpred().

