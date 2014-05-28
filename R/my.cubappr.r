### Special note for cubappr().
### n.G and phi.DrawScale are actually the same role of
### n.G.pred and phi.pred.DrawScale in cubpred().
###

### Function to run MCMC inference for following model:
### Given phi, n, b; across n.G genes:
###   phi ~ rlnorm(n.G, mu.Phi, sigma.Phi)
###   y   ~ for(aa) rmultinom(n.G, invmlogit( phi * b[[aa]] ), n[[aa]] )
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
    phi.pred.DrawScale = .CF.CONF$phi.pred.DrawScale,
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
  my.check.data(phi.Init = phi.Init)

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
  if(.CF.CONF$estimate.bias.Phi){
    nPrior <- nPrior + 1                    # one more for bias.Phi
  }

  ### Storages for saving posterior samples.
  b.Mat <- my.generate.list(NA, nBparams, nSave)     # log(mu) and Delta.t
  p.Mat <- my.generate.list(NA, nPrior, nSave)       # prior parameters
  phi.pred.Mat <- my.generate.list(NA, n.G, nSave)   # E[Phi]

### Initial Parameters ###
  ### Initial values for p first since scaling may change phi.Obs.
  p.Init <- my.pInit(p.Init, phi.Init, model.Phi[1],
                     p.nclass = p.nclass, cub.method = "appr")

  ### Initial values for b.
  if(is.null(bInit)){
    bInitList <- .cubfitsEnv$my.fitMultinomAll(reu13.df.obs, phi.Init, y, n)
    bRInitList <- lapply(bInitList, function(B){ B$R })
    bInit <- lapply(bInitList,
               function(B){
                 B$coefficients +
                 init.b.Scale * backsolve(B$R, rnorm(nrow(B$R)))
               })
  } else{
    bRInitList <- lapply(bInit, function(B){ B$R })
    bInit <- lapply(bInit, function(B){ B$coefficients })
  }
  bInitVec <- unlist(bInit)

### Set current step ###
  ### Set current step for b.
  b.Mat[[1]] <- bInitVec
  b.Curr <- bInit

  ### Set current step for p.
  p.Mat[[1]] <- p.Init 
  p.Curr <- p.Init

  ### Set current step for phi.
  phi.pred.Mat[[1]] <- phi.Init
  phi.Curr <- phi.Init

  ### For hyper-prior parameters.
  hp.param <- list(log.phi.Obs.mean = mean(log(phi.Init)),
                   # hp.sigma.Phi = 1 / sqrt(var(log(phi.Init))),
                   hp.Init = p.Init)

### MCMC here ###
  ### Get length for acceptance and adaptive storage.
  n.p <- 1
  if(.CF.CONF$estimate.bias.Phi){
    n.p <- n.p + 1
  }

  ### Set acceptance rate storage.
  my.set.acceptance(nSave, n.aa, n.p = n.p, n.G.pred = n.G)

  ### Set adaptive storage.
  if(.CF.CONF$estimate.bias.Phi){
    ### Bias of phi is coupled with p parameters.
    p.DrawScale <- c(p.DrawScale, .CF.CONF$bias.Phi.DrawScale)
  }
  my.set.adaptive(nSave,
                  n.aa = n.aa, b.DrawScale = b.DrawScale,
                  n.p = n.p, p.DrawScale = p.DrawScale,
                  n.G.pred = n.G, phi.pred.DrawScale = phi.pred.DrawScale,
                  adaptive = adaptive[1])

  ### Run MCMC iterations.
  my.verbose(verbose, 0, report)
  .cubfitsEnv$my.dump(0, list = c("b.Mat", "p.Mat", "phi.pred.Mat"))

  ### MCMC start.
  for(iter in 1:(nIter + burnin)){
    ### Step 1: Update b using M-H step.
    bUpdate <- .cubfitsEnv$my.drawBConditionalAll(
                 b.Curr, phi.Curr, y, n, reu13.df.obs,
                 bRInitList = bRInitList)
    b.Curr <- lapply(bUpdate, function(U){ U$bNew })

    ### Step 2: Draw other parameters.
    p.Curr <- .cubfitsEnv$my.pPropTypeNoObs(
                n.G, phi.Curr, p.Curr, hp.param)

    ### Step 3: Predict phi using M-H step.
    ###         This is different to cubfits() and cubpred().
    phi.Curr <- my.drawPhiConditionalAllPred(
                  phi.Curr, y, n, b.Curr, p.Curr,
                  reu13.df = reu13.df.obs)

    ### Step A: Update scaling factor.
    if(iter %/% .CF.AC$renew.iter + 1 == .cubfitsEnv$curr.renew){
      my.copy.adaptive()
    } else{
      .cubfitsEnv$my.update.DrawScale("b", update.curr.renew = FALSE)
      .cubfitsEnv$my.update.DrawScale("p", update.curr.renew = FALSE)
      .cubfitsEnv$my.update.DrawScale("phi.pred", update.curr.renew = TRUE)
    }

    ### Dump parameters out.
    if((iter %% iterThin) == 0){
      thinnedIter <- iter / iterThin + 1
      b.Mat[[thinnedIter]] <- do.call("c", b.Curr)
      p.Mat[[thinnedIter]] <- p.Curr
      phi.pred.Mat[[thinnedIter]] <- phi.Curr
    }
    my.verbose(verbose, iter, report)
    .cubfitsEnv$my.dump(iter, list = c("b.Mat", "p.Mat", "phi.pred.Mat"))
  } ### MCMC end.

### Return ###
  ret <- list(b.Mat = b.Mat, p.Mat = p.Mat, phi.pred.Mat = phi.pred.Mat)
  ret
} # End of my.cubappr().

