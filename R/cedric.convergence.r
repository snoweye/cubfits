
appendCUBresults <- function(res, to)
{
  res.list.length <- length(res$b.Mat)
  init.list.length <- length(to$b.Mat)
  if(init.list.length == 0)
  {
    to <- list()
    to$b.Init <- res$b.Init
    to$b.RInit <- res$b.RInit
    to$p.Init <- res$p.Init
    to$phi.Init <- res$phi.Init
  }
  
  for(i in 1:res.list.length)
  {
    to$b.Mat[init.list.length + i] <- res$b.Mat[i]
    to$p.Mat[init.list.length + i] <- res$p.Mat[i]
    
    if("phi.Mat" %in% names(res)){
      to$phi.Mat[init.list.length + i] <- res$phi.Mat[i]
    }
    if("phi.pred.Mat" %in% names(res)){
      to$phi.pred.Mat[init.list.length + i] <- res$phi.pred.Mat[i]
    }
  }
  return(to)
}

isConverged <- function(chains, niter, epsilon=0.05, thin=10, teston=c("phi", "sphi"))
{
  result <- FALSE
  dataobj <- 0
  start <- 0
  
  list <- mcmc.list()
  for(i in 1:length(chains))
  {  
    if(teston[1] == "sphi")
    {
      p.mat <- do.call("rbind", chains[[i]]$p.Mat)
      index <- 3 # index of sphi for with x obs
      if(dim(p.mat)[2] == 2){index <- 2} # index of sphi for without x obs
      dataobj <- p.mat[, index]
    }else if(teston[1] == "phi"){ # test on phis
      if("phi.Mat" %in% names(chains[[i]])) # with x obs 
      {
        phi.mat <- do.call("rbind", chains[[i]]$phi.Mat)
      }
      if("phi.pred.Mat" %in% names(chains[[i]])) # without x obs
      {
        phi.mat <- do.call("rbind", chains[[i]]$phi.pred.Mat)
      }
      dataobj <- phi.mat
    }else{
      stop("convergence test can not be perfomed on choosen data\n")
    }
#     if(length(dataobj) > 20000)
#     {
#       cat(paste(summary(dataobj), "\n"))
#     }
#     cat(paste(length(dataobj), " ", any(is.na(data)), " "))
    start <- length(dataobj) - niter
    mcmcobj <- mcmc(data=dataobj, start=start, thin=thin)
    list[[i]] <- mcmcobj
  }
  
  #perform Gelman and Rubin's convergence test
#   cat(paste(length(list), " "))
  
  diag <- gelman.diag(list, autoburnin=FALSE)
  
  if(teston[1] == "sphi") # scalar test on s phi
  {
    cat(paste("Gelman score: ", diag[[1]][1], "\n"))
    result <- abs(diag[[1]][1] - 1) < epsilon 
    ret <- list(isConverged=result, gelman=diag[[1]][1])
  }else # else is enough here. The correctness of the method was determined above and leaves only two options
  { # multivariate test on all phi values
    cat(paste("Gelman score: ", diag$mpsrf, "\n"))
    result <- abs(diag$mpsrf - 1) < epsilon
    ret <- list(isConverged=result, gelman=diag$mpsrf)
  }
  
  return(ret)
}

cubmultichain <- function(cubmethod, niter, reset.qr, seeds, teston=c("phi", "sphi"),
                          min=0, max=160000, nchains=2, thin=10, eps=0.05, ncores=2, ...)
{
  cl <- makeCluster(ncores)
  registerDoSNOW(cl)
  
  ########################
  ## checking arguments ##
  ########################
  if(min > max){
    warning("min iterations > max iterations. setting max = min\n")
    max <- min
  }
  if(nchains < 2){
    warning("number of chains not sufficient. setting chains to 2\n")
    nchains <- 2
  }
  if(!cubmethod %in% c("cubfits", "cubappr", "cubpred")){
    stop(paste("Unkown method: ", cubmethod, "!\n"))
  }
  
  ## arguments for cub methods
  input_list <- as.list(list(...))
  
  if("p.Init" %in% names(input_list)){
    p.init <- input_list$p.Init
    if(length(p.init) < nchains){
      stop("p.Init does not contain initial p vectors for every chain!\n")
    }
    input_list$p.Init <- NULL
  }else{
    p.init <- list(NULL)
    length(p.init) <- nchains  
  }
  if("b.RInit" %in% names(input_list)){
    b.rinit <- input_list$b.RInit
    input_list$b.RInit <- NULL
  }else{
    b.rinit <- list(NULL)
    length(b.rinit) <- nchains  
  }
  
  init.phi <- list()
  init.pred.phi <- list()
  if("phi.Init" %in% names(input_list)){
    init.phi <- input_list$phi.Init
    input_list$phi.Init <- NULL
  }
  if("phi.pred.Init" %in% names(input_list)){
    init.pred.phi <- input_list$phi.pred.Init
    input_list$phi.pred.Init <- NULL
  }
  if(".CF.CT" %in% names(input_list)){
    .CF.CT <- input_list$.CF.CT
    input_list$.CF.CT <- NULL
  }else{
    .CF.CT <- eval(parse(text = "cubfits::.CF.CT"))
  } 
  if(".CF.CONF" %in% names(input_list)){
    .CF.CONF <- input_list$.CF.CONF
    input_list$.CF.CONF <- NULL
  }else{
    .CF.CONF <- eval(parse(text = "cubfits::.CF.CONF"))
  }   
  results <- list()
  length(results) <- nchains
  if(is.null(seeds)){
    seeds <- round(runif(nchains, 1, 100000))
  }
  
  #############################################################
  ## running chains in parralel and checking for convergence ##
  #############################################################
  j <- 1
  gel.res <- 0
  converged <- FALSE
  while(!converged)
  { 
    ## run chains in parallel
    #for(i in nchains:1) # for debuging
    res <- foreach(i = 1:nchains) %dopar%
    {
      suppressMessages(library(cubfits, quietly = TRUE))
      
      .GlobalEnv$.CF.CT <- .CF.CT
      .GlobalEnv$.CF.CONF <- .CF.CONF
      set.seed(seeds[i])
      if(cubmethod == "cubfits"){
        do.call(cubfits, c(input_list, list(phi.Init = init.phi[[i]]), list(p.Init = p.init[[i]]), list(b.RInit = b.rinit[[i]])))
      }else if(cubmethod == "cubappr"){
        do.call(cubappr, c(input_list, list(phi.pred.Init = init.pred.phi[[i]]), list(p.Init = p.init[[i]]), list(b.RInit = b.rinit[[i]])))
      }else if(cubmethod == "cubpred"){
        do.call(cubpred, c(input_list, list(phi.Init = init.phi[[i]]), list(phi.pred.Init = init.pred.phi[[i]]), list(p.Init = p.init[[i]]), list(b.RInit = b.rinit[[i]])))
      }
    }
    ## append chains and get new initial values for restart
    for(i in 1:nchains)
    {
      
      if(cubmethod == "cubfits" | cubmethod == "cubpred")
      {
        init.phi[[i]] <- normalize.data.set(res[[i]]$phi.Mat[[length(res[[i]]$phi.Mat)]])
      }
      if(cubmethod == "cubappr" | cubmethod == "cubpred")
      {
        init.pred.phi[[i]] <- normalize.data.set(res[[i]]$phi.pred.Mat[[length(res[[i]]$phi.pred.Mat)]])
      }
      p.init[[i]] <- res[[i]]$p.Mat[[length(res[[i]]$p.Mat)]]
      results[[i]] <- appendCUBresults(res[[i]], results[[i]])
      
      if(length(results[[i]]$p.Mat) < reset.qr) # reset the "cov" only in the begining
      {
        b.rinit[i] <- list(NULL)
      }else{ # use the same matrix every time after some "burnin"
        b.rinit[[i]] <- res[[i]]$b.RInit
      }
    }
    curiter <- length(results[[1]]$p.Mat)
    
    ## Do convergence test
    if(curiter > niter){ #if there are not enough iterations, just keep goint until we have enough for a convergence test
      gelman <- isConverged(results, niter, epsilon=eps, thin=thin, teston=teston)
      gel.res[j] <- gelman$gelman
      converged <- gelman$isConverged
      j <- j + 1
    }
    
    #check if we have at least min iterations
    if(curiter < min){converged <- FALSE}
    #check if max iteration limit is reached
    if(curiter > max){converged <- TRUE}
  }
## return full length chains
return(list(chains=results, convergence=gel.res))
}
