plotPTraces <- function(pMat, ...)
{
  input_list <- as.list(list(...))
  pMat <- do.call("rbind", pMat)
  
  if("main" %in% names(input_list)){
    main <- input_list$main
    input_list$main <- NULL
  }else{
    main <- "Hyperparameter Traces"
  }
  if("xlab" %in% names(input_list)){
    xlab <- input_list$xlab
    input_list$xlab <- NULL
  }else{
    xlab <- "Iteration"
  }
  if("ylab" %in% names(input_list)){
    ylab <- input_list$ylab
    input_list$ylab <- NULL
  }else{
    n.traces <- dim(pMat)[2]
    if(n.traces==2)
    {
      ylabs <- c("M", expression("s"[phi]))
    }else{
      ylabs <- c(expression("s"[epsilon]), "M", expression("s"[phi]), "K")
    }
  }
  if("type" %in% names(input_list)){
    type <- input_list$type
    input_list$type <- NULL
  }else{
    type <- "l"
  }
  
  par(oma=c(1,1,2,1), mgp=c(2,1,0), mar = c(3,4,2,1), mfrow=c(2, ceiling(n.traces/2) ))
  for(i in 1:n.traces)
  {
    do.call(plot, c(input_list, list(x=pMat[, i]), list(xlab=xlab), list(ylab=ylabs[i]), list(type=type)) )
    #plot(pMat[, i], xlab=xlab, ylab=ylabs[i], type=type)
  }
  title(main=main, outer=T)
}

plotExpectedPhiTrace <- function(phiMat, ...)
{
  input_list <- as.list(list(...))
  
  if("main" %in% names(input_list)){
    main <- input_list$main
    input_list$main <- NULL
  }else{
    main <- expression(paste("Trace E[", phi, "]", sep=""))
  }
  if("xlab" %in% names(input_list)){
    xlab <- input_list$xlab
    input_list$xlab <- NULL
  }else{
    xlab <- "Iteration"
  }
  if("ylab" %in% names(input_list)){
    ylab <- input_list$ylab
    input_list$ylab <- NULL
  }else{
    ylab <- expression(paste("E[", phi, "]", sep=""))
  }
  if("type" %in% names(input_list)){
    type <- input_list$type
    input_list$type <- NULL
  }else{
    type <- "l"
  }
  phiMat <- do.call("cbind", phiMat)
  phiMat <- colMeans(phiMat)
  
  
  do.call( plot, c(input_list, list(x=phiMat), list(xlab=xlab), list(ylab=ylab), list(main=main), list(type=type)) )
  #plot(phi.Mat, xlab=xlab, ylab=ylab, main=main, type=type)  
}

plotCUB <- function(reu13.df.obs, bMat, phi.bin, phiMat, n.use.samples=2000, rescale=F,
                     main="CUB", model.label=c("True Model"), model.lty=1)
{
  ### Arrange data.
  aa.names <- names(reu13.df.obs)
  #phi.bin <- phi.bin * phi.scale
  phi.bin.lim <- range(c(phi.bin, phiMat))
  
  lbound <- length(bMat)-n.use.samples
  ubound <- length(bMat)
  b.mat <- do.call(cbind, bMat[lbound:ubound])
  Eb <- rowMeans(b.mat)
  Eb <- convert.bVec.to.b(Eb, aa.names)
  
  ### Compute.
  ret.phi.bin <- prop.bin.roc(reu13.df.obs, phi.bin)
  predict.roc <- prop.model.roc(Eb, phi.bin.lim)
  
  ### Fix xlim at log10 scale. 
  lim.bin <- range(log10(ret.phi.bin[[1]]$center))
  lim.model <- range(log10(predict.roc[[1]]$center))
  xlim <- c(lim.bin[1] - (lim.bin[2] - lim.bin[1]) / 4,
            max(lim.bin[2], lim.model[2]))
  
  
  mat <- matrix(c(rep(1, 5), 2:21, rep(22, 5)),
                nrow = 6, ncol = 5, byrow = TRUE)
  mat <- cbind(rep(23, 6), mat, rep(24, 6))
  nf <- layout(mat, c(3, rep(8, 5), 2), c(3, 8, 8, 8, 8, 3), respect = FALSE)
  ### Plot title.
  par(mar = c(0, 0, 0, 0))
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.6, main)
  text(0.5, 0.4, date(), cex = 0.6)
  
  ### Plot results.
  for(i.aa in 1:length(aa.names))
  {
    tmp.obs <- ret.phi.bin[[i.aa]]
    tmp.roc <- predict.roc[[i.aa]]
    
    plotbin(tmp.obs, tmp.roc, main = "", xlab = "", ylab = "",
            lty = model.lty, axes = FALSE, xlim = xlim)
    box()
    text(0, 1, aa.names[i.aa], cex = 1.5)
    if(i.aa %in% c(1, 6, 11, 16)){
      axis(2)
    }
    if(i.aa %in% 15:19){
      axis(1)
    }
    if(i.aa %in% 1:5){
      axis(3)
    }
    if(i.aa %in% c(5, 10, 15)){
      axis(4)
    }
    axis(1, tck = 0.02, labels = FALSE)
    axis(2, tck = 0.02, labels = FALSE)
    axis(3, tck = 0.02, labels = FALSE)
    axis(4, tck = 0.02, labels = FALSE)
  }
  
  ## adding a histogram of phi values to plot
  #hist(log10(phiMat), xlab=expression(phi), main="")
    
  ### Add label.
  
  plot(NULL, NULL, axes = FALSE, main = "", xlab = "", ylab = "",
       xlim = c(0, 1), ylim = c(0, 1))
  legend(0.1, 0.8, model.label, lty = model.lty, box.lty = 0)
  
  ### Plot xlab.
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.5, "Production Rate (log10)")
  
  ### Plot ylab.
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.5, "Propotion", srt = 90)
}


plotBPosterior <- function(bMat, names.aa, interval, param = c("logmu", "deltat"), main="AA parameter posterior", nclass=100)
{
  bmat <- convert.bVec.to.b(bMat[[1]], names.aa)
  bmat <- convert.b.to.bVec(bmat)
  names.b <- names(bmat)
  id.intercept <- grep("Intercept", names.b)
  id.slope <- 1:length(names.b)
  id.slope <- id.slope[-id.intercept]
  
  
  id.plot <- rep(0, length(names.b))
  if(param[1] == "logmu"){
    xlab <- expression(paste("log ( ", mu, " )"))
    id.plot[id.intercept] <- id.intercept
  } else if(param[1] == "deltat"){
    xlab <- expression(paste(Delta, "t"))
    id.plot[id.slope] <- id.slope
  }  
  
  nf <- layout(matrix(c(rep(1, 5), 2:21), nrow = 5, ncol = 5, byrow = TRUE),
               rep(1, 5), c(2, 8, 8, 8, 8), respect = FALSE)
  ### Plot title.
  par(mar = c(0, 0, 0, 0))
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.6, main)
  text(0.5, 0.4, date(), cex = 0.6)
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  
  
  ### Plot by aa.
  for(i.aa in names.aa){
    id.tmp <- grepl(i.aa, names.b) & id.plot
    trace <- lapply(1:length(bMat), function(i){ bMat[[i]][id.tmp] })
    trace <- do.call("rbind", trace)
    if(length(trace) == 0) next
    
    ymax <- vector(mode = "numeric", length = length(id.tmp))
    for(i in 1:sum(id.tmp))
    {
      ymax[i] <- max(hist(trace[interval, i], plot=F, nclass=nclass)$counts)
    }
    ylim <- c(0, max(ymax))
    xlim <- range(trace[interval, ])
    plot(NULL, NULL, xlim = xlim, ylim = ylim,
         xlab = xlab, ylab = "Frequency", main = i.aa)
    plot.order <- order(apply(trace, 2, sd), decreasing = TRUE)
    for(i.codon in plot.order){
      hist(trace[interval, i.codon], add=T, nclass=nclass, col=.CF.PT$color[i.codon], lty=0)
    } 
  }
  trace <- lapply(1:length(bMat), function(i){
    bMat[[i]][id.plot]
  })
  trace <- do.call("rbind", trace)
  ylim = c(0, max(hist(trace[interval, ], plot=F, nclass=5*nclass)$counts) )
  xlim <- range(trace[interval, ])
  plot(NULL, NULL, xlim = xlim, ylim = ylim,
       xlab = xlab, ylab = "Frequency", main = "Combined")
  
  hist(trace[interval, ], nclass=500, col="blue", add=T, lty=0, xlab=xlab, ylab=Frequency)
}

plotTraces <- function(bMat, names.aa, param = c("logmu", "deltat"), main="AA parameter trace")
{  
  
  bmat <- convert.bVec.to.b(bMat[[1]], names.aa)
  bmat <- convert.b.to.bVec(bmat)
  names.b <- names(bmat)
  id.intercept <- grep("Intercept", names.b)
  id.slope <- 1:length(names.b)
  id.slope <- id.slope[-id.intercept]
  
  
  id.plot <- rep(0, length(names.b))
  if(param[1] == "logmu"){
    ylab <- expression(paste("log ( ", mu, " )"))
    id.plot[id.intercept] <- id.intercept
  } else if(param[1] == "deltat"){
    ylab <- expression(paste(Delta, "t"))
    id.plot[id.slope] <- id.slope
  }  
  
  x <- 1:length(bMat)
  xlim <- range(x)
  
  ### Trace plot.
  nf <- layout(matrix(c(rep(1, 5), 2:21), nrow = 5, ncol = 5, byrow = TRUE),
               rep(1, 5), c(2, 8, 8, 8, 8), respect = FALSE)
  ### Plot title.
  par(mar = c(0, 0, 0, 0))
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.6, main)
  text(0.5, 0.4, date(), cex = 0.6)
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  
  ### Plot by aa.
  for(i.aa in names.aa){
    id.tmp <- grepl(i.aa, names.b) & id.plot
    trace <- lapply(1:length(bMat), function(i){ bMat[[i]][id.tmp] })
    trace <- do.call("rbind", trace)
    if(length(trace) == 0) next
    
    ylim <- range(trace, na.rm=T)
    plot(NULL, NULL, xlim = xlim, ylim = ylim,
         xlab = "Iterations", ylab = ylab, main = i.aa)
    plot.order <- order(apply(trace, 2, sd), decreasing = TRUE)
    for(i.codon in plot.order){
      lines(x = x, y = trace[, i.codon], col = .CF.PT$color[i.codon])
    } 
  }
}
