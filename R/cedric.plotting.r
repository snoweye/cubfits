plotPTraces <- function(ret.p.Mat)
{
  ret.p.Mat <- do.call("rbind", ret.p.Mat)
  n.traces <- dim(ret.p.Mat)[2]
  if(n.traces==2)
  {
    ylabs <- c("log(M)", expression("log(s"[phi]*")"))
  }else{
    ylabs <- c(expression("log(s"[epsilon]*")"), "log(M)", expression("log(s"[phi]*")"), "log(K)")
  }
  
  par(oma=c(1,1,2,1), mgp=c(2,1,0), mar = c(3,4,2,1), mfrow=c(2, ceiling(n.traces/2) ))
  for(i in 1:n.traces)
  {
    plot(ret.p.Mat[, i], xlab="iteration", ylab=ylabs[i])
  }
}

plotExpectedPhiTrace <- function(phi.Mat, 
                                 main=expression(paste("Trace E[", phi, "]", sep="")), 
                                 xlab="iteration", ylab=expression(paste("E[", phi, "]", sep="")))
{
  
  phi.Mat <- do.call("cbind", phi.Mat)
  phi.Mat <- colMeans(phi.Mat)
  
  plot(phi.Mat, xlab=xlab, ylab=ylab, main=main)  
}

plotCUB <- function(reu13.df.obs, ret.b.Mat, phi.bin, estim.phi, n.use.iter=2000, rescale=F,
                     main="CUB", model.label=c("True Model"), model.lty=1)
{
  ### Arrange data.
  aa.names <- names(reu13.df.obs)
  #phi.bin <- phi.bin * phi.scale
  phi.bin.lim <- range(c(phi.bin, estim.phi))
  
  lbound <- length(ret.b.Mat)-n.use.iter
  ubound <- length(ret.b.Mat)
  b.mat <- do.call(cbind, ret.b.Mat[lbound:ubound])
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




plotTraces <- function(ret.b.Mat, names.aa, param = c("logmu", "deltat"), main="AA parameter trace")
{
  b.mat <- convert.bVec.to.b(ret.b.Mat[[1]], names.aa)
  b.mat <- convert.b.to.bVec(b.mat)
  names.b <- names(b.mat)
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
  
  x <- 1:length(ret.b.Mat)
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
    trace <- lapply(1:length(ret.b.Mat), function(i){ ret.b.Mat[[i]][id.tmp] })
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
