plottrace.param <- function(ret.b.Mat, names.b,
    names.aa, id.intercept, workflow.name, i.case, model, bInit = NULL,
    param = c("logmu", "deltat")){
  if(param[1] == "logmu"){
    ylab <- "log(mu)"
    id.plot <- 1:nrow(names.b) %in% id.intercept
  } else if(param[1] == "deltat"){
    ylab <- "Delta.t"
    id.plot <- !(1:nrow(names.b) %in% id.intercept)
  } else{
    stop("param is not found.")
  }

  x <- 1:length(ret.b.Mat)
  xlim <- range(x)

  # Trace plot.
  nf <- layout(matrix(c(rep(1, 5), 2:21), nrow = 5, ncol = 5, byrow = TRUE),
               rep(1, 5), c(1, 8, 8, 8, 8), respect = FALSE)
  # Plot title.
  par(mar = c(0, 0, 0, 0))
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.5,
       paste(workflow.name, ", ", get.case.main(i.case, model), sep = ""))
  par(mar = c(5.1, 4.1, 4.1, 2.1))

  # Plot by aa.
  for(i.aa in names.aa){
    id.tmp <- (names.b == i.aa) & id.plot
    trace <- lapply(1:length(ret.b.Mat), function(i){ ret.b.Mat[[i]][id.tmp] })
    trace <- do.call("rbind", trace)

    ylim <- range(trace)
    plot(NULL, NULL, xlim = xlim, ylim = ylim,
         xlab = "Iterations", ylab = ylab, main = i.aa)
    plot.order <- order(apply(trace, 2, sd), decreasing = TRUE)
    for(i.codon in plot.order){
      lines(x = x, y = trace[, i.codon], col = .CF.PT$color[i.codon])
    }

    if(!is.null(bInit)){
      for(i.b in bInit[id.tmp]){
        abline(h = i.b, lty = 2)
      }
    }
  }
} # End of plottrace.param().


plottrace.meanEPhi <- function(ret.phi.Mat, workflow.name, i.case, model){
  x <- 1:length(ret.phi.Mat)
  xlim <- range(x)

  # Trace of mean of expected expression.
  trace <- lapply(1:length(ret.phi.Mat), function(i){ mean(ret.phi.Mat[[i]]) })
  trace <- do.call("c", trace)

  ylim <- range(c(range(trace), 1))
  plot(NULL, NULL, xlim = xlim, ylim = ylim,
       xlab = "Iterations", ylab = "Mean of EPhi")
  mtext(paste(workflow.name, ", ", get.case.main(i.case, model), sep = ""),
        line = 3, cex = 0.6)
  lines(x = x, y = trace)
  abline(h = 1, col = 2)
} # End of plottrace.meanEPhi().

