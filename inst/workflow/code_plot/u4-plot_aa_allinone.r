### This is for simulation only, ploting delta.t * phi aganst true values.

plot.aa.allinone <- function(
    i.aa, id.label, tl.codon,                    # For AA.
    label.negsel.true,                           # For label.
    bInit, EPhi,                                 # For true.
    b.mcmc, t.phi.mcmc,                          # For unscaled results.
    b.negsel.PM, phi.PM,                         # For scaled results.
    workflow.name, i.case, model){
  # Get right b for AA.
  b.AA.true <- bInit[id.label]
  b.AA.PM <- matrix(b.mcmc[id.label,], nrow = tl.codon)
  b.AA.PM.scaled <- b.negsel.PM[id.label]

  # Get true deltat * phi
  deltat.phi.true <- NULL
  for(i.codon in 1:tl.codon){
    deltat.phi.true <- cbind(deltat.phi.true,
                             b.AA.true[i.codon] * EPhi)
  }

  # Get mean of unscaled deltat * phi
  deltat.phi.PM <- NULL
  for(i.codon in 1:tl.codon){
    deltat.phi.PM <- cbind(deltat.phi.PM,
                           colMeans(b.AA.PM[i.codon,] * t.phi.mcmc))
  }

  # Get scaled deltat * scaled phi
  deltat.phi.PM.scaled <- NULL
  for(i.codon in 1:tl.codon){
    deltat.phi.PM.scaled <- cbind(deltat.phi.PM.scaled,
                                  b.AA.PM.scaled[i.codon] * phi.PM)
  }

  ### AA plot.
  nf <- layout(rbind(rep(1, 5), matrix(1 + 1:(5 * 2), nrow = 2)),
               rep(1, 5), c(1, 8, 8), respect = FALSE)

  # Plot title.
  par(mar = c(0, 0, 0, 0))
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.5,
       paste(workflow.name, ", ", get.case.main(i.case, model),
             sep = ""))
  par(mar = c(5.1, 4.1, 4.1, 2.1))

  # Plot codon with large range first.
  order.plot <- apply(deltat.phi.true, 2, function(x){ abs(diff(range(x))) })
  u.codon <- gsub(".\\.(.*)", "\\1", id.label)
  u.codon.sorted <- sort(u.codon)
  color <- get.color(u.codon.sorted)

  # Plot top row.
  lim <- range(range(deltat.phi.true), range(deltat.phi.PM))
  plot(NULL, NULL,
       xlim = lim, ylim = lim,
       xlab = "true Delta.t * true Phi", ylab = "mean (Delta.t * Phi)",
       main = i.aa)
  for(i.codon in order.plot){
    points(deltat.phi.true[, i.codon], deltat.phi.PM[, i.codon],
           col = color[u.codon.sorted == u.codon[i.codon]])
  }

  # Plot bottom row.
  lim <- range(range(deltat.phi.true), range(deltat.phi.PM.scaled))
  plot(NULL, NULL,
       xlim = lim, ylim = lim,
       xlab = "true Delta.t * true Phi",
       ylab = "mean scaled Delta.t * mean scaled Phi",
       main = i.aa)
  for(i.codon in order.plot){
    points(deltat.phi.true[, i.codon], deltat.phi.PM.scaled[, i.codon],
           col = color[u.codon.sorted == u.codon[i.codon]])
  }

  invisible()
} # End of plot.aa.allinone().

