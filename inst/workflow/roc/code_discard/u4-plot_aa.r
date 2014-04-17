### This is for simulation only, ploting delta.t * phi aganst true values.

plot.aa <- function(i.aa, id.label, tl.codon,    # For AA.
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
  # nf <- layout(matrix(c(rep(1, tl.codon), 1 + 1:(tl.codon * 2)),
  #                     nrow = 3, ncol = tl.codon, byrow = TRUE),
  #              rep(1, tl.codon), c(1, 8, 8), respect = FALSE)
  nf <- layout(rbind(rep(1, 5), matrix(1 + 1:(5 * 2), nrow = 2)),
               rep(1, 5), c(1, 8, 8), respect = FALSE)

  # Plot title.
  par(mar = c(0, 0, 0, 0))
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.5,
       paste(workflow.name, ", ", get.case.main(i.case, model),
             sep = ""))
  par(mar = c(5.1, 4.1, 4.1, 2.1))

  i.plot.count <- 0
  for(i.codon in 1:tl.codon){
    # Plot top row.
    lim <- range(range(deltat.phi.true[, i.codon]),
                 range(deltat.phi.PM[, i.codon]))
    plotprxy(deltat.phi.true[, i.codon], deltat.phi.PM[, i.codon],
             log10.x = FALSE, log10.y = FALSE,
             xlim = lim, ylim = lim,
             xlab = "true Delta.t * true Phi", ylab = "mean (Delta.t * Phi)",
             main = label.negsel.true[id.label[i.codon]])
    title(sub = sprintf("mean Delta.t=%.4f, true Delta.t=%.4f",
                        b.AA.PM[i.codon], b.AA.true[i.codon]))

    # Plot bottom row.
    lim <- range(range(deltat.phi.true[, i.codon]),
                 range(deltat.phi.PM.scaled[, i.codon]))
    plotprxy(deltat.phi.true[, i.codon], deltat.phi.PM.scaled[, i.codon],
             log10.x = FALSE, log10.y = FALSE,
             xlim = lim, ylim = lim,
             xlab = "true Delta.t * true Phi",
             ylab = "mean scaled Delta.t * mean scaled Phi",
             main = label.negsel.true[id.label[i.codon]])
    title(sub = sprintf("scaled Delta.t=%.4f, true Delta.t=%.4f",
                        b.AA.PM.scaled[i.codon], b.AA.true[i.codon]))
  }
} # End of plot.aa().

