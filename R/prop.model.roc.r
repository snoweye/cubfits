# Wrapper of find.prop.model.roc().
prop.model.roc <- function(bInit, phi.Obs.lim = c(0.01, 10), phi.Obs.scale = 1,
    nclass = 40, x.log10 = TRUE){
  aa.list <- names(bInit)
  if("Z" %in% aa.list){
    synonymous.codon <- .CF.GV$synonymous.codon.split[aa.list]
  } else{
    synonymous.codon <- .CF.GV$synonymous.codon[aa.list]
  }

  # Compute mutation and elong
  for(i.aa in aa.list){
    bInit[[i.aa]]$u.codon <- synonymous.codon[[i.aa]]
  }

  # For SCUO w/o phi.Obs simulations, the prior fixes the mean of phi.Obs at 1, but
  # lets std drift only. Therefore, E[Phi] may or may not be distributed in
  # the same range of phi.Obs.
  # Convert phi.Obs.lim to SCUO's scale, and evaluate the proportional frequences
  # of synonymous codons, then convert back to phi.Obs's scale (for plotting)
  # inside the function find.prop.model.roc().
  phi.Obs.lim <- phi.Obs.lim / phi.Obs.scale

  # For better plottling
  if(x.log10){
    phi.bin <- 10^seq(log10(phi.Obs.lim[1]), log10(phi.Obs.lim[2]), length = nclass)
  } else{
    phi.bin <- seq(phi.Obs.lim[1], phi.Obs.lim[2], length = nclass)
  }

  # Call find.prop.model.roc().
  ret <- find.prop.model.roc(bInit, phi.bin, phi.Obs.scale)
  ret
} # End of prop.model.roc().


### Summarize by amino acid for bInint could be from MCMC outputs.
find.prop.model.roc <- function(bInit, phi.bin, phi.Obs.scale = 1){
  u.aa <- unique(names(bInit))
  x <- cbind(1, phi.bin)

  predict.roc <- list()
  for(aa in u.aa){
    exponent <- x %*% bInit[[aa]]$coef.mat
    scodon.prob <- my.inverse.mlogit(exponent)
    predict.roc[[aa]] <- cbind(scodon.prob, phi.bin * phi.Obs.scale)
    predict.roc[[aa]] <- as.data.frame(predict.roc[[aa]],
                                       stringsAsFactors = FALSE)
    colnames(predict.roc[[aa]]) <- c(bInit[[aa]]$u.codon, "center")
  }

  names(predict.roc) <- u.aa
  predict.roc
} # End of find.prop.model.roc().

