
### This plots one amino acid for ROC or NSE model.
plotbin <- function(ret.bin, ret.model = NULL, main = NULL,
    xlab = "Production Rate (log10)", ylab = "Proportion",
    xlim = NULL, lty = 1, x.log10 = TRUE, stderr = FALSE, ...){
  if(x.log10){
    ret.bin$center <- log10(ret.bin$center)
    if(!is.null(ret.model)){
      ret.model$center <- log10(ret.model$center)
    }
  }

  ### Observed dots and whiskers.
  if(is.null(xlim)){
    if(!is.null(ret.model)){
      # x.lim <- range(c(ret.bin$center, ret.model$center))
      lim.bin <- range(ret.bin$center)
      lim.model <- range(ret.model$center)
      x.lim <- c((lim.bin[1] + lim.model[1]) / 2,
                 max(lim.bin[2], lim.model[2]))
    } else{
      x.lim <- range(ret.bin$center)
    }
  } else{
    x.lim <- xlim
  }
  y.lim <- c(0, 1) + c(-0.05, 0.05)

  u.codon <- sort(unique(as.character(ret.bin$codon)))
  u.center <- unique(ret.bin$center)
  color <- get.color(u.codon)
  ### Reorder R for better legend.
  if(all(u.codon %in% .CF.GV$synonymous.codon$R)){
    u.codon <- u.codon[c(3:6, 1:2)]
    color <- color[c(3:6, 1:2)]
  }

  ### Make an empty plot
  plot(NULL, NULL, xlim = x.lim, ylim = y.lim,
       main = main, xlab = xlab, ylab = ylab, ...)

  ### Add observed dots for means and whiskers for standard deviations.
  for(i.codon in 1:length(u.codon)){
    ret.bin.codon <- ret.bin[ret.bin$codon == u.codon[i.codon],]
    points(ret.bin.codon$center,
           ret.bin.codon$freq.mean,
           pch = 19, col = color[i.codon], cex = 0.5)
    for(i.center in 1:nrow(ret.bin.codon)){
      if(!is.na(ret.bin.codon$freq.std[i.center])){
        if(stderr){
          freq.bar <- ret.bin.codon$freq.stderr[i.center]
        } else{
          freq.bar <- ret.bin.codon$freq.std[i.center]
        }
        lines(
          list(x = rep(ret.bin.codon$center[i.center], 2),
               y = ret.bin.codon$freq.mean[i.center] + c(-1, 1) * freq.bar),
          col = color[i.codon], lwd = 0.8)
      }
    }
  }

  ### Add modeled lines.
  u.codon.star <- u.codon
  if(!is.null(ret.model)){
    plotaddmodel(ret.model, lty, u.codon, color, x.log10 = FALSE)

    ### Add focal codon.
    u.codon.star <- attr(ret.model, "u.codon.star")
    if(!is.null(u.codon.star)){
      if(all(u.codon %in% .CF.GV$synonymous.codon$R)){
        u.codon.star <- u.codon.star[c(3:6, 1:2)]
      }
    }
  }

  ### Add legends.
  legend(x.lim[1], y.lim[2], u.codon.star, col = color, 
         box.lty = 0, lty = 1, pch = 19, cex = 0.8)
} # End of plotbin().

