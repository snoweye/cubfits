### For plotting b corr.

my.range <- function(x){
  lim <- range(x)
  lim + c(-1, 1) * (lim[2] - lim[1]) * 0.05
} # End of my.range().

plot.b.corr <- function(x, y, label, x.ci = NULL, y.ci = NULL,
    xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, main = NULL,
    workflow.name = NULL, add.lm = FALSE, add.ci = FALSE){
  if(!is.null(x.ci)){
    x.ci <- matrix(x.ci, ncol = 2)
  }
  if(!is.null(y.ci)){
    y.ci <- matrix(y.ci, ncol = 2)
  }
  if(is.null(xlim)){
    xlim <- my.range(x)
  }
  if(is.null(ylim)){
    ylim <- my.range(y)
  }

  ### main part.
  plot(NULL, NULL, xlim = xlim, ylim = ylim,
       xlab = xlab, ylab = ylab, main = main)
  mtext(workflow.name, line = 3, cex = 0.6)

  ### Add lm.
  if(add.lm){
    m.1 <- try(lm(y ~ x), silent = TRUE)
    if(class(m.1) != "try-error"){
      a <- m.1$coef[1]
      b <- m.1$coef[2]
      R2 <- summary(m.1)$r.squared
      abline(a = a, b = b, col = 2)

      width <- xlim[2] - xlim[1]
      height <- ylim[2] - ylim[1]
      text(xlim[1] - width * 0.04, ylim[2] - height * 0.05,
           parse(text = paste("y == ", a, " + ", b, " * x", sep = "")),
           pos = 4, cex = 0.5)

      text(xlim[1] - width * 0.04, ylim[2] - height * 0.10,
           parse(text = paste("R^2 == ",
                              sprintf("%.4f", R2), sep = "")),
           pos = 4, cex = 0.5)

      ### Add 95% CI.
      if(add.ci){
        text(xlim[1] - width * 0.04, ylim[2] - height * 0.15,
             "95% CI ",
             pos = 4, cex = 0.5)

        tmp <- summary(m.1)

        ### For intercept.
        CI.95 <- tmp$coefficients[1, 1] +
                 qt(c(0.025, 0.975), tmp$df[2]) * tmp$coefficients[1, 2]
        text(xlim[1] - width * 0.02, ylim[2] - height * 0.20,
             paste("intercept: ",
                   sprintf("(%.4f, %.4f)", CI.95[1], CI.95[2]), sep = ""),
             pos = 4, cex = 0.5)

        ### For slop.
        CI.95 <- tmp$coefficients[2, 1] +
                 qt(c(0.025, 0.975), tmp$df[2]) * tmp$coefficients[2, 2]
        text(xlim[1] - width * 0.02, ylim[2] - height * 0.25,
             paste("slop: ",
                   sprintf("(%.4f, %.4f)", CI.95[1], CI.95[2]), sep = ""),
             pos = 4, cex = 0.5)

        ### For LRT.
        # lrt.pvalue <- 0
        # text(xlim[1] - width * 0.02, ylim[2] - height * 0.30,
        #      paste("LRT p-value: "
        #            sprintf("%.4f", lrt.pvalue), sep = ""),
        #      pos = 4, cex = 0.5)
      }
    }
  }

  ### Add one-to-one.
  abline(a = 0, b = 1, col = 4, lty = 2)

  ### Add main lines.
  for(i in 1:length(x)){
    if(!is.null(x.ci)){
      lines(x = x.ci[i,], y = rep(y[i], 2), col = 1)
    }
    if(!is.null(y.ci)){
      lines(x = rep(x[i], 2), y = y.ci[i,], col = 1)
    }
  }
  if(is.null(x.ci) || is.null(y.ci)){
    points(x, y, pch = 20, cex = 0.5, col = 2)
  }

  ### Add label.
  x.split <- xlim[1] + (xlim[2] - xlim[1]) / 2
  tmp.id <- order(x)
  tmp.x <- x[tmp.id]
  tmp.y <- y[tmp.id]
  tmp.label <- label[tmp.id]
  x.adj <- (xlim[2] - xlim[1]) * 0.1 *
           rep(c(1, -1), c(sum(tmp.x <= x.split), sum(tmp.x > x.split))) *
           ((1:length(tmp.x) - 1) %% 3 + 1)
  tmp.x <- tmp.x + x.adj
  text(tmp.x, tmp.y, labels = tmp.label, cex = 0.4)
} # End of plot.b.corr().
