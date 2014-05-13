### This script plot binning results and predicted curves from multinomial
### logistic regression without measurement errors.

rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

### Load environment and set data.
source("00-set_env.r")
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)
fn.in <- paste(prefix$data, "init_", model, ".rda", sep = "")
load(fn.in)

### Arrange data.
phi.Obs.lim <- range(phi.Obs)
aa.names <- names(reu13.df.obs)
ret.phi.Obs <- prop.bin.roc(reu13.df.obs, phi.Obs)
predict.roc <- prop.model.roc(fitlist, phi.Obs.lim)

### Plot bin and model for measurements.
fn.out <- paste(prefix$plot.diag, "bin_est_phiObs.pdf", sep = "")
pdf(fn.out, width = 16, height = 11)
  mat <- matrix(c(rep(1, 5), 2:21, rep(22, 5)),
                nrow = 6, ncol = 5, byrow = TRUE)
  mat <- cbind(rep(23, 6), mat, rep(24, 6))
  nf <- layout(mat, c(3, rep(8, 5), 2), c(3, 8, 8, 8, 8, 3), respect = FALSE)
  ### Plot title.
  par(mar = c(0, 0, 0, 0))
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.5, workflow.name)
  text(0.5, 0.2, "bin: observed phi")
  par(mar = c(0, 0, 0, 0))

  ### Plot results.
  for(i.aa in 1:length(aa.names)){
    tmp.obs <- ret.phi.Obs[[i.aa]]
    tmp.roc <- predict.roc[[i.aa]]
    plotbin(tmp.obs, tmp.roc, main = "", xlab = "", ylab = "",
            lty = 2, axes = FALSE)
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
  model.label <- c("Logistic Regression")
  model.lty <- 2
  plot(NULL, NULL, axes = FALSE, main = "", xlab = "", ylab = "",
       xlim = c(0, 1), ylim = c(0, 1))
  legend(0.1, 0.8, model.label, lty = model.lty, box.lty = 0)

  ### Plot xlab.
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.5, "Estimated Production Rate (log10)")

  ### Plot ylab.
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.5, "Propotion", srt = 90)
dev.off()
