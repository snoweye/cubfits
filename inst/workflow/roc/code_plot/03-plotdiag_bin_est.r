### This script plot binning results and predicted curves from multinomial
### logistic regression without measurement errors.

rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

# Load environment and set data.
source("00-set_env.r")
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)
fn.in <- paste(prefix$data, "init_", model, ".rda", sep = "")
load(fn.in)

# Arrange data.
phi.Obs.lim <- range(phi.Obs)
aa.names <- names(reu13.df.obs)
ret.phi.Obs <- prop.bin.roc(reu13.df.obs, phi.Obs)
predict.roc <- prop.model.roc(fitlist, phi.Obs.lim)

# Plot bin and model for measurements.
fn.out <- paste(prefix$plot.diag, "bin_est_phiObs.pdf", sep = "")
pdf(fn.out, width = 12, height = 11)
  nf <- layout(matrix(c(rep(1, 5), 2:21), nrow = 5, ncol = 5, byrow = TRUE),
               rep(1, 5), c(1, 8, 8, 8, 8), respect = FALSE)
  # Plot title.
  par(mar = c(0, 0, 0, 0))
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.8, workflow.name)
  text(0.5, 0.4, "bin: observed phi")
  par(mar = c(5.1, 4.1, 4.1, 2.1))

  # Plot results.
  for(i.aa in 1:length(aa.names)){
    tmp.obs <- ret.phi.Obs[[i.aa]]
    tmp.roc <- predict.roc[[i.aa]]
    plotbin(tmp.obs, tmp.roc, main = aa.names[i.aa], lty = 2)
  }
  model.label <- c("Logistic Regression")
  model.lty <- 2
  plot(NULL, NULL, axes = FALSE, main = "", xlab = "", ylab = "",
       xlim = c(0, 1), ylim = c(0, 1))
  legend(0, 0.9, model.label, lty = model.lty, box.lty = 0)
dev.off()

