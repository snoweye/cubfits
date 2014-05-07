### Since this is a simulation, we can plot against true values.

rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

### Preload environment and set data.
source("00-set_env.r")
fn.in <- paste(prefix$data, "simu_true_", model, ".rda", sep = "")
if(file.exists(fn.in)){
  load(fn.in)
} else{
  stop(paste(fn.in, " is not found.", sep = ""))
}
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)

### Arrange data.
aa.names <- names(reu13.df.obs)
phi.Obs <- phi.Obs * phi.scale
phi.Obs.lim <- range(c(phi.Obs, EPhi))

### Compute.
ret.phi.Obs <- prop.bin.roc(reu13.df.obs, phi.Obs)
predict.roc <- prop.model.roc(Eb, phi.Obs.lim)

# Plot bin and model.
fn.out <- paste(prefix$plot.diag, "bin_true_phiObs.pdf", sep = "")
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
    plotbin(tmp.obs, tmp.roc, main = aa.names[i.aa], lty = 3)
  }
  model.label <- c("True Model")
  model.lty <- 3
  plot(NULL, NULL, axes = FALSE, main = "", xlab = "", ylab = "",
       xlim = c(0, 1), ylim = c(0, 1))
  legend(0, 0.9, model.label, lty = model.lty, box.lty = 0)
dev.off()

### Compute.
ret.EPhi <- prop.bin.roc(reu13.df.obs, EPhi)

# Plot bin and model.
fn.out <- paste(prefix$plot.diag, "bin_true_EPhi.pdf", sep = "")
pdf(fn.out, width = 12, height = 11)
  nf <- layout(matrix(c(rep(1, 5), 2:21), nrow = 5, ncol = 5, byrow = TRUE),
               rep(1, 5), c(1, 8, 8, 8, 8), respect = FALSE)
  # Plot title.
  par(mar = c(0, 0, 0, 0))
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.8, workflow.name)
  text(0.5, 0.4, "bin: true phi")
  par(mar = c(5.1, 4.1, 4.1, 2.1))

  # Plot results.
  for(i.aa in 1:length(aa.names)){
    tmp.obs <- ret.EPhi[[i.aa]]
    tmp.roc <- predict.roc[[i.aa]]
    plotbin(tmp.obs, tmp.roc, main = aa.names[i.aa], lty = 3)
  }
  model.label <- c("True Model")
  model.lty <- 3
  plot(NULL, NULL, axes = FALSE, main = "", xlab = "", ylab = "",
       xlim = c(0, 1), ylim = c(0, 1))
  legend(0, 0.9, model.label, lty = model.lty, box.lty = 0)
dev.off()
