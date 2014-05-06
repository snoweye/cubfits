rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

# Load environment and set data.
source("00-set_env.r")
source(paste(prefix$code.plot, "u0-get_case_main.r", sep = ""))
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)

# Arrange data.
phi.Obs.lim <- range(phi.Obs)
aa.names <- names(reu13.df.obs)
ret.phi.Obs <- prop.bin.roc(reu13.df.obs, phi.Obs)

for(i.case in case.names){
  # Subset of mcmc output.
  fn.in <- paste(prefix$subset, i.case, "_PM.rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }
  load(fn.in)
  fn.in <- paste(prefix$subset, i.case, "_PM_scaling.rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }
  load(fn.in)

  # To adjust to similar range of phi.Obs.
  ret.EPhi <- prop.bin.roc(reu13.df.obs, phi.PM)
  b.PM <- convert.bVec.to.b(b.PM, aa.names)
  predict.roc <- prop.model.roc(b.PM, phi.Obs.lim)

  # Plot bin and model for measurements.
  fn.out <- paste(prefix$plot.diag, "bin_pred_phiObs_", i.case, ".pdf", sep = "")
  pdf(fn.out, width = 12, height = 11)
    nf <- layout(matrix(c(rep(1, 5), 2:21), nrow = 5, ncol = 5, byrow = TRUE),
                 rep(1, 5), c(1, 8, 8, 8, 8), respect = FALSE)
    # Plot title.
    par(mar = c(0, 0, 0, 0))
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.8,
         paste(workflow.name, ", ", get.case.main(i.case, model), sep = ""))
    text(0.5, 0.4, "bin: observed phi")
    par(mar = c(5.1, 4.1, 4.1, 2.1))

    # Plot results.
    for(i.aa in 1:length(aa.names)){
      tmp.obs <- ret.phi.Obs[[i.aa]]
      tmp.roc <- predict.roc[[i.aa]]
      plotbin(tmp.obs, tmp.roc, main = aa.names[i.aa])
    }
    model.label <- c("MCMC Posterior")
    model.lty <- 1
    plot(NULL, NULL, axes = FALSE, main = "", xlab = "", ylab = "",
         xlim = c(0, 1), ylim = c(0, 1))
    legend(0, 0.9, model.label, lty = model.lty, box.lty = 0)
  dev.off()

  # Plot bin and model for predictions.
  fn.out <- paste(prefix$plot.diag, "bin_pred_EPhi_", i.case, ".pdf", sep = "")
  pdf(fn.out, width = 12, height = 11)
    nf <- layout(matrix(c(rep(1, 5), 2:21), nrow = 5, ncol = 5, byrow = TRUE),
                 rep(1, 5), c(1, 8, 8, 8, 8), respect = FALSE)
    # Plot title.
    par(mar = c(0, 0, 0, 0))
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.8,
         paste(workflow.name, ", ", get.case.main(i.case, model), sep = ""))
    text(0.5, 0.4, "bin: posterior mean of Phi")
    par(mar = c(5.1, 4.1, 4.1, 2.1))

    # Plot results.
    for(i.aa in 1:length(aa.names)){
      tmp.obs <- ret.EPhi[[i.aa]]
      tmp.roc <- predict.roc[[i.aa]]
      plotbin(tmp.obs, tmp.roc, main = aa.names[i.aa])
    }
    model.label <- c("MCMC Posterior")
    model.lty <- 1
    plot(NULL, NULL, axes = FALSE, main = "", xlab = "", ylab = "",
         xlim = c(0, 1), ylim = c(0, 1))
    legend(0, 0.9, model.label, lty = model.lty, box.lty = 0)
  dev.off()
}
