### This script plots binning and model predictions from MCMC with measurement
### errors.

rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

# Load environment and set data.
source("00-set_env.r")
source(paste(prefix$code.plot, "u0-get_case_main.r", sep = ""))
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)
fn.in <- paste(prefix$data, "init_", model, ".rda", sep = "")
load(fn.in)
fn.in <- paste(prefix$data, "simu_true_", model, ".rda", sep = "")
load(fn.in)
EPhi.true <- EPhi

# Arrange data.
EPhi.true.lim <- range(EPhi.true)
aa.names <- names(reu13.df.obs)
ret.EPhi.true <- prop.bin.roc(reu13.df.obs, EPhi.true)
noerror.roc <- prop.model.roc(fitlist, EPhi.true.lim)

tmp <- convert.b.to.bVec(fitlist)
id.slop <- grep("omega", names(tmp))

### Since Eb is not generated in scale of mean 1, but phi.Obs was already
### scaled in mean 1. I have to scale b.true accordingly.
b.true <- convert.b.to.bVec(Eb)
# b.true[id] <- b.true[id] * phi.scale
 b.true[id.slop] <- b.true[id.slop] * mean(EPhi) 
b.true <- convert.bVec.to.b(b.true, aa.names)
true.roc <- prop.model.roc(b.true, EPhi.true.lim)


# Load each chain.
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

  b.PM <- convert.bVec.to.b(b.PM, aa.names)
  predict.roc <- prop.model.roc(b.PM, EPhi.true.lim)

  # Plot bin and model for measurements.
  fn.out <- paste(prefix$plot.single, "bin_merge_true_",
                  i.case, ".pdf", sep = "")
  pdf(fn.out, width = 12, height = 11)
    nf <- layout(matrix(c(rep(1, 5), 2:21), nrow = 5, ncol = 5, byrow = TRUE),
                 rep(1, 5), c(1, 8, 8, 8, 8), respect = FALSE)
    # Plot title.
    par(mar = c(0, 0, 0, 0))
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.8,
         paste(workflow.name, ", ", get.case.main(i.case, model), sep = ""))
    text(0.5, 0.4, "bin: true Phi")
    par(mar = c(5.1, 4.1, 4.1, 2.1))

    # Plot results.
    for(i.aa in 1:length(aa.names)){
      tmp.obs <- ret.EPhi.true[[i.aa]]
      tmp.roc <- predict.roc[[i.aa]]
      plotbin(tmp.obs, tmp.roc, main = aa.names[i.aa])

      u.codon <- sort(unique(tmp.obs$codon))
      color <- cubfits:::get.color(u.codon)

      if(length(grep("_wophi_", i.case)) == 0){  # wphi case, add regression.
        tmp.roc <- noerror.roc[[i.aa]]
        plotaddmodel(tmp.roc, 2, u.codon, color)
      }

      if(exists("Eb")){
        tmp.roc <- true.roc[[i.aa]]
        plotaddmodel(tmp.roc, 3, u.codon, color)
      }
    }

    # Add label.
    model.label <- "MCMC Posterior"
    model.lty <- 1
    if(length(grep("_wophi_", i.case)) == 0){  # wphi case, add regression.
      model.label <- c(model.label, "Logistic Regression")
      model.lty <- c(model.lty, 2)
    }
    if(exists("Eb")){
      model.label <- c(model.label, "True Model")
      model.lty <- c(model.lty, 3)
    }
    plot(NULL, NULL, axes = FALSE, main = "", xlab = "", ylab = "",
         xlim = c(0, 1), ylim = c(0, 1))
    legend(0, 0.9, model.label, lty = model.lty, box.lty = 0)
  dev.off()
}
