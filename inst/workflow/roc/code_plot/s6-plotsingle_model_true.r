### This script plots binning and model predictions from MCMC with measurement
### errors.

rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

### Load environment and set data.
source("00-set_env.r")
source(paste(prefix$code.plot, "u0-get_case_main.r", sep = ""))
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)
fn.in <- paste(prefix$data, "init_", model, ".rda", sep = "")
load(fn.in)
fn.in <- paste(prefix$data, "simu_true_", model, ".rda", sep = "")
load(fn.in)
EPhi.true <- EPhi

### Arrange data.
EPhi.true.lim <- range(EPhi.true)
aa.names <- names(reu13.df.obs)
ret.EPhi.true <- prop.bin.roc(reu13.df.obs, EPhi.true)
noerror.roc <- prop.model.roc(fitlist, EPhi.true.lim)

tmp <- convert.b.to.bVec(fitlist)
id.slop <- grep("Delta.t", names(tmp), invert = TRUE)

### Since Eb is not generated in scale of mean 1, but phi.Obs was already
### scaled in mean 1. I have to scale b.true accordingly.
b.true <- convert.b.to.bVec(Eb)
# b.true[id.slop] <- b.true[id.slop] * phi.scale
 b.true[id.slop] <- b.true[id.slop] * mean(EPhi) 
b.true <- convert.bVec.to.b(b.true, aa.names)
true.roc <- prop.model.roc(b.true, EPhi.true.lim)


### Load each chain.
for(i.case in case.names){
  ### Subset of mcmc output.
  fn.in <- paste(prefix$subset, i.case, "_PM.rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }
  load(fn.in)

  ### Subset of mcmc output with scaling.
  fn.in <- paste(prefix$subset, i.case, "_PM_scaling.rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }
  load(fn.in)

  b.PM <- convert.bVec.to.b(b.PM, aa.names)
  predict.roc <- prop.model.roc(b.PM, EPhi.true.lim)

  ### Plot bin and model for measurements.
  fn.out <- paste(prefix$plot.single, "bin_merge_true_",
                  i.case, ".pdf", sep = "")
  pdf(fn.out, width = 16, height = 11)
    mat <- matrix(c(rep(1, 5), 2:21, rep(22, 5)),
                  nrow = 6, ncol = 5, byrow = TRUE)
    mat <- cbind(rep(23, 6), mat, rep(24, 6))
    nf <- layout(mat, c(3, rep(8, 5), 2), c(3, 8, 8, 8, 8, 3), respect = FALSE)
    ### Plot title.
    par(mar = c(0, 0, 0, 0))
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.5,
         paste(workflow.name, ", ", get.case.main(i.case, model), sep = ""))
    text(0.5, 0.2, "bin: true Phi")
    par(mar = c(0, 0, 0, 0))

    ### Plot results.
    for(i.aa in 1:length(aa.names)){
      tmp.obs <- ret.EPhi.true[[i.aa]]
      tmp.roc <- predict.roc[[i.aa]]
      plotbin(tmp.obs, tmp.roc, main = "", xlab = "", ylab = "",
              lty = 1, axes = FALSE)
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

      ### Add true model if it is available.
      u.codon <- sort(unique(tmp.obs$codon))
      color <- cubfits:::get.color(u.codon)

      if(length(grep("_wophi_", i.case)) == 0){  ### wphi case, add regression.
        tmp.roc <- noerror.roc[[i.aa]]
        plotaddmodel(tmp.roc, 2, u.codon, color)
      }

      if(exists("Eb")){
        tmp.roc <- true.roc[[i.aa]]
        plotaddmodel(tmp.roc, 3, u.codon, color)
      }
    }

    ### Add label.
    model.label <- "MCMC Posterior"
    model.lty <- 1
    if(length(grep("_wophi_", i.case)) == 0){  ### wphi case, add regression.
      model.label <- c(model.label, "Logistic Regression")
      model.lty <- c(model.lty, 2)
    }
    if(exists("Eb")){
      model.label <- c(model.label, "True Model")
      model.lty <- c(model.lty, 3)
    }
    plot(NULL, NULL, axes = FALSE, main = "", xlab = "", ylab = "",
         xlim = c(0, 1), ylim = c(0, 1))
    legend(0.1, 0.8, model.label, lty = model.lty, box.lty = 0)

    ### Plot xlab.
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    if(exists("Eb")){
      text(0.5, 0.5, "True Production Rate (log10)")
    } else{
      text(0.5, 0.5, "Estimated Production Rate (log10)")
    }

    ### Plot ylab.
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.5, "Propotion", srt = 90)
  dev.off()
}
