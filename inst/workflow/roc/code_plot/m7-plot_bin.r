rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

### Load environment and set data.
source("00-set_env.r")
source(paste(prefix$code.plot, "u0-get_case_main.r", sep = ""))
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)

### Arrange data.
phi.Obs.lim <- range(phi.Obs)
aa.names <- names(reu13.df.obs)

### Load all data.
ret.all <- NULL
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

  ### To adjust to similar range of phi.Obs.
  ret.EPhi <- prop.bin.roc(reu13.df.obs, phi.PM)
  b.PM <- convert.bVec.to.b(b.PM, aa.names)
  predict.roc <- prop.model.roc(b.PM, phi.Obs.lim)

  ret.all[[i.case]] <- list(ret.EPhi = ret.EPhi, predict.roc = predict.roc)
}

### Get possible match cases. wophi fits vs wphi fits on wphi EPhi.
match.case <- rbind(
  c("ad_wophi_pm", "ad_wphi_pm"),
  c("ad_wophi_scuo", "ad_wphi_scuo"),
  c("ad_wophi_true", "ad_wphi_true"),
  c("ad_wphi_wophi_pm", "ad_wphi_pm"),
  c("ad_wphi_wophi_scuo", "ad_wphi_scuo")
)
match.case <- matrix(paste(model, match.case, sep = "_"), ncol = 2)

### Plot matched cases.
for(i.match in 1:nrow(match.case)){
  if(!match.case[i.match, 1] %in% names(ret.all) ||
     !match.case[i.match, 2] %in% names(ret.all)){
    next
  }

  ### Dispatch.
  # ret.EPhi.1 <- ret.all[[match.case[i.match, 1]]]$ret.EPhi
  ret.EPhi.2 <- ret.all[[match.case[i.match, 2]]]$ret.EPhi
  predict.roc.1 <- ret.all[[match.case[i.match, 1]]]$predict.roc
  predict.roc.2 <- ret.all[[match.case[i.match, 2]]]$predict.roc

  ### Fix xlim at log10 scale.
  lim.bin <- range(log10(ret.EPhi.2[[1]]$center))
  lim.model <- range(log10(predict.roc.1[[1]]$center),
                     log10(predict.roc.2[[1]]$center))
  xlim <- c(lim.bin[1] - (lim.bin[2] - lim.bin[1]) / 8,
            lim.bin[2] + (lim.bin[2] - lim.bin[1]) / 8)

  ### Plot bin and model for measurements.
  fn.out <- paste(prefix$plot.match, "bin_", match.case[i.match, 1], "_",
                  match.case[i.match, 2], ".pdf", sep = "")
  pdf(fn.out, width = 16, height = 11)
    mat <- matrix(c(rep(1, 5), 2:21, rep(22, 5)),
                  nrow = 6, ncol = 5, byrow = TRUE)
    mat <- cbind(rep(23, 6), mat, rep(24, 6))
    nf <- layout(mat, c(3, rep(8, 5), 2), c(3, 8, 8, 8, 8, 3), respect = FALSE)
    ### Plot title.
    par(mar = c(0, 0, 0, 0))
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.6,
         paste(workflow.name, ", ", match.case[i.match, 1], " vs ",
               match.case[i.match, 2], sep = ""))
    text(0.5, 0.4, "bin: posterior mean of Phi")
    par(mar = c(0, 0, 0, 0))

    ### Plot results.
    for(i.aa in 1:length(aa.names)){
      tmp.obs <- ret.EPhi.2[[i.aa]]
      tmp.roc <- predict.roc.2[[i.aa]]
      plotbin(tmp.obs, tmp.roc, main = "", xlab = "", ylab = "",
              lty = 2, axes = FALSE, xlim = xlim)
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

      #### Add the first model.
      u.codon <- sort(unique(tmp.obs$codon))
      color <- cubfits:::get.color(u.codon)

      tmp.roc <- predict.roc.1[[i.aa]]
      plotaddmodel(tmp.roc, 1, u.codon, color)
    }

    ### Add label.
    model.label <- paste(model, c("without phi", "with phi"), sep = " ")
    model.lty <- 2:1
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
}
