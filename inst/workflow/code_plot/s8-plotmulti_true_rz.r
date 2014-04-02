### This is for simulation only, ploting correlation aganst true values.

rm(list = ls())

library(cubfits, quiet = TRUE)

source("00-set_env.r")
source(paste(prefix$code.plot, "u0-get_case_main.r", sep = ""))
source(paste(prefix$code, "u1-get_negsel.r", sep = ""))
source(paste(prefix$code.plot, "u2-plot_b_corr.r", sep = ""))
source(paste(prefix$code.plot, "u5-new_page.r", sep = ""))

# Load true Phi.
fn.in <- paste(prefix$data, "simu_true_", model, ".rda", sep = "")
if(file.exists(fn.in)){
  load(fn.in)
} else{
  stop(paste(fn.in, " is not found.", sep = ""))
}

bInit <- convert.b.to.bVec(Eb)

# Load data.
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)

# Get AA and synonymous codons.
aa.list <- names(reu13.df.obs)
label <- NULL
for(i.aa in aa.list){
  tmp <- sort(unique(reu13.df.obs[[i.aa]]$Codon))
  tmp <- tmp[-length(tmp)]
  label <- c(label, paste(i.aa, tmp, sep = "."))
}

# Get true values.
all.names <- names(bInit)
id.slop <- grep("(Intercept)", all.names, invert = TRUE)
scale.EPhi <- mean(EPhi)
bInit[id.slop] <- bInit[id.slop] * scale.EPhi

id.slop <- grep("(Intercept)", all.names, invert = TRUE)
bInit.negsel <- get.negsel(bInit, id.slop, aa.list, label)


# Plot by case.
for(i.case in case.names){
  # Check files.
  fn.in <- paste(prefix$subset, i.case, ".rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }
  fn.in <- paste(prefix$subset, i.case, "_PM.rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }


  # Subset of mcmc output.
  fn.in <- paste(prefix$subset, i.case, ".rda", sep = "")
  load(fn.in)

  # Subset of mcmc output.
  fn.in <- paste(prefix$subset, i.case, "_PM.rda", sep = "")
  load(fn.in)
  fn.in <- paste(prefix$subset, i.case, "_PM_scaling.rda", sep = "")
  load(fn.in)


  # Set layout.
  fn.out <- paste(prefix$plot.multi, "rz_", i.case, ".pdf", sep = "")
  workflow.name <- paste("RZ ", workflow.name, sep = "")
  pdf(fn.out, width = 6, height = 10)
#New page
    new.page(workflow.name, i.case, model)

    # Plot org Delta.t.
    x <- bInit.negsel$b.negsel.PM
    y <- b.negsel.PM
    y.ci <- b.negsel.ci.PM
    x.label <- bInit.negsel$b.negsel.label
    plot.b.corr(x, y, x.label, y.ci = y.ci,
                xlab = "True", ylab = "Estimated",
                main = "Delta.t", add.lm = TRUE)

    # Plot RZ Delta.t.
    fn.in <- paste(prefix$subset, i.case, "_PM_scaling_rz.rda", sep = "")
    load(fn.in)
    y <- b.negsel.PM
    y.ci <- b.negsel.ci.PM
    plot.b.corr(x, y, x.label, y.ci = y.ci,
                xlab = "True", ylab = "Estimated",
                main = "Delta.t RZ", add.lm = TRUE)

  dev.off()
}
