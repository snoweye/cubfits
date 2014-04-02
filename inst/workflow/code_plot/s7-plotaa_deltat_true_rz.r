### This is for simulation only, ploting delta.t * phi aganst true values.

rm(list = ls())

library(cubfits)

source("00-set_env.r")
source(paste(prefix$code.plot, "u0-get_case_main.r", sep = ""))
source(paste(prefix$code, "u1-get_negsel.r", sep = ""))
source(paste(prefix$code.plot, "u4-plot_aa_allinone.r", sep = ""))

### Load true.
fn.in <- paste(prefix$data, "simu_true_", model, ".rda", sep = "")
if(file.exists(fn.in)){
  load(fn.in)
} else{
  stop(paste(fn.in, " is not found.", sep = ""))
}

bInit <- convert.b.to.bVec(Eb)

### Load initial.
fn.in <- paste(prefix$data, "/pre_process.rda", sep = "")
load(fn.in)

# Get AA and synonymous codons.
aa.list <- names(reu13.df.obs)
label <- NULL
for(i.codon in aa.list){
  tmp <- sort(unique(reu13.df.obs[[i.codon]]$Codon))
  tmp <- tmp[-length(tmp)]
  label <- c(label, paste(i.codon, tmp, sep = "."))
}

# Get id.deltat.
all.names <- names(bInit)
id.deltat <- grep("(Intercept)", all.names, invert = TRUE)

### Convert true to negsel and delta.t only.
tmp <- get.negsel(bInit, id.deltat, aa.list, label)
bInit <- tmp$b.negsel.PM
label.negsel.true <- tmp$b.negsel.label


# Plot by case.
for(i.case in case.names){
  ### Load subset mcmc run.
  fn.in <- paste(prefix$subset, i.case, ".rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }
  load(fn.in)

  ### Convert unscaled result to negsel and delta.t only.
  tmp <- lapply(1:ncol(b.mcmc),
           function(i.iter){
             tmp <- get.negsel(b.mcmc[, i.iter], id.deltat, aa.list, label)
             tmp$b.negsel.PM
           })
  b.mcmc <- do.call("cbind", tmp)

  ### Load scaled posterior mean (negsel).
  fn.in <- paste(prefix$subset, i.case, "_PM_scaling.rda", sep = "")
  load(fn.in)
  fn.in <- paste(prefix$subset, i.case, "_PM_scaling_rz.rda", sep = "")
  load(fn.in)

  ### Plot.
  t.phi.mcmc <- t(phi.mcmc)
  fn.out <- paste(prefix$plot.AA, "rz_", i.case, "_deltat.phi.pdf", sep = "")
  workflow.name <- paste("RZ ", workflow.name, sep = "")
  # pdf(fn.out, width = 3 * 5, height = 7)
  pdf(fn.out, width = 4, height = 9)
    for(i.aa in aa.list){
      id.label <- grep(paste("^", i.aa, "\\.", sep = ""), label)
      tl.codon <- length(id.label)

      plot.aa.allinone(i.aa, id.label, tl.codon,        # For AA.
                       label.negsel.true,               # For label.
                       bInit, EPhi,                     # For true.
                       b.mcmc, t.phi.mcmc,              # For unscaled results.
                       b.negsel.PM, phi.PM,             # For scaled results.
                       workflow.name, i.case, model)
    }
  dev.off()
}
