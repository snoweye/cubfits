### This is for simulation only, plotting against true values.

rm(list = ls())

library(cubfits, quiet = TRUE)

source("00-set_env.r")
source(paste(prefix$code, "u0-get_case_main.r", sep = ""))

# Load true Phi.
fn.in <- paste(prefix$data, "simu_true_", model, ".rda", sep = "")
if(file.exists(fn.in)){
  load(fn.in)
} else{
  stop(paste(fn.in, " is not found.", sep = ""))
}
Phi <- EPhi

# pre processed phi.Obs.
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)

for(i.case in case.names){
  # subset of mcmc output.
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

  # plot posterior mean.
  fn.out <- paste(prefix$plot.single,
                  "qq_true_", i.case, ".pdf", sep = "")
  pdf(fn.out, width = 5, height = 5)
    # x-axis: predicted, y-axis: true.
    x <- log10(Phi / mean(Phi))
    y <- log10(phi.PM / mean(phi.PM))
    qqplot(x, y,
           xlab = "True Production Rate (log10)",
           ylab = "Predicted Production Rate (log10)",
           main = paste(i.case, " posterior mean", sep = ""),
           pch = 20, cex = 0.6)
    abline(a = 0, b = 1, col = 4, lty = 2)
  dev.off()

  # plot posterior median.
  fn.out <- paste(prefix$plot.single,
                  "qq_true_med_", i.case, ".pdf", sep = "")
  pdf(fn.out, width = 5, height = 5)
    x <- log10(Phi / mean(Phi))
    y <- log10(phi.MED / mean(phi.MED))
    qqplot(x, y,
           xlab = "True Production Rate (log10)",
           ylab = "Predicted Production Rate (log10)",
           main = paste(i.case, " posterior median", sep = ""),
           pch = 20, cex = 0.6)
    abline(a = 0, b = 1, col = 4, lty = 2)
  dev.off()
}

