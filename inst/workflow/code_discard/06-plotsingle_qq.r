### This script plot X-Y protein production rates.

rm(list = ls())

library(cubfits, quiet = TRUE)

source("00-set_env.r")
source(paste(prefix$code, "u0-get_case_main.r", sep = ""))

# Pre processed phi.Obs.
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)

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

  # Plot posterior mean.
  fn.out <- paste(prefix$plot.single,
                  "qq_", i.case, ".pdf", sep = "")
  pdf(fn.out, width = 5, height = 5)
    # x-axis: predicted, y-axis: observed.
    x <- log10(phi.PM / mean(phi.PM))
    y <- log10(phi.Obs / mean(phi.Obs))
    qqplot(x, y,
           xlab = "Predicted Production Rate (log10)",
           ylab = "Observed Production Rate (log10)",
           main = paste(i.case, " posterior mean", sep = ""),
           pch = 20, cex = 0.6)
    abline(a = 0, b = 1, col = 4, lty = 2)
  dev.off()

  # Plot posterior median.
  fn.out <- paste(prefix$plot.single,
                  "qq_med_", i.case, ".pdf", sep = "")
  pdf(fn.out, width = 5, height = 5)
    # x-axis: predicted, y-axis: observed.
    x <- log10(phi.MED / mean(phi.MED))
    y <- log10(phi.Obs / mean(phi.Obs))
    qqplot(x, y,
           xlab = "Predicted Production Rate (log10)",
           ylab = "Observed Production Rate (log10)",
           main = paste(i.case, " posterior median", sep = ""),
           pch = 20, cex = 0.6)
    abline(a = 0, b = 1, col = 4, lty = 2)
  dev.off()
}

