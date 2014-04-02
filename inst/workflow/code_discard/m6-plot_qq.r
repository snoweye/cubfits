### This script plots comparisions across cases.

rm(list = ls())

library(cubfits, quiet = TRUE)

source("00-set_env.r")

if(length(case.names) < 4){
  stop("Need 4 cases to match with.")
}

# Ordered by "ad_wophi_pm", "ad_wophi_scuo", "ad_wphi_pm", and "ad_wphi_scuo".
phi.mean <- list()
phi.median <- list()
for(i.case in 1:4){
  # Subset of mcmc output.
  fn.in <- paste(prefix$subset, case.names[i.case], "_PM_scaling.rda", sep = "")
  if(!file.exists(fn.in)){
    stop(paste(fn.in, " is not found.", sep = ""))
  }
  load(fn.in)
  phi.mean[[i.case]] <- phi.PM
  phi.median[[i.case]] <- phi.MED
}

# Plot posterior mean.
fn.out <- paste(prefix$plot.match, "qq_pm.pdf", sep = "")
pdf(fn.out, width = 5, height = 5)
  # x-axis: with phi, y-axis: without phi.
  x <- log10(phi.mean[[3]] / mean(phi.mean[[3]]))
  y <- log10(phi.mean[[1]] / mean(phi.mean[[1]]))
  qqplot(x, y,
         xlab = "Production Rate with phi (log10)",
         ylab = "Production Rate without phi (log10)",
         main = "fits vs appr (pm, posterior mean)",
         pch = 20, cex = 0.6)
  abline(a = 0, b = 1, col = 4, lty = 2)
dev.off()

fn.out <- paste(prefix$plot.match, "qq_scuo.pdf", sep = "")
pdf(fn.out, width = 5, height = 5)
  # x-axis: with phi, y-axis: without phi.
  x <- log10(phi.mean[[4]] / mean(phi.mean[[4]]))
  y <- log10(phi.mean[[2]] / mean(phi.mean[[2]]))
  qqplot(x, y,
         xlab = "Production Rate with phi (log10)", 
         ylab = "Production Rate without phi (log10)",
         main = "fits vs appr (scuo, posterior mean)",
         pch = 20, cex = 0.6)
  abline(a = 0, b = 1, col = 4, lty = 2)
dev.off()

# Plot posterior median.
fn.out <- paste(prefix$plot.match, "qq_med_pm.pdf", sep = "")
pdf(fn.out, width = 5, height = 5)
  # x-axis: with phi, y-axis: without phi.
  x <- log10(phi.median[[3]] / mean(phi.median[[3]]))
  y <- log10(phi.median[[1]] / mean(phi.median[[1]]))
  qqplot(x, y,
         xlab = "Production Rate with phi (log10)",
         ylab = "Production Rate without phi (log10)",
         main = "fits vs appr (pm, posterior median)",
         pch = 20, cex = 0.6)
  abline(a = 0, b = 1, col = 4, lty = 2)
dev.off()

fn.out <- paste(prefix$plot.match, "qq_med_scuo.pdf", sep = "")
pdf(fn.out, width = 5, height = 5)
  # x-axis: with phi, y-axis: without phi.
  x <- log10(phi.median[[4]] / mean(phi.median[[4]]))
  y <- log10(phi.median[[2]] / mean(phi.median[[2]]))
  qqplot(x, y,
         xlab = "Production Rate with phi (log10)", 
         ylab = "Production Rate without phi (log10)",
         main = "fits vs appr (scuo, posterior median)",
         pch = 20, cex = 0.6)
  abline(a = 0, b = 1, col = 4, lty = 2)
dev.off()

