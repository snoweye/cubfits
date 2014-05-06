### This script plots comparisions across cases.

rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

source("00-set_env.r")

if(length(case.names) < 4){
  stop("Need 4 cases to match with.")
}

# Ordered by "ad_wophi_pm", "ad_wophi_scuo", "ad_wphi_pm", and "ad_wphi_scuo".
phi.mean <- list()
phi.median <- list()
phi.std <- list()
for(i.case in 1:4){
  # Subset of mcmc output.
  fn.in <- paste(prefix$subset, case.names[i.case], "_PM_scaling.rda", sep = "")
  if(!file.exists(fn.in)){
    stop(paste(fn.in, " is not found.", sep = ""))
  }
  load(fn.in)
  phi.mean[[i.case]] <- phi.PM
  phi.median[[i.case]] <- phi.MED
  phi.std[[i.case]] <- phi.STD.log10
}

# Plot posterior mean.
fn.out <- paste(prefix$plot.match, "prxy_pm.pdf", sep = "")
pdf(fn.out, width = 5, height = 5)
  # x-axis: with phi, y-axis: without phi.
  plotprxy(phi.mean[[3]], phi.mean[[1]], weights = 1 / phi.std[[1]],
           xlab = "Production Rate with phi (log10)",
           ylab = "Production Rate without phi (log10)",
           main = "fits vs appr (pm, posterior mean)")
  mtext(workflow.name, line = 3, cex = 0.6)
dev.off()

fn.out <- paste(prefix$plot.match, "prxy_scuo.pdf", sep = "")
pdf(fn.out, width = 5, height = 5)
  # x-axis: with phi, y-axis: without phi.
  plotprxy(phi.mean[[4]], phi.mean[[2]], weights = 1 / phi.std[[2]],
           xlab = "Production Rate with phi (log10)", 
           ylab = "Production Rate without phi (log10)",
           main = "fits vs appr (scuo, posterior mean)")
  mtext(workflow.name, line = 3, cex = 0.6)
dev.off()

# Plot posterior median.
fn.out <- paste(prefix$plot.match, "prxy_med_pm.pdf", sep = "")
pdf(fn.out, width = 5, height = 5)
  # x-axis: with phi, y-axis: without phi.
  plotprxy(phi.median[[3]], phi.median[[1]], weights = 1 / phi.std[[1]],
           xlab = "Production Rate with phi (log10)",
           ylab = "Production Rate without phi (log10)",
           main = "fits vs appr (pm, posterior median)")
  mtext(workflow.name, line = 3, cex = 0.6)
dev.off()

fn.out <- paste(prefix$plot.match, "prxy_med_scuo.pdf", sep = "")
pdf(fn.out, width = 5, height = 5)
  # x-axis: with phi, y-axis: without phi.
  plotprxy(phi.median[[4]], phi.median[[2]], weights = 1 / phi.std[[2]],
           xlab = "Production Rate with phi (log10)", 
           ylab = "Production Rate without phi (log10)",
           main = "fits vs appr (scuo, posterior median)")
  mtext(workflow.name, line = 3, cex = 0.6)
dev.off()

