### This script plot correlation of omega with adjusted reference codon
### such that all selection coefficients are negative.

rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

source("00-set_env.r")
source(paste(prefix$code.plot, "u2-plot_b_corr.r", sep = ""))

if(length(case.names) < 4){
  stop("Need 4 cases to match with.")
}

# Load data.
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)

# Ordered by "ad_wophi_pm", "ad_wophi_scuo", "ad_wphi_pm", and "ad_wphi_scuo".
b.ci.org <- list()
b.mean.org <- list()
label.org <- list()
b.ci <- list()
b.mean <- list()
label <- list()
for(i.case in 1:4){
  # Subset of mcmc output.
  fn.in <- paste(prefix$subset, case.names[i.case], "_PM.rda", sep = "")
  if(!file.exists(fn.in)){
    stop(paste(fn.in, " is not found.", sep = ""))
  }
  load(fn.in)

  b.ci.org[[i.case]] <- b.negsel.ci.PM
  b.mean.org[[i.case]] <- b.negsel.PM
  label.org[[i.case]] <- b.negsel.label

  # Subset of mcmc output with scaling.
  fn.in <- paste(prefix$subset, case.names[i.case], "_PM_scaling.rda", sep = "")
  load(fn.in)

  b.ci[[i.case]] <- b.negsel.ci.PM
  b.mean[[i.case]] <- b.negsel.PM
  label[[i.case]] <- b.negsel.label
}


# Plot omega.
x.pm <- b.mean[[3]]
y.pm <- b.mean[[1]]
x.pm.label <- label[[3]]
x.scuo <- b.mean[[4]]
y.scuo <- b.mean[[2]]
x.scuo.label <- label[[4]]
xlim <- my.range(c(x.pm, x.scuo))
ylim <- my.range(c(y.pm, y.scuo))
x.pm.ci <- b.ci[[3]]
y.pm.ci <- b.ci[[1]]
x.scuo.ci <- b.ci[[4]]
y.scuo.ci <- b.ci[[2]]

# Convert selection to negative value by changing the relative base. 
# Note that this can fail in some inconsistent cases.
if(any((x.pm < 0 & y.pm > 0) | (x.pm > 0 & y.pm < 0) |
        label[[1]] != label[[3]])){
  stop("Inconsistent cases (PM).")
}
if(any((x.scuo < 0 & y.scuo > 0) | (x.scuo > 0 & y.scuo < 0) |
        label[[2]] != label[[4]])){
  stop("Inconsistent cases (SCUO).")
}

# Plot omega.
fn.out <- paste(prefix$plot.match, "corr_negsel_omega_pm.pdf", sep = "")
pdf(fn.out, width = 5, height = 5)
  plot.b.corr(x.pm, y.pm, x.pm.label,
              x.ci = x.pm.ci, y.ci = y.pm.ci,
              xlim = xlim, ylim = ylim,
              xlab = "omega with phi", ylab = "omega without phi",
              main = "roc_ad_pm", workflow.name = workflow.name)
dev.off()

fn.out <- paste(prefix$plot.match, "corr_negsel_omega_scuo.pdf", sep = "")
pdf(fn.out, width = 5, height = 5)
  plot.b.corr(x.scuo, y.scuo, x.scuo.label,
              x.ci = x.scuo.ci, y.ci = y.scuo.ci,
              xlim = xlim, ylim = ylim,
              xlab = "omega with phi", ylab = "omega without phi",
              main = "roc_ad_scuo", workflow.name = workflow.name)
dev.off()


# Plot omega original (no scaling by x to mean = 1).
x.pm <- b.mean.org[[3]]
y.pm <- b.mean.org[[1]]
x.pm.label <- label.org[[3]]
x.scuo <- b.mean.org[[4]]
y.scuo <- b.mean.org[[2]]
x.scuo.label <- label.org[[4]]
xlim <- my.range(c(x.pm, x.scuo))
ylim <- my.range(c(y.pm, y.scuo))
x.pm.ci <- b.ci.org[[3]]
y.pm.ci <- b.ci.org[[1]]
x.scuo.ci <- b.ci.org[[4]]
y.scuo.ci <- b.ci.org[[2]]

# Convert selection to negative value by changing the relative base. 
# Note that this can fail in some inconsistent cases.
if(any((x.pm < 0 & y.pm > 0) | (x.pm > 0 & y.pm < 0) |
        label[[1]] != label[[3]])){
  stop("Inconsistent cases (PM).")
}
if(any((x.scuo < 0 & y.scuo > 0) | (x.scuo > 0 & y.scuo < 0) |
        label[[2]] != label[[4]])){
  stop("Inconsistent cases (SCUO).")
}

fn.out <- paste(prefix$plot.match, "corr_nonscale_negsel_omega_pm.pdf", sep = "")
pdf(fn.out, width = 5, height = 5)
  plot.b.corr(x.pm, y.pm, x.pm.label,
              x.ci = x.pm.ci, y.ci = y.pm.ci,
              xlim = xlim, ylim = ylim,
              xlab = "omega with phi", ylab = "omega without phi",
              main = "roc_ad_pm_nonscale", workflow.name = workflow.name)
dev.off()

fn.out <- paste(prefix$plot.match, "corr_nonscale_negsel_omega_scuo.pdf", sep = "")
pdf(fn.out, width = 5, height = 5)
  plot.b.corr(x.scuo, y.scuo, x.scuo.label,
              x.ci = x.scuo.ci, y.ci = y.scuo.ci,
              xlim = xlim, ylim = ylim,
              xlab = "omega with phi", ylab = "omega without phi",
              main = "roc_ad_scuo_nonscale", workflow.name = workflow.name)
dev.off()

