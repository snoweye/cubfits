### This script plots correlation of b (log(mu), omega).

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

  b.ci.org[[i.case]] <- b.ci.PM
  b.mean.org[[i.case]] <- b.PM
  label.org[[i.case]] <- b.label

  # Subset of mcmc output with scaling.
  fn.in <- paste(prefix$subset, case.names[i.case], "_PM_scaling.rda", sep = "")
  load(fn.in)

  b.ci[[i.case]] <- b.ci.PM
  b.mean[[i.case]] <- b.PM
  label[[i.case]] <- b.label
}


# Plot logmu.
all.names <- names(b.PM)
id.intercept <- grep("log.mu", all.names)

x.pm <- b.mean[[3]][id.intercept]
y.pm <- b.mean[[1]][id.intercept]
x.pm.label <- label[[3]]
x.scuo <- b.mean[[4]][id.intercept]
y.scuo <- b.mean[[2]][id.intercept]
x.scuo.label <- label[[4]]
xlim <- my.range(c(x.pm, x.scuo))
ylim <- my.range(c(y.pm, y.scuo))
x.pm.ci <- b.ci[[3]][id.intercept,]
y.pm.ci <- b.ci[[1]][id.intercept,]
x.scuo.ci <- b.ci[[4]][id.intercept,]
y.scuo.ci <- b.ci[[2]][id.intercept,]

fn.out <- paste(prefix$plot.match, "corr_logmu_pm.pdf", sep = "")
pdf(fn.out, width = 5, height = 5)
  plot.b.corr(x.pm, y.pm, x.pm.label,
              x.ci = x.pm.ci, y.ci = y.pm.ci,
              xlim = xlim, ylim = ylim,
              xlab = "log(mu) with phi", ylab = "log(mu) without phi",
              main = "roc_ad_pm", workflow.name = workflow.name)
dev.off()

fn.out <- paste(prefix$plot.match, "corr_logmu_scuo.pdf", sep = "")
pdf(fn.out, width = 5, height = 5)
  plot.b.corr(x.scuo, y.scuo, x.scuo.label,
              x.ci = x.scuo.ci, y.ci = y.scuo.ci,
              xlim = xlim, ylim = ylim,
              xlab = "log(mu) with phi", ylab = "log(mu) without phi",
              main = "roc_ad_scuo", workflow.name = workflow.name)
dev.off()


# Plot omega.
id.slop <- grep("omega", all.names)

x.pm <- b.mean[[3]][id.slop]
y.pm <- b.mean[[1]][id.slop]
x.pm.label <- label[[3]]
x.scuo <- b.mean[[4]][id.slop]
y.scuo <- b.mean[[2]][id.slop]
x.scuo.label <- label[[4]]
xlim <- my.range(c(x.pm, x.scuo))
ylim <- my.range(c(y.pm, y.scuo))
x.pm.ci <- b.ci[[3]][id.slop,]
y.pm.ci <- b.ci[[1]][id.slop,]
x.scuo.ci <- b.ci[[4]][id.slop,]
y.scuo.ci <- b.ci[[2]][id.slop,]

fn.out <- paste(prefix$plot.match, "corr_omega_pm.pdf", sep = "")
pdf(fn.out, width = 5, height = 5)
  plot.b.corr(x.pm, y.pm, x.pm.label,
              x.ci = x.pm.ci, y.ci = y.pm.ci,
              xlim = xlim, ylim = ylim,
              xlab = "omega with phi", ylab = "omega without phi",
              main = "roc_ad_pm", workflow.name = workflow.name)
dev.off()

fn.out <- paste(prefix$plot.match, "corr_omega_scuo.pdf", sep = "")
pdf(fn.out, width = 5, height = 5)
  plot.b.corr(x.scuo, y.scuo, x.scuo.label,
              x.ci = x.scuo.ci, y.ci = y.scuo.ci,
              xlim = xlim, ylim = ylim,
              xlab = "omega with phi", ylab = "omega without phi",
              main = "roc_ad_scuo", workflow.name = workflow.name)
dev.off()


# Plot omega original (no scaling by x to mean = 1).
id.slop <- grep("omega", all.names)

x.pm <- b.mean.org[[3]][id.slop]
y.pm <- b.mean.org[[1]][id.slop]
x.pm.label <- label.org[[3]]
x.scuo <- b.mean.org[[4]][id.slop]
y.scuo <- b.mean.org[[2]][id.slop]
x.scuo.label <- label.org[[4]]
xlim <- my.range(c(x.pm, x.scuo))
ylim <- my.range(c(y.pm, y.scuo))
x.pm.ci <- b.ci.org[[3]][id.slop,]
y.pm.ci <- b.ci.org[[1]][id.slop,]
x.scuo.ci <- b.ci.org[[4]][id.slop,]
y.scuo.ci <- b.ci.org[[2]][id.slop,]

fn.out <- paste(prefix$plot.match, "corr_nonscale_omega_pm.pdf", sep = "")
pdf(fn.out, width = 5, height = 5)
  plot.b.corr(x.pm, y.pm, x.pm.label,
              x.ci = x.pm.ci, y.ci = y.pm.ci,
              xlim = xlim, ylim = ylim,
              xlab = "omega with phi", ylab = "omega without phi",
              main = "roc_ad_pm_nonscale", workflow.name = workflow.name)
dev.off()

fn.out <- paste(prefix$plot.match, "corr_nonscale_omega_scuo.pdf", sep = "")
pdf(fn.out, width = 5, height = 5)
  plot.b.corr(x.scuo, y.scuo, x.scuo.label,
              x.ci = x.scuo.ci, y.ci = y.scuo.ci,
              xlim = xlim, ylim = ylim,
              xlab = "omega with phi", ylab = "omega without phi",
              main = "roc_ad_scuo_nonscale", workflow.name = workflow.name)
dev.off()

