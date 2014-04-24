### This script extracts MCMC results and converts them to parameters.

rm(list = ls())

suppressMessages(library(cubfits))

workflow.in <- "org-02-fitsappr_yassour"
workflow.out <- "41-simu_yassour"

# Load sequence information.
fn.in <- paste("./", workflow.in, "/all.out/data/pre_process.rda", sep = "")
load(fn.in)

# Dump sequence length.
gene.length <- lapply(n.list, function(i.n){ do.call("sum", i.n) })
gene.length <- do.call("c", gene.length)
names(gene.length) <- names(n.list)
fn.out <- paste("./", workflow.out, "/param/small_length.rda", sep = "")
save(gene.length, file = fn.out)

# Load workflow subset results.
fn.in <- paste("./", workflow.in, "/all.out/subset/roc_ad_fits_pm_PM.rda",
               sep = "")
load(fn.in)
fn.in <- paste("./", workflow.in, "/all.out/subset/roc_ad_fits_pm_PM_scaling.rda",
               sep = "")
load(fn.in)

# Dump sequence information.
AA.prob <- do.call("c", lapply(n, sum))
AA.prob <- AA.prob / sum(AA.prob)
names(AA.prob) <- names(n)
phi.Obs <- phi.PM
meanlog <- mean(log(phi.PM))
sdlog <- sd(log(phi.PM))
sigmaW <- p.PM[1]
fn.out <- paste("./", workflow.out, "/param/small_train.rda", sep = "")
save(AA.prob, phi.Obs, meanlog, sdlog, sigmaW, file = fn.out)

# Dump parameters.
bInitList.roc <- convert.bVec.to.b(b.PM, names(n))
bInit.roc <- lapply(bInitList.roc, function(x){ x$coefficients })
fn.out <- paste("./", workflow.out, "/param/small_bInit.rda", sep = "")
save(bInitList.roc, bInit.roc, file = fn.out)
