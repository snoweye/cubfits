### This script collect the poster means of MCMC runs.

rm(list = ls())

source("00-set_env.r")

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

# Get all cases.
for(i.case in case.names){
  # Check first.
  fn.in <- paste(prefix$output, i.case, "/output_mcmc.rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }

  # Load MCMC output.
  load(fn.in)

  # Obtain the last 5000 iterations.
  b.mcmc <- do.call("cbind", ret$b.Mat[range$subset])
  p.mcmc <- do.call("cbind", ret$p.Mat[range$subset])
  phi.mcmc <- do.call("cbind", ret$phi.Mat[range$subset])

  rownames(b.mcmc) <- names(ret$b.Mat[[1]])
  rownames(p.mcmc) <- names(ret$p.Mat[[1]])
  rownames(phi.mcmc) <- names(ret$phi.Mat[[1]])

  # Dump posterior distributions.
  fn.out <- paste(prefix$subset, i.case, ".rda", sep = "")
  save(b.mcmc, p.mcmc, phi.mcmc, file = fn.out)

  # Obtain posterior means.
  b.PM <- rowMeans(b.mcmc)
  b.ci.PM <- apply(b.mcmc, 1, quantile, prob = ci.prob)
  p.PM <- rowMeans(p.mcmc)
  phi.PM <- rowMeans(phi.mcmc)

  # Negative selection.
  ret <- get.negsel(b.PM, b.ci.PM, id.slop, aa.list, label)
  b.negsel.PM <- ret$b.negsel.PM
  b.negsel.ci.PM <- ret$b.negsel.ci.PM
  label.negsel <- ret$label.negsel

  # Dump summarized results.
  fn.out <- paste(prefix$subset, i.case, "_PM.rda", sep = "")
  save(b.PM, b.ci.PM, p.PM, phi.PM, b.negsel.PM, b.negsel.ci.PM,
       label, label.negsel,
       file = fn.out)

  ### Thinning.
  # id.save <- (range$subset %% range$thinning) == 1
  # b.mcmc <- b.mcmc[, id.save]
  # p.mcmc <- p.mcmc[, id.save]
  # phi.mcmc <- phi.mcmc[, id.save]

  ### Dump posterior distributions.
  # fn.out <- paste(prefix$subset, i.case, "_thinning.rda", sep = "")
  # save(b.mcmc, p.mcmc, phi.mcmc, file = fn.out)

  ### Obtain posterior means.
  # b.PM <- rowMeans(b.mcmc)
  # p.PM <- rowMeans(p.mcmc)
  # phi.PM <- rowMeans(phi.mcmc)

  ### Dump posterior means.
  # fn.out <- paste(prefix$subset, i.case, "_PM_thinning.rda", sep = "")
  # save(b.PM, p.PM, phi.PM, file = fn.out)

  # Scaling.
  scale.EPhi <- colMeans(phi.mcmc)
  phi.mcmc <- t(t(phi.mcmc) / scale.EPhi)
  all.names <- rownames(b.mcmc)
  id.slop <- grep("(Intercept)", all.names, invert = TRUE)
  b.mcmc[id.slop,] <- t(t(b.mcmc[id.slop,]) * scale.EPhi)
  b.PM <- rowMeans(b.mcmc)
  b.ci.PM <- apply(b.mcmc, 1, quantile, prob = ci.prob)
  phi.PM <- rowMeans(phi.mcmc)

  # Negative selection.
  ret <- get.negsel(b.PM, b.ci.PM, id.slop, aa.list, label)
  b.negsel.PM <- ret$b.negsel.PM
  b.negsel.ci.PM <- ret$b.negsel.ci.PM
  label.negsel <- ret$label.negsel

  # Dump summarized results.
  fn.out <- paste(prefix$subset, i.case, "_PM_scaling.rda", sep = "")
  save(b.PM, b.ci.PM, phi.PM, b.negsel.PM, b.negsel.ci.PM, label.negsel,
       label, label.negsel,
       file = fn.out)
}

