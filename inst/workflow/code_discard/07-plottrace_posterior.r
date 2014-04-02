rm(list = ls())

library(cubfits, quiet = TRUE)

# Load environment and set data.
source("00-set_env.r")
source(paste(prefix$code, "u0-get_case_main.r", sep = ""))
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)

### Initial
init.function(model = model)

# Trace each run.
for(i.case in case.names){
  # All mcmc outputs
  fn.in <- paste(prefix$output, i.case, "/output_mcmc.rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }
  load(fn.in)

  if("sigmaW" %in% rownames(ret$p.Mat[[1]])){
    posterior <- lapply(1:length(ret$phi.Mat),
                   function(i.iter){
                     xx <- ret$phi.Mat[[i.iter]]
                     bInit <- convert.bVec.to.b(ret$b.Mat[[i.iter]],
                                                names(reu13.df.obs),
                                                model = model)
                     bInit <- lapply(bInit, function(B) B$coefficients)
                     sigmaWsq <- ret$p.Mat[[i.iter]][1]^2
                     mu.Phi <- ret$p.Mat[[i.iter]][2]
                     sigma.Phi.sq <- ret$p.Mat[[i.iter]][3]^2
                     tmp <- .cubfitsEnv$my.logPosteriorAll(xx, phi.Obs, y, n,
                               bInit, sigmaWsq, mu.Phi = mu.Phi, sigma.Phi.sq = sigma.Phi.sq,
                               reu13.df = reu13.df.obs)
                     sum(tmp)
                  })
  } else{
    posterior <- lapply(1:length(ret$phi.Mat),
                   function(i.iter){
                     xx <- ret$phi.Mat[[i.iter]]
                     bInit <- convert.bVec.to.b(ret$b.Mat[[i.iter]],
                                                names(reu13.df.obs),
                                                model = model)
                     bInit <- lapply(bInit, function(B) B$coefficients)
                     mu.Phi <- ret$p.Mat[[i.iter]][1]
                     sigma.Phi.sq <- ret$p.Mat[[i.iter]][2]^2
                     tmp <- .cubfitsEnv$my.logPosteriorAllPred(xx, y, n,
                               bInit, mu.Phi = mu.Phi, sigma.Phi.sq = sigma.Phi.sq,
                               reu13.df = reu13.df.obs)
                     sum(tmp)
                  })
  }

  # Load posterior mean results.
  fn.in <- paste(prefix$subset, i.case, "_PM.rda", sep = "")
  load(fn.in)
  fn.in <- paste(prefix$subset, i.case, "_PM_scaling.rda", sep = "")
  load(fn.in)

  bInit <- convert.bVec.to.b(b.PM, names(reu13.df.obs), model = model)
  bInit <- lapply(bInit, function(B) B$coefficients)

  if("sigmaW" %in% rownames(p.PM)){
    sigmaWsq <- p.PM[1]^2
    mu.Phi <- p.PM[2]
    sigma.Phi.sq <- p.PM[3]^2
    tmp <- .cubfitsEnv$my.logPosteriorAll(phi.PM, phi.Obs, y, n,
              bInit, sigmaWsq, mu.Phi = mu.Phi, sigma.Phi.sq = sigma.Phi.sq,
              reu13.df = reu13.df.obs)
  } else{
    mu.Phi <- p.PM[1]
    sigma.Phi.sq <- p.PM[2]^2
    tmp <- .cubfitsEnv$my.logPosteriorAllPred(phi.PM, y, n,
              bInit, mu.Phi = mu.Phi, sigma.Phi.sq = sigma.Phi.sq,
              reu13.df = reu13.df.obs)
  }
  posterior.PM <- sum(tmp)

  x <- 1:length(ret$phi.Mat)
  xlim <- range(x)
  ylim <- range(range(posterior), posterior.PM)

  # Trace of posterior
  fn.out <- paste(prefix$plot.trace, "posterior_", i.case,
                  ".pdf", sep = "")
  pdf(fn.out, width = 6, height = 4)
    plot(NULL, NULL, xlim = xlim, ylim = ylim,
         main = paste("PM of prop. posterior = ",
                      sprintf("%.4f", posterior.PM), sep = ""),
         xlab = "Iterations", ylab = "Prop. Posterior")
    mtext(paste(workflow.name, ", ", get.case.main(i.case, model), sep = ""),
          line = 3)
    lines(x = x, y = posterior)
    abline(h = posterior.PM, col = 2)
  dev.off()

  # Dump posterior
  fn.out <- paste(prefix$subset, "trace_posterior_", i.case, ".rda", sep = "")
  save(posterior, posterior.PM, file = fn.out)
}

