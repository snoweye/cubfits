rm(list = ls())

source("00-set_env.r")

b.mcmc.ci <- list()
x.mcmc.ci <- list()

for(i.case in 1:length(case.names)){
  # subset of mcmc output
  fn.in <- paste(prefix$subset, case.names[i.case], ".rda", sep = "")
  if(! file.exists(fn.in)){
    next
  }
  load(fn.in)

  all.names <- rownames(b.mcmc)

  id <- grep("(Intercept)", all.names, invert = TRUE)
  scale.x <- mean(colMeans(x.mcmc))

  b.mcmc[id,] <- b.mcmc[id,] * scale.x
  x.mcmc <- x.mcmc / scale.x

  tmp <- apply(b.mcmc, 1, quantile, prob = ci.prob)
  b.mcmc.ci[[i.case]] <- abs(tmp[2,] - tmp[1,])

  x.mcmc <- log10(x.mcmc)
  tmp <- apply(x.mcmc, 1, quantile, prob = ci.prob)
  x.mcmc.ci[[i.case]] <- abs(tmp[2,] - tmp[1,])
}

# plot M
id <- grep("(Intercept)", all.names, invert = FALSE)
m.mcmc.ci <- lapply(b.mcmc.ci, function(x) x[id])
xlim <- range(do.call("c", lapply(m.mcmc.ci, range)))
for(i.case in 1:length(case.names)){
  fn.out <- paste(prefix$plot.diag, "ci_M_", case.names[i.case],
                  ".pdf", sep = "")
  pdf(fn.out, width = 5, height = 5)
    hist(m.mcmc.ci[[i.case]], xlim = xlim, nclass = 20,
         xlab = "M CI width",
         main = case.names[i.case])
    title(sub = workflow.name, cex.sub = 0.6)
  dev.off()
}

# plot S
id <- grep("(Intercept)", all.names, invert = TRUE)
s.mcmc.ci <- lapply(b.mcmc.ci, function(x) x[id])
xlim <- range(do.call("c", lapply(s.mcmc.ci, range)))
for(i.case in 1:length(case.names)){
  fn.out <- paste(prefix$plot.diag, "ci_S_", case.names[i.case],
                  ".pdf", sep = "")
  pdf(fn.out, width = 5, height = 5)
    hist(s.mcmc.ci[[i.case]], xlim = xlim, nclass = 20,
         xlab = "S CI width",
         main = case.names[i.case])
    title(sub = workflow.name, cex.sub = 0.6)
  dev.off()
}

# plot EX_g
xlim <- range(do.call("c", lapply(x.mcmc.ci, range)))
for(i.case in 1:length(case.names)){
  fn.out <- paste(prefix$plot.diag, "ci_EXg_", case.names[i.case],
                  ".pdf", sep = "")
  pdf(fn.out, width = 5, height = 5)
    hist(x.mcmc.ci[[i.case]], xlim = xlim, nclass = 40,
         xlab = "log10(X_g) CI width",
         main = case.names[i.case])
    title(sub = workflow.name, cex.sub = 0.6)
  dev.off()
}

