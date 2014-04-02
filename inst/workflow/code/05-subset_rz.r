### This script collect the poster means of MCMC runs.

rm(list = ls())

library(cubfits, quiet = TRUE)
source("00-set_env.r")
source(paste(prefix$code, "u1-get_negsel.r", sep = ""))

# Load data.
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)

# Get AA and synonymous codons.
aa.list <- names(reu13.df.obs)
b.label <- NULL
for(i.aa in aa.list){
  tmp <- sort(unique(reu13.df.obs[[i.aa]]$Codon))
  tmp <- tmp[-length(tmp)]
  b.label <- c(b.label, paste(i.aa, tmp, sep = "."))
}

# Get all cases.
for(i.case in case.names){
  # Check first.
  fn.in <- paste(prefix$subset, i.case, ".rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }

  # Load MCMC output.
  load(fn.in)
  all.names <- rownames(b.mcmc)
  id.slop <- grep("(Intercept)", all.names, invert = TRUE)
  b.PM <- rowMeans(b.mcmc)
  b.ci.PM <- t(apply(b.mcmc, 1, quantile, prob = ci.prob))

  # Scaling.
  t.phi.mcmc <- t(phi.mcmc)
  x.phi.mcmc <- as.vector(t.phi.mcmc / colMeans(phi.mcmc))
  n.1 <- length(x.phi.mcmc) - 1
  for(i.b in id.slop){
    y.phi.mcmc <- as.vector(t.phi.mcmc * b.mcmc[i.b,])
    ### Slower.
    # m.1 <- lm(y.phi.mcmc ~ -1 + x.phi.mcmc)
    # b.PM[i.b] <- as.vector(m.1$coefficients)
    # sm.1 <- summary(m.1)
    # b.ci.PM[i.b,] <- sm.1$coefficients[1, 1] +
    #                  c(-1, 1) * sm.1$coefficients[1, 2] * 1.96
    ### Tests.
    # n.1 <- 5000 - 1
    # x.phi.mcmc <- rnorm(n.1 + 1)
    # y.phi.mcmc <- rnorm(n.1 + 1)
    # summary(lm(y.phi.mcmc ~ -1 + x.phi.mcmc))$coefficients
    x.tx <- sum(x.phi.mcmc * x.phi.mcmc)
    beta1.hat <- sum(x.phi.mcmc * y.phi.mcmc) / x.tx
    se1.hat <- sqrt(sum((y.phi.mcmc - x.phi.mcmc * beta1.hat)^2) /
                    (n.1 * x.tx))
    b.PM[i.b] <- beta1.hat
    b.ci.PM[i.b,] <- beta1.hat + qt(ci.prob, n.1) * se1.hat * sqrt(n.1)
  }

  # Negative selection.
  ret <- get.negsel(b.PM, id.slop, aa.list, b.label, b.ci.PM = b.ci.PM)
  b.negsel.PM <- ret$b.negsel.PM
  b.negsel.ci.PM <- ret$b.negsel.ci.PM
  b.negsel.label <- ret$b.negsel.label

  # Dump summarized results.
  fn.out <- paste(prefix$subset, i.case, "_PM_scaling_rz.rda", sep = "")
  save(b.PM, b.ci.PM, b.negsel.PM, b.negsel.ci.PM,
       b.label, b.negsel.label,
       file = fn.out)
} # End of all.jobs().

