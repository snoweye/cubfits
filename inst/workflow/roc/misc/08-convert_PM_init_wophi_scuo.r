### This script extracts MCMC results and converts them to parameters.

rm(list = ls())

suppressMessages(library(cubfits))

workflow.in <- c("sdlog_1.5", "sdlog_1.0", "sdlog_2.0")
workflow.out <- "./output_collected"

b.mcmc.all <- NULL
phi.mcmc.all <- NULL
p.mcmc.all <- NULL
id.thin <- seq(1, 5000, by = 10)
for(i.workflow in workflow.in){
  if(i.workflow == workflow.in[1]){
    # Load sequence information.
    fn.in <- paste("./", i.workflow, "/all.out/data/pre_process.rda", sep = "")
    load(fn.in)
  }

  # Load workflow subset results.
  fn.in <- paste("./", i.workflow, "/all.out/subset/roc_ad_wophi_scuo.rda",
                 sep = "")
  load(fn.in)

  b.mcmc.all <- cbind(b.mcmc.all, b.mcmc[, id.thin])
  phi.mcmc.all <- cbind(phi.mcmc.all, phi.mcmc[, id.thin])
  p.mcmc.all <- cbind(p.mcmc.all, p.mcmc[, id.thin])
}
rownames(b.mcmc.all) <- rownames(b.mcmc)
rownames(phi.mcmc.all) <- rownames(phi.mcmc)
rownames(p.mcmc.all) <- rownames(p.mcmc)

# Dump
fn.out <- paste("./", workflow.out, "/output_mcmc_all.rda", sep = "")
save(b.mcmc.all, phi.mcmc.all, p.mcmc.all, file = fn.out)


# Summarize collected thinning mcmc results.
tmp <- colMeans(phi.mcmc.all)

all.names <- rownames(b.mcmc.all)
# id.slop <- grep("(Intercept)", all.names, invert = TRUE)
id.slop <- grep("Delta.t", all.names)
b.mcmc.all[id.slop,] <- t(t(b.mcmc.all[id.slop,]) * tmp)
b.PM <- rowMeans(b.mcmc.all)

phi.PM <- colMeans(t(phi.mcmc.all) / tmp)
phi.Obs <- phi.PM

fn.out <- paste("./", workflow.out, "/small_train.rda", sep = "")
save(b.PM, phi.PM, phi.Obs, file = fn.out)

# Write genome.phi.tsv
fn.out <- paste("./", workflow.out, "/genome.phi.tsv", sep = "")
da <- data.frame(ORF = names(phi.Obs), phi = phi.Obs)
write.table(da, file = fn.out, quote = FALSE, sep = "\t", row.names = FALSE)

# Dump parameters.
bInitList.roc <- convert.bVec.to.b(b.PM, names(n))
bInit.roc <- lapply(bInitList.roc, function(x){ x$coefficients })
fn.out <- paste("./", workflow.out, "/small_bInit.rda", sep = "")
save(bInitList.roc, bInit.roc, file = fn.out)

