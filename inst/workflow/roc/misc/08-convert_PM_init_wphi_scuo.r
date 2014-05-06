### This script extracts MCMC results and converts them to parameters.

rm(list = ls())

suppressMessages(suppressMessages(library(cubfits, quietly = TRUE)))

workflow.in <- c("pm", "sdlog_1.5", "sdlog_1.0", "sdlog_2.0")
workflow.out <- "./param_collected"

b.mcmc.all <- NULL
phi.mcmc.all <- NULL
sigmaW.all <- NULL
id.thin <- seq(1, 5000, by = 10)
for(i.workflow in workflow.in){
  if(i.workflow == workflow.in[1]){
    # Load sequence information.
    fn.in <- paste("./", i.workflow, "/all.out/data/pre_process.rda", sep = "")
    load(fn.in)

    # Obtain sequence information.
    AA.prob <- do.call("c", lapply(n, sum))
    AA.prob <- AA.prob / sum(AA.prob)
    names(AA.prob) <- names(n)

    # Dump sequence length.
    gene.length <- lapply(n.list, function(i.n){ do.call("sum", i.n) })
    gene.length <- do.call("c", gene.length)
    names(gene.length) <- names(n.list)
    fn.out <- paste("./", workflow.out, "/small_length.rda", sep = "")
    save(gene.length, file = fn.out)

    # Load workflow subset results.
    fn.in <- paste("./", i.workflow, "/all.out/subset/roc_ad_wphi_pm.rda",
                   sep = "")
    load(fn.in)

    b.mcmc.all <- cbind(b.mcmc.all, b.mcmc[, id.thin])
    phi.mcmc.all <- cbind(phi.mcmc.all, phi.mcmc[, id.thin])
    sigmaW.all <- c(sigmaW.all, p.mcmc[1, id.thin])
  } else{
    # Load workflow subset results.
    fn.in <- paste("./", i.workflow, "/all.out/subset/roc_ad_wphi_scuo.rda",
                   sep = "")
    load(fn.in)

    b.mcmc.all <- cbind(b.mcmc.all, b.mcmc[, id.thin])
    phi.mcmc.all <- cbind(phi.mcmc.all, phi.mcmc[, id.thin])
    sigmaW.all <- c(sigmaW.all, p.mcmc[1, id.thin])
  }
}

# Summarize collected thinning mcmc results.
tmp <- colMeans(phi.mcmc.all)

all.names <- rownames(b.mcmc.all)
# id.slop <- grep("(Intercept)", all.names, invert = TRUE)
id.slop <- grep("Delta.t", all.names)
b.mcmc.all[id.slop,] <- t(t(b.mcmc.all[id.slop,]) * tmp)
b.PM <- rowMeans(b.mcmc.all)

phi.PM <- colMeans(t(phi.mcmc.all) / tmp)
phi.Obs <- phi.PM

meanlog <- mean(log(phi.PM))
sdlog <- sd(log(phi.PM))
sigmaW <- mean(sigmaW.all)
print(c(meanlog, sdlog, sigmaW))

fn.out <- paste("./", workflow.out, "/small_train.rda", sep = "")
save(AA.prob, phi.Obs, meanlog, sdlog, sigmaW, file = fn.out)

# Write genome.phi.tsv
fn.out <- paste("./", workflow.out, "/genome.phi.tsv", sep = "")
da <- data.frame(ORF = names(phi.Obs), phi = phi.Obs)
write.table(da, file = fn.out, quote = FALSE, sep = "\t", row.names = FALSE)

# Dump parameters.
bInitList.roc <- convert.bVec.to.b(b.PM, names(n))
bInit.roc <- lapply(bInitList.roc, function(x){ x$coefficients })
fn.out <- paste("./", workflow.out, "/small_bInit.rda", sep = "")
save(bInitList.roc, bInit.roc, file = fn.out)

