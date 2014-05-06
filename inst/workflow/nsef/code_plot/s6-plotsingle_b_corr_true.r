### This is for simulation only, ploting correlation aganst true values.

rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

source("00-set_env.r")
source(paste(prefix$code.plot, "u0-get_case_main.r", sep = ""))
source(paste(prefix$code.plot, "u2-plot_b_corr.r", sep = ""))

# Load true Phi.
fn.in <- paste(prefix$data, "simu_true_", model, ".rda", sep = "")
if(file.exists(fn.in)){
  load(fn.in)
} else{
  stop(paste(fn.in, " is not found.", sep = ""))
}

bInit <- convert.b.to.bVec(Eb)

# Load data.
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)

# Get AA and synonymous codons.
aa.names <- names(reu13.df.obs)
coef.names <- cubfits:::get.my.coefnames(model)
b.names <- NULL
for(i.aa in aa.names){
  tmp <- sort(unique(reu13.df.obs[[i.aa]]$Codon))
  tmp <- tmp[-length(tmp)]
  b.names <- c(b.names, rep(coef.names, each = length(tmp)))
}

# Get true values.
all.names <- b.names
id.slop <- grep("omega", all.names)
scale.EPhi <- mean(EPhi)
bInit[id.slop] <- bInit[id.slop] * scale.EPhi

for(i.case in case.names){
  title <- paste(workflow.name, ", ", get.case.main(i.case, model), sep = "")

  # Subset of mcmc output.
  fn.in <- paste(prefix$subset, i.case, "_PM_scaling.rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }
  load(fn.in)

  # Plot log(mu)
  id.intercept <- grep("log.mu", all.names)
  x <- bInit[id.intercept]
  y <- b.PM[id.intercept]
  xlim <- my.range(x)
  ylim <- my.range(y)
  y.ci <- b.ci.PM[id.intercept,]

  fn.out <- paste(prefix$plot.single, "corr_true_logmu_", i.case,
                  ".pdf", sep = "")
  pdf(fn.out, width = 5, height = 5)
    plot.b.corr(x, y, b.label,
                y.ci = y.ci,
                xlim = xlim, ylim = ylim,
                xlab = "True", ylab = "Estimated",
                main = "log(mu)", workflow.name = title)
  dev.off()

  # Plot omega
  id.slop <- grep("omega", all.names)
  x <- bInit[id.slop]
  y <- b.PM[id.slop]
  xlim <- my.range(x)
  ylim <- my.range(y)
  y.ci <- b.ci.PM[id.slop,]

  fn.out <- paste(prefix$plot.single, "corr_true_omega_", i.case,
                  ".pdf", sep = "")
  pdf(fn.out, width = 5, height = 5)
    plot.b.corr(x, y, b.label,
                y.ci = y.ci,
                xlim = xlim, ylim = ylim,
                xlab = "True", ylab = "Estimated",
                main = "omega", workflow.name = title)
  dev.off()
}
