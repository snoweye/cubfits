rm(list = ls())
start.time <- proc.time()

# Initial
library(cubfits, quiet = TRUE)
set.seed(1234)

# Get xInit.appr
xInit.appr <- ex.test$phi.Obs / mean(ex.test$phi.Obs)

# Run
init.function(model = "roc")
.CF.AC$renew.iter <- 3
ret.time <- system.time({
  ret <- cubappr(ex.test$reu13.df, xInit.appr, ex.test$y, ex.test$n,
                 nIter = 10, burnin = 10,
                 verbose = TRUE, report = 5,
                 model = "roc", adaptive = "simple")
})
print(ret.time)

# Report
x <- rowMeans(do.call("cbind", ret$phi.Mat)[, 11:20])
y <- ex.test$phi.Obs
x <- log10(x / mean(x))
y <- log10(y / mean(y))
cat("mean of x  : ", mean(x), "\n", sep = "")
cat("R^2 (log10): ", summary(lm(y ~ x))$r.squared, "\n", sep = "")

print(proc.time() - start.time)
