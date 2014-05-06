start.time <- proc.time()

# Initial
suppressMessages(library(cubfits, quietly = TRUE))
set.seed(1234)

# Convert data.
reu13.list <- convert.reu13.df.to.list(ex.test$reu13.df)
y.list <- convert.y.to.list(ex.test$y)
n.list <- convert.n.to.list(ex.test$n)

# Get phi.Init.pred
init.function(model = "roc")
phi.Obs <- ex.train$phi.Obs / mean(ex.train$phi.Obs)
phi.Obs.test <- ex.test$phi.Obs / mean(ex.test$phi.Obs)
fitlist <- fitMultinom(ex.train$reu13.df, phi.Obs, ex.train$y, ex.train$n)
phi.Init.pred <- estimatePhi(fitlist, reu13.list, y.list, n.list,
                             E.Phi = median(phi.Obs.test),
                             lower.optim = min(phi.Obs.test) * 0.9,
                             upper.optim = max(phi.Obs.test) * 1.1)
phi.Init.pred <- phi.Init.pred / mean(phi.Init.pred)

# Run
.CF.AC$renew.iter <- 3
ret.time <- system.time({
  ret <- cubpred(ex.train$reu13.df, phi.Obs, ex.train$y, ex.train$n,
                 ex.test$reu13.df, ex.test$y, ex.test$n,
                 nIter = 10, burnin = 10,
                 phi.Init.pred = phi.Init.pred,
                 verbose = TRUE, report = 5,
                 model = "roc", adaptive = "simple")
})
print(ret.time)

# Report
x <- rowMeans(do.call("cbind", ret$phi.Mat.pred)[, 11:20])
y <- ex.test$phi.Obs
x <- log10(x / mean(x))
y <- log10(y / mean(y))
print(mean(x))
print(summary(lm(y ~ x))$r.squared)
# warning: iterations terminated because half-step sizes are very small

print(proc.time() - start.time)
