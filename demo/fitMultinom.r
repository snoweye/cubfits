rm(list = ls())

library(cubfits, quiet = TRUE)

# fit Shah & Gilchrist (2011)
init.function(model = "roc")
fitlist <- fitMultinom(ex.train$reu13.df, ex.train$phi.Obs,
                       ex.train$y, ex.train$n)
ret.fit <- prop.model.roc(fitlist, phi.Obs.lim = range(ex.train$phi.Obs))
aa.list <- names(ex.train$reu13.df)

# plot.
par(mfrow = c(1, 3))
for(i.aa in 1:length(aa.list)){
  plotmodel(ret.model = ret.fit[[i.aa]], main = aa.list[i.aa])
}
