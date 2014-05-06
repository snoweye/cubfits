suppressMessages(library(cubfits, quietly = TRUE))
set.seed(1234)

# Normalize phi.Obs
phi.Obs <- ex.train$phi.Obs / mean(ex.train$phi.Obs)

# Generate sequences.
da.nsef <- simu.orf(length(phi.Obs), bInit$nsef,
                    phi.Obs = phi.Obs, model = "nsef")
names(da.nsef) <- names(phi.Obs)
write.seq(da.nsef, "toy_nsef.fasta")

# Read seqeuences back.
seq.nsef <- read.seq("toy_nsef.fasta")
seqstring.nsef <- convert.seq.data.to.string(seq.nsef)
phi <- data.frame(ORF = names(phi.Obs), phi.value = phi.Obs)

# Generate data structures from sequences.
aa.names <- names(bInit$nsef)
reu13.df <- gen.reu13.df(seqstring.nsef, phi, aa.names = aa.names)
n <- gen.n(seqstring.nsef, aa.names = aa.names)
y <- gen.y(seqstring.nsef, aa.names = aa.names)

# Run codon fits.
.CF.AC$renew.iter <- 3
ret.time <- system.time({
  ret <- cubfits(reu13.df, phi.Obs, y, n,
                 nIter = 10, burnin = 10,
                 phi.DrawScale = 0.01,
                 verbose = TRUE, report = 5,
                 model = "nsef", adaptive = "simple")
})

x <- rowMeans(do.call("cbind", ret$phi.Mat)[, 11:20])
y <- phi.Obs
plotprxy(x, y)

x <- log10(x / mean(x))
y <- log10(y / mean(y))
print(mean(x))
print(summary(lm(y ~ x))$r.squared)

