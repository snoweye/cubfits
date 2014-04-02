rm(list = ls())
library(cubfits, quiet = TRUE)

seq.data <- read.seq(get.expath("seq_200.fasta"))
phi.df <- read.phi.df(get.expath("phi_200.tsv"))
aa.list <- c("A", "C", "D")

# Read in from FASTA file.
seq.string <- convert.seq.data.to.string(seq.data)
reu13.df <- gen.reu13.df(seq.string, phi.df, aa.list)
reu13.list.new <- gen.reu13.list(seq.string, aa.list)
y <- gen.y(seq.string, aa.list)
n <- gen.n(seq.string, aa.list)
scuo <- gen.scuo(seq.string, aa.list)

# Convert to list format.
reu13.list <- convert.reu13.df.to.list(reu13.df)
y.list <- convert.y.to.list(y)
n.list <- convert.n.to.list(n)
