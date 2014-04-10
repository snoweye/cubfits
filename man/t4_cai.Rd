\name{Codon Adaptation Index}
\alias{calc_cai_values}
\title{Function for Codon Adaptation Index (CAI)}
\description{
  Calculate the Codon Adaptation Index (CAI) for each gene.
  Used as a substitute for expression in cases of without expression
  measurements.
}
\usage{
  calc_cai_values(y, y.list, w = NULL)
}
\arguments{
  \item{y}{an object of format \code{\link{y}}.}
  \item{y.list}{an object of format \code{\link{y.list}}.}
  \item{w}{a specified relative frequency of synonymous codons.}
}
\details{
  This function computes CAI for each gene. Typically, this method
  is completely based on entropy and information theory to estimate
  expression values of sequences according to their codon information.

  If the input \code{w} is \code{NULL}, then empirical values are computed.
}
\value{
  A list with two named elements \code{CAI} and \code{w} are returned
  where \code{CAI} are CAI of input sequences (\code{y} and \code{y.list})
  and \code{w} are the relative frequencey used to computed those CAI's.
}
\references{
  Sharp P.M. and Li W.-H.
    ``The codon Adaptation Index -- a measure of directional synonymous codon
      usage bias, and its potential applications''
    Nucleic Acids Res. 15 (3): 1281-1295, 1987.
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}.
}
\seealso{
  \code{\link{calc_scuo_values}()},
  \code{\link{calc_scu_values}()}.
}
\examples{
\dontrun{
rm(list = ls())
library(cubfits, quietly = TRUE)

y <- ex.train$y
y.list <- convert.y.to.list(y)
CAI <- calc_cai_values(y, y.list)$CAI
plot(CAI, log10(ex.train$phi.Obs), main = "Expression vs CAI",
     xlab = "CAI", ylab = "Expression (log10)")

### Verify with the seqinr example.
library(seqinr, quietly = TRUE)
inputdatfile <- system.file("sequences/input.dat", package = "seqinr")
input <- read.fasta(file = inputdatfile, forceDNAtolower = FALSE)
names(input)[65] <- paste(names(input)[65], ".1", sep = "") # name duplicated.
input <- input[order(names(input))]

### Convert to cubfits format.
seq.string <- convert.seq.data.to.string(input)
new.y <- gen.y(seq.string)
new.y.list <- convert.y.to.list(new.y)
ret <- calc_cai_values(new.y, new.y.list)

### Rebuild w.
w <- rep(1, 64)
names(w) <- codon.low2up(rownames(caitab))
for(i in 1:64){
  id <- which(names(ret$w) == names(w)[i])
  if(length(id) == 1){
    w[i] <- ret$w[id]
  }
}
CAI.res <- sapply(input, seqinr::cai, w = w)

### Plot.
plot(CAI.res, ret$CAI,
     main = "Comparison of seqinR and cubfits results",
     xlab = "CAI from seqinR", ylab = "CAI from cubfits", las = 1)
abline(c(0, 1))
}
}
\keyword{tool}
