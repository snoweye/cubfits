\name{Selection on Codon Usage}
\alias{calc_scu_values}
\title{Function for Selection on Codon Usage (SCU)}
\description{
  Calculate the average translational selection per transcript
  include mSCU and SCU (if gene expression is provided) for each gene.
}
\usage{
  calc_scu_values(b, y.list, phi.Obs = NULL)
}
\arguments{
  \item{b}{an object of format \code{\link{b}}.}
  \item{y.list}{an object of format \code{\link{y.list}}.}
  \item{phi.Obs}{an object of format \code{\link{phi.Obs}}, for
                 SCU only.}
}
\details{
  This function computes SCU and mSCU for each gene. Typically, this method
  is completely based on estimated parameters of mutation and selection
  such as outputs of MCMC or \code{\link{fitMultinom}()}.
}
\value{
  A list with two named elements \code{SCU} and \code{mSCU} are returned.
}
\references{
  Wallace E.W.J., Airoldi E.M., and Drummond D.A.
    ``Estimating Selection on Synonymous Codon Usage from Noisy Experimental
      Data''
    Mol Biol Evol (2013) 30(6):1438--1453.
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}.
}
\seealso{
  \code{\link{calc_scuo_values}()},
  \code{\link{calc_cai_values}()}.
}
\examples{
\dontrun{
library(cubfits, quietly = TRUE)

b <- b.Init$roc
phi.Obs <- ex.train$phi.Obs
y <- ex.train$y
y.list <- convert.y.to.list(y)
mSCU <- calc_scu_values(b, y.list, phi.Obs)$mSCU
plot(mSCU, log10(phi.Obs), main = "Expression vs mSCU",
     xlab = "mSCU", ylab = "Expression (log10)")

### Compare with CAI with weights seqinr::cubtab$sc.
library(seqinr, quietly = TRUE)
w <- caitab$sc
names(w) <- codon.low2up(rownames(caitab))
CAI <- calc_cai_values(y, y.list, w = w)$CAI

plot(mSCU, CAI, main = "CAI vs mSCU",
     xlab = "mSCU", ylab = "CAI")
}
}
\keyword{tool}
