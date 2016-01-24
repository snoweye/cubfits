\name{Cedric Convergence Utilities}
\alias{cubmultichain}
\alias{cubsinglechain}
\alias{isConverged}
\title{Cedric Convergence Utilities}
\description{
  This utility function provides convergence related functions by Cedric.
}
\usage{
  cubmultichain(cubmethod, reset.qr, seeds=NULL,
                teston=c("phi", "sphi"), swap=0, swapAt=0.05, monitor=NULL, 
                min=0, max=160000, nchains=2, conv.thin=10,
                eps=0.1, ncores=2, ...)
                          
  cubsinglechain(cubmethod, frac1=0.1, frac2=0.5, reset.qr,
                 seed=NULL, teston=c("phi", "sphi"), monitor=NULL, 
                 min=0, max=160000, conv.thin=10, eps=1, ...)
}
\arguments{
  \item{cubmethod}{String to choose method. Options are "cubfits", "cubappr", "cubpred"}
  \item{reset.qr}{recalculate QR decomposition matrix of covariance matrix until reset.qr samples are reached}
  \item{swap}{proportion of b matrix parameters to be swaped between convergence checks}
  \item{swapAt}{difference (L1-norm) between two consequtive convergence test leading to a swap in the b matrix}
  \item{seeds}{Vector of seed for random number generation}
  \item{seed}{Seed for random number generation}
  \item{teston}{Select data to test convergence on}
  \item{monitor}{A function to monitor the progress of the MCMC. The fucntions expects the result object and for cubmultichain an index i. 
  (cubmultichain call: monitor(x,i), cubsinglechain call: monitor(x))}
  \item{min}{Minimum samples to be obtained. eps is ignored until number of samples reaches min}
  \item{max}{Maximum samples to be obtained. eps is ignored after max samples is obtained}
  \item{eps}{Convergence criterium}
  \item{conv.thin}{thinning of samples before performing convergence test}
  \item{nchains}{number of chains to run in parallel}
  \item{ncores}{number of cores to use for parallel execution of chains}
  \item{frac1}{fraction of samples at the beginning of set for Geweke test}
  \item{frac2}{fraction of samples at the end of set for Geweke test}
  \item{...}{named arguments for functions "cubfits", "cubappr" or "cubpred"}
}
\details{
  under development
}
\value{
  under development
}
\references{
  \url{https://github.com/clandere/cubfits/}
}
\author{
  Cedric Landerer \email{cedric.landerer@gmail.com}.
}
\seealso{
 cubfits, cubappr, cubpred
}
\examples{
\dontrun{
suppressMessages(library(cubfits, quietly = TRUE))

}
}
\keyword{cedric}
