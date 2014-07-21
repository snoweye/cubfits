\name{Cedric Convergence Utilities}
\alias{cubmultichain}
\alias{cubsinglechain}
\alias{isConverged}
\title{Cedric Convergence Utilities}
\description{
  This utility function provides convergence related functions by Cedric.
}
\usage{
  cubmultichain(cubmethod, nsamples, reset.qr, seeds, teston = c("phi", "sphi"),
                min = 0, max = 160000, nchains = 2, thin = 10, eps = 0.05,
                ncores = 2, ...)
  cubsinglechain(cubmethod, nsamples, frac1=0.1, frac2=0.5, reset.qr, 
                seed=NULL, teston=c("phi", "sphi"), monitor=NULL, 
                growthfactor=0, min=0, max=160000, conv.thin=10, eps=0.05, ...)
  isConverged(chains, nsamples, eps = 0.1, thin = 10, frac1 = 0.1,
              frac2 = 0.5, teston = c("phi", "sphi"),
              test = c("gelman", "geweke"))
}
\arguments{
  \item{cubmethod}{String to choose method. Options are "cubfits", "cubappr", "cubpred"}
  \item{nsamples}{Number of samples to obtain before checking for convergence}
  \item{seeds}{Set seed for random number generation}
  \item{seed}{seed for random number generation}
  \item{teston}{select data to test convergence on}
  \item{min}{Minimum samples to be obtained. eps is ignored until number of samples reaches min}
  \item{max}{Maximum samples to be obtained. eps is ignored after max samples is obtained}
  \item{nchains}{number of chains to run in parallel}
  \item{thin}{thining of samples before convergence test}
  \item{eps}{Convergence criterium}
  \item{ncores}{number of cores to use for parallel execution of chains}
  \item{frac1}{fraction of samples at the beginning of set for Geweke test}
  \item{frac2}{fraction of samples at the end of set for Geweke test}
  \item{reset.qr}{recalculate QR decomposition matrix of covariance matrix until reset.qr samples are reached}
  \item{chains}{chains to test convergence on}
  \item{test}{select convergence test. Gelman for multichain convergence and Geweke for Singlechain convergence}
  \item{...}{named arguments for functions "cubfits", "cubappr" or "cubpred"}
}
\details{
  under development
}
\value{
  under development
}
\references{
  \url{https://github.com/snoweye/cubfits/}
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
