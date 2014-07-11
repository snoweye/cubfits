\name{Cedric Convergence Utilities}
\alias{cubmultichain}
\alias{cubsinglechain}
\alias{isConverged}
\title{Cedric Convergence Utilities}
\description{
  This utility function provides convergence related functions by Cedric.
}
\usage{
  cubmultichain(cubmethod, niter, reset.qr, seeds, teston = c("phi", "sphi"),
                min = 0, max = 160000, nchains = 2, thin = 10, eps = 0.05,
                ncores = 2, ...)
  cubsinglechain(cubmethod, niter, frac1 = 0.1, frac2 = 0.5, reset.qr,
                 seed = NULL, teston = c("phi", "sphi"),
                 min = 0, max = 160000, thin = 10, eps = 0.05, ...)
  isConverged(chains, niter, epsilon = 0.1, thin = 10, frac1 = 0.1,
              frac2 = 0.5, teston = c("phi", "sphi"),
              test = c("gelman", "geweke"))
}
\arguments{
  \item{cubmethod}{under development}
  \item{niter}{under development}
  \item{seeds}{under development}
  \item{teston}{under development}
  \item{min}{under development}
  \item{max}{under development}
  \item{nchains}{under development}
  \item{thin}{under development}
  \item{eps}{under development}
  \item{ncores}{under development}
  \item{...}{under development}
  \item{frac1}{under development}
  \item{frac2}{under development}
  \item{reset.qr}{under development}
  \item{seed}{under development}
  \item{chains}{under development}
  \item{epsilon}{under development}
  \item{test}{under development}
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
  under development
}
\examples{
\dontrun{
suppressMessages(library(cubfits, quietly = TRUE))

}
}
\keyword{cedic}
