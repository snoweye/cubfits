\name{Cedric Plot Utilities}
\alias{plotCUB}
\alias{plotTraces}
\alias{plotExpectedPhiTrace}
\alias{plotPTraces}
\title{Cedric Plot Utilities}
\description{
  This utility function provides basic plots by Cedric.
}
\usage{
  plotCUB(reu13.df.obs, ret.b.Mat, phi.bin, estim.phi, n.use.iter=2000,
          rescale=F, main="CUB", model.label=c("True Model"), model.lty=1)

  plotPTraces(ret.p.Mat, main="Hyperparameter Traces", lty=1)

  plotExpectedPhiTrace(phi.Mat, 
                       main=expression(paste("Trace E[", phi, "]", sep="")),
                       xlab="Iteration",
                       ylab=expression(paste("E[", phi, "]", sep="")), lty=1)

  plotTraces(ret.b.Mat, names.aa, param = c("logmu", "deltat"),
             main="AA parameter trace")
}
\arguments{
  \item{reu13.df.obs}{NTD}
  \item{ret.b.Mat}{NTD}
  \item{phi.bin}{NTD}
  \item{estim.phi}{NTD}
  \item{n.use.iter}{NTD}
  \item{rescale}{NTD}
  \item{main}{NTD}
  \item{xlab}{NTD}
  \item{ylab}{NTD}
  \item{model.label}{NTD}
  \item{model.lty}{NTD}
  \item{names.aa}{NTD}
  \item{param}{NTD}
  \item{phi.Mat}{NTD}
  \item{ret.p.Mat}{NTD}
  \item{lty}{NTD}
}
\details{
  NTD
}
\value{
  NTD
}
\references{
  \url{https://github.com/snoweye/cubfits/}
}
\author{
  Cedric Landerer \email{cedric.landerer@gmail.com}.
}
\seealso{
  NTD
}
\examples{
\dontrun{
suppressMessages(library(cubfits, quietly = TRUE))

}
}
\keyword{cedric}
