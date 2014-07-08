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
  \item{reu13.df.obs}{under development}
  \item{ret.b.Mat}{under development}
  \item{phi.bin}{under development}
  \item{estim.phi}{under development}
  \item{n.use.iter}{under development}
  \item{rescale}{under development}
  \item{main}{under development}
  \item{xlab}{under development}
  \item{ylab}{under development}
  \item{model.label}{under development}
  \item{model.lty}{under development}
  \item{names.aa}{under development}
  \item{param}{under development}
  \item{phi.Mat}{under development}
  \item{ret.p.Mat}{under development}
  \item{lty}{under development}
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
\keyword{cedric}
