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
  plotPTraces(pMat, ...)
  
  plotExpectedPhiTrace(phiMat, ...)
  
  plotCUB(reu13.df.obs, bMat = NULL, bVec = NULL, phi.bin,
          n.use.samples = 2000, main = "CUB", model.label = c("True Model"),
          model.lty = 1, weightedCenters = TRUE)

  plotTraces(bMat, names.aa, param = c("logmu", "deltaeta", "deltat"),
             main = "AA parameter trace")
}
\arguments{
  \item{reu13.df.obs}{under development}
  \item{bVec}{a parameter vector}
  \item{phi.bin}{phi values to bin for comparison}
  \item{n.use.samples}{under development}
  \item{main}{Main name for plotTraces}
  \item{model.label}{Name of model}
  \item{model.lty}{line type for model}
  \item{weightedCenters}{if centers are weighted.}

  \item{names.aa}{List of amino acids used for estimation}
  \item{param}{select to plot parameter trace for either log(mu) values or
                delta t}

  \item{phiMat}{phi matrix from the output of "cubmultichain",
                "cubsinglechain", "cubfits", "cubappr", or "cubpred"}
  \item{bMat}{b matrix from the output of "cubmultichain", "cubsinglechain",
              "cubfits", "cubappr", or "cubpred"}
  \item{pMat}{p matrix from the output of "cubmultichain", "cubsinglechain",
              "cubfits", "cubappr", or "cubpred"}
  \item{...}{other ploting options}
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
  plot
}
\examples{
\dontrun{
suppressMessages(library(cubfits, quietly = TRUE))

}
}
\keyword{cedric}
