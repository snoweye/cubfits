\name{Posterior Results of Yassour2009}
\alias{yassour.PM.fits}
\alias{yassour.PM.appr}
\alias{yassour.info}
\docType{data}
\title{Posterior Results of Yassour 2009 Yeast Experiment Dataset}
\description{
  Output summarized from MCMC posterior results
  analyzing Yassour 2009 data.
}
\usage{
  yassour.PM.fits
  yassour.PM.appr
  yassour.info
}
\format{
  These are \code{list}'s containing several posterior means:
  \code{E.Phi} for expected expression, \code{b.InitList.roc} for parameters,
  \code{AA.prob} for proportion of amino acids, \code{sigmaW} for
  standard error of measure errors, and \code{gene.length} for
  gene length.
}
\details{
  \code{yassour.PM.fits} and \code{yassour.PM.appr} are the MCMC output
  of with/without observed expression, respectively.
  Both contain posterior means of expected expressions and coefficient
  parameters: \code{E.Phi} and \code{b.InitList.roc} are
  scaled results such that each MCMC iteration has mean 1 at \code{E.Phi}.

  \code{yassour.info} contains sequences information (Yeast):
  \code{AA.prob} and \code{gene.length} are summarized
  from corresponding genes in the analysis.

  Note that some of genes may not have good quality of expression or sequence
  information, so those genes are dropped from \code{\link{yassour}} dataset.
}
\references{
  \url{https://github.com/snoweye/cubfits/}
}
\seealso{
  \code{\link{yassour}}
}
\examples{
\dontrun{
str(yassour.PM.fits)
str(yassour.PM.appr)
str(yassour.PM.info)
}
}
\keyword{datasets}
