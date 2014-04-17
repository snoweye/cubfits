\name{CUB Model Fits}
\alias{cubfits}
\title{Codon Usage Bias Fits for Observed ORFs and Expression}
\description{
  This function provides codon usage bias fits with observed ORFs and
  expressions which possibly contains measurement errors.
}
\usage{
  cubfits(reu13.df.obs, phi.Obs, y, n,
          nIter = 1000, burnin = 100,
          pInit = NULL,
          bInit = NULL, init.b.Scale = 1, b.DrawScale = 1,
          phi.Init = NULL, init.Phi.Scale = 1, phi.DrawScale = 1,
          model = .CF.CT$model[1], model.Phi = .CF.CT$model.Phi[1],
          adaptive = .CF.CT$adaptive[1], scale.Phi = .CF.CT$scale.Phi[1],
          verbose = .CF.DP$verbose,
          iterThin = .CF.DP$iterThin, report = .CF.DP$report)
}
\arguments{
  \item{reu13.df.obs}{a \code{\link{reu13.df}} object, ORFs information.}
  \item{phi.Obs}{a \code{\link{phi.Obs}} object, expression with measurement errors.}
  \item{y}{a \code{\link{y}} object, codon counts.}
  \item{n}{a \code{\link{n}} object, total codon counts.}
  \item{nIter}{number of iterations after burn-in iterations.}
  \item{burnin}{number of burn-in iterations.}
  \item{pInit}{initial values for hyper-parameters.}
  \item{bInit}{initial values for parameters \code{\link{b}}.}
  \item{init.b.Scale}{for initial \code{\link{b}} if \code{bInit = NULL}.}
  \item{b.DrawScale}{scaling factor for adaptive MCMC with random walks
                     when drawing new \code{\link{b}}.}
  \item{phi.Init}{initial values for Phi.}
  \item{init.Phi.Scale}{for initial Phi if \code{phi.Init = NULL}.}
  \item{phi.DrawScale}{scaling factor for adaptive MCMC with random walks
                       when drawing new Phi.}
  \item{model}{model to be fitted, currently "roc" only.}
  \item{model.Phi}{prior model for Phi, currently "lognormal".}
  \item{adaptive}{adaptive method of MCMC for proposing new \code{\link{b}} and
                  Phi.}
  \item{scale.Phi}{method of checking pre-scaled Phi.}
  \item{verbose}{print iteration messages.}
  \item{iterThin}{thinning iterations.}
  \item{report}{number of iterations to report more information.}
}
\details{
  This function correctly and carefully implements a combining version of
  Shah and Gilchrist (2011) and Wallace et al. (2013).

  Total number of MCMC iterations is \code{burnin + nIter + 1}, but the
  outputs may be thinned to \code{(burnin + nIter) / iterThin + 1}
  iterations.

  Temporary result dumping may be controlled by \code{\link{.CF.DP}}.
}
\value{
  A list contains three big lists of MCMC traces including:
  \code{b.Mat} for mutation and selection coefficients of \code{\link{b}},
  \code{p.Mat} for hyper-parameters, and
  \code{phi.Mat} for expected expression values Phi.
  All lists are of length \code{(burnin + nIter) / iterThin + 1} and
  each element contains the output of each iteration.

  All lists also can be binded as trace matrices, such as via
  \code{do.call("rbind", b.Mat)} yielding a matrix of dimension number of
  iterations by number of parameters. Then, those traces can be analyzed
  further via other MCMC packages such as \pkg{coda}.
}
\note{
  Note that \code{phi.Obs} need to be normalized to mean 1.
}
\references{
  \url{http://www.tiem.utk.edu/~mikeg/}

  Shah P. and Gilchrist MA.
    ``Explaining complex codon usage patterns
      with selection for translational efficiency, mutation bias, and genetic
      drift''
    Proc Natl Acad Sci USA (2011) 108:10231--10236.

  Wallace E.W.J., Airoldi E.M., and Drummond D.A.
    ``Estimating Selection on Synonymous Codon Usage from Noisy Experimental
      Data''
    Mol Biol Evol (2013) 30(6):1438--1453.
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}.
}
\seealso{
  \link{DataIO}, \link{DataConverting},
  \code{\link{cubappr}()} and \code{\link{cubpred}()}.
}
\examples{
\dontrun{
library(cubfits)

demo(roc.train, 'cubfits', ask = F, echo = F)
}
}
\keyword{main function}
