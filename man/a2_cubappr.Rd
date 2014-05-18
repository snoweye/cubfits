\name{CUB Model Approximation}
\alias{cubappr}
\title{Codon Usage Bias Approximation for ORFs without Expression}
\description{
  This function provides codon usage bias approximation with observed ORFs
  but without any expressions.
}
\usage{
  cubappr(reu13.df.obs, phi.Init, y, n,
          nIter = 1000, burnin = 100,
          bInit = NULL, init.b.Scale = .CF.CONF$init.b.Scale,
              b.DrawScale = .CF.CONF$b.DrawScale,
          p.Init = NULL, p.nclass = .CF.CONF$p.nclass,
              p.DrawScale = .CF.CONF$p.DrawScale,
          phi.DrawScale.pred = .CF.CONF$phi.DrawScale.pred,
          model = .CF.CT$model[1], model.Phi = .CF.CT$model.Phi[1],
          adaptive = .CF.CT$adaptive[1],
          verbose = .CF.DP$verbose,
          iterThin = .CF.DP$iterThin, report = .CF.DP$report)
}
\arguments{
  \item{reu13.df.obs}{a \code{\link{reu13.df}} object, ORFs information.}
  \item{phi.Init}{a \code{\link{phi.Obs}} object, temporarily initial of expression
               without measurement errors.}
  \item{y}{a \code{\link{y}} object, codon counts.}
  \item{n}{a \code{\link{n}} object, total codon counts.}
  \item{nIter}{number of iterations after burn-in iterations.}
  \item{burnin}{number of burn-in iterations.}
  \item{bInit}{initial values for parameters \code{\link{b}}.}
  \item{init.b.Scale}{for initial \code{\link{b}} if \code{bInit = NULL}.}
  \item{b.DrawScale}{scaling factor for adaptive MCMC with random walks
                     when drawing new \code{\link{b}}.}
  \item{p.Init}{initial values for hyper-parameters.}
  \item{p.nclass}{number of components for \code{model.Phi = "logmixture"}.}
  \item{p.DrawScale}{scaling factor for adaptive MCMC with random walks
                     when drawing new \code{sigma.Phi}.}
  \item{phi.DrawScale.pred}{scaling factor for adaptive MCMC with random walks
                            when drawing new Phi of predicted set.}
  \item{model}{model to be fitted, currently "roc" only.}
  \item{model.Phi}{prior model for Phi, currently "lognormal".}
  \item{adaptive}{adaptive method of MCMC for proposing new \code{\link{b}} and
                  Phi.}
  \item{verbose}{print iteration messages.}
  \item{iterThin}{thinning iterations.}
  \item{report}{number of iterations to report more information.}
}
\details{
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
  Note that \code{phi.Init} need to be normalized to mean 1.

  \code{p.DrawScale} may cause scaling prior if adaptive MCMC is used, and
  it can result in non-exits of equilibrium distribution.
}
\references{
  \url{https://github.com/snoweye/cubfits/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}.
}
\seealso{
  \link{DataIO}, \link{DataConverting},
  \code{\link{cubfits}()} and \code{\link{cubpred}()}.
}
\examples{
\dontrun{
suppressMessages(library(cubfits, quietly = TRUE))

demo(roc.appr, 'cubfits', ask = F, echo = F)
}
}
\keyword{main function}
