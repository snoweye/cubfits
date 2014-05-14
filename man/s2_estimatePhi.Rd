\name{Estimate Phi}
\alias{estimatePhi}
\title{Initialization of Phi (Generic)}
\description{
  This generic function estimates Phi (expression value) either by posterior
  mean (PM) or by maximum likelihood estimator (MLE) depending on options set
  by \code{\link{init.function}()}.
}
\usage{
  estimatePhi(fitlist, reu13.list, y.list, n.list,
    E.Phi = .CF.OP$E.Phi, lower.optim = .CF.OP$lower.optim,
    upper.optim = .CF.OP$upper.optim,
    lower.integrate = .CF.OP$lower.integrate,
    upper.integrate = .CF.OP$upper.integrate, control = list())
}
\arguments{
  \item{fitlist}{an object of format \code{\link{b}}.}
  \item{reu13.list}{an object of format \code{\link{reu13.list}}.}
  \item{y.list}{an object of format \code{\link{y.list}}.}
  \item{n.list}{an object of format \code{\link{n.list}}.}
  \item{E.Phi}{potential expected value of Phi.}
  \item{lower.optim}{lower bound to \code{optim()}.}
  \item{upper.optim}{upper bound to \code{optim()}.}
  \item{lower.integrate}{lower bound to \code{integrate()}.}
  \item{upper.integrate}{upper bound to \code{integrate()}.}
  \item{control}{control options to \code{optim()}.}
}
\details{
  \code{estimatePhi()} is a generic function first initialized by
  \code{\link{init.function}()}, then it estimates Phi accordingly.
  By default, \code{.CF.CT$init.Phi} sets the method \code{PM} for the
  posterior mean.

  \code{PM} uses a flat prior and \code{integrate()} to estimate
  Phi. While, \code{MLE} uses \code{optim()} to estimate Phi which
  may have boundary solutions for some sequences.
}
\value{
  Estimated Phi for every sequence is returned.
}
\references{
  \url{https://github.com/snoweye/cubfits/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}.
}
\seealso{
  \code{\link{init.function}()} and \code{\link{fitMultinom}()}.
}
\examples{
\dontrun{
suppressMessages(library(cubfits, quietly = TRUE))
set.seed(1234)

# Convert data.
reu13.list <- convert.reu13.df.to.list(ex.test$reu13.df)
y.list <- convert.y.to.list(ex.test$y)
n.list <- convert.n.to.list(ex.test$n)

# Get phi.Init.pred
init.function(model = "roc")
fitlist <- fitMultinom(ex.train$reu13.df, ex.train$phi.Obs, ex.train$y, ex.train$n)
phi.Init.pred <- estimatePhi(fitlist, reu13.list, y.list, n.list,
                         E.Phi = median(ex.test$phi.Obs),
                         lower.optim = min(ex.test$phi.Obs) * 0.9,
                         upper.optim = max(ex.test$phi.Obs) * 1.1)
}
}
\keyword{tool}
