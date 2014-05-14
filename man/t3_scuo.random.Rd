\name{Randomize SCUO Index}
\alias{scuo.random}
\title{Generate Randomized SCUO Index}
\description{
  Generate randomized SCUO indices in log normal distribution, but
  provided original unchanged SCUO order.
}
\usage{
  scuo.random(SCUO, phi.Obs = NULL, meanlog = .CF.PARAM$phi.meanlog,
              sdlog = .CF.PARAM$phi.sdlog)
}
\arguments{
  \item{SCUO}{SCUO index returned from \code{\link{calc_scuo_values}()}.}
  \item{phi.Obs}{optional object of format \code{\link{phi.Obs}}.}
  \item{meanlog}{mean of log normal distribution.}
  \item{sdlog}{std of log normal distribution.}
}
\details{
  This function takes \code{SCUO} indices (outputs of
  \code{\link{calc_scuo_values}()})
  computes the rank of them, generates log normal random variables, and
  replaces \code{SCUO} indices by those variables in the same rank orders.
  Typically, these random variables are used to replace expression values
  when either no expression is observed or for the purpose of model validation.

  If \code{phi.Obs} is provided, the mean and std of \code{log(phi.Obs)} are used
  for log normal random variables. Otherwise, \code{menalog} and \code{sdlog}
  are used.

  The default \code{meanlog} and \code{sdlog} was estimated from
  \code{yassour} dataset.
}
\value{
  A vector of log normal random variables is returned.
}
\references{
  \url{https://github.com/snoweye/cubfits/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}.
}
\seealso{
  \code{\link{calc_scuo_values}()}, \code{yassour}.
}
\examples{
\dontrun{
suppressMessages(library(cubfits, quietly = TRUE))

### example dataset.
y.scuo <- convert.y.to.scuo(ex.train$y)
SCUO <- calc_scuo_values(y.scuo)$SCUO
plotprxy(ex.train$phi.Obs, SCUO)

### yassour dataset.
GM <- apply(yassour[, -1], 1, function(x) exp(mean(log(x[x != 0]))))
phi.Obs <- GM / sum(GM) * 15000
mean(log(phi.Obs))
sd(log(phi.Obs))
ret <- scuo.random(SCUO, meanlog = -0.441473, sdlog = 1.393285)
plotprxy(ret, SCUO)
}
}
\keyword{tool}
