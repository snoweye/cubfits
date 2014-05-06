\name{Asymmetric Laplace Distribution}
\alias{dasl}
\alias{dasla}
\alias{pasl}
\alias{pasla}
\alias{qasl}
\alias{qasla}
\alias{rasl}
\alias{rasla}
\alias{asl.optim}
\title{The Asymmetric Laplace Distribution}
\description{
  Density, probability, quantile, random number generation, and MLE functions
  for the asymmetric Laplace distribution
  with parameters either in \eqn{ASL(\theta, \mu, \sigma)}{ASL(theta, mu, sigma)}
  or the alternative
  \eqn{ASL^*(\theta, \kappa, \sigma)}{ASL*(theta, kappa, sigma)}.
}
\usage{
  dasl(x, theta = 0, mu = 0, sigma = 1, log = FALSE)
  dasla(x, theta = 0, kappa = 1, sigma = 1, log = FALSE)

  pasl(q, theta = 0, mu = 0, sigma = 1, lower.tail = TRUE,
          log.p = FALSE)
  pasla(q, theta = 0, kappa = 1, sigma = 1, lower.tail = TRUE,
           log.p = FALSE)

  qasl(p, theta = 0, mu = 0, sigma = 1, lower.tail = TRUE,
          log.p = FALSE)
  qasla(p, theta = 0, kappa = 1, sigma = 1, lower.tail = TRUE,
           log.p = FALSE)

  rasl(n, theta = 0, mu = 0, sigma = 1)
  rasla(n, theta = 0, kappa = 1, sigma = 1)

  asl.optim(x)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
           taken to be the number required.}
  \item{theta}{center parameter.}
  \item{mu, kappa}{location parameters.}
  \item{sigma}{shape parameter.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities p are given as
                    log(p).}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
                    P[X <= x] otherwise, P[X > x].}
}
\details{
  The density \eqn{f(x)}{f(x)} of
  \eqn{ASL^*(\theta, \kappa, \sigma)}{ASL*(theta, kappa, sigma)}
  is given as
  \eqn{\frac{\sqrt{2}}{\sigma}\frac{\kappa}{1 + \kappa^2}
       exp(- \frac{\sqrt{2}\kappa}{\sigma} |x - \theta|)}{
       sqrt(2) / sigma kappa / (1 + \kappa^2)
       exp(- sqrt(2) kappa / sigma |x - \theta|)
      }
  if \eqn{x \ge \theta}{x >= theta}, and
  \eqn{\frac{\sqrt{2}}{\sigma}\frac{\kappa}{1 + \kappa^2}
       exp(- \frac{\sqrt{2}}{\sigma\kappa} |x - \theta|)}{
       sqrt(2) / sigma kappa / (1 + \kappa^2)
       exp(- sqrt(2) / (sigma kappa) |x - \theta|)
      }
  if \eqn{x < \theta}{x < theta}.

  The parameter domains of ASL and ASL* are
  \eqn{\theta \in R}{theta in real},
  \eqn{\sigma > 0}{sigma > 0},
  \eqn{\kappa > 0}{kappa > 0}, and
  \eqn{\mu \in R}{mu in real}.
  The relation of \eqn{\mu}{mu} and \eqn{\kappa}{kappa} are
  \eqn{\kappa = \frac{\sqrt{2\sigma^2 + \mu^2}-\mu}{\sqrt{2\sigma}}}{
       kappa = (sqrt(2 sigma^2 + mu^2) - mu) / sqrt(2 sigma)}
  or
  \eqn{\mu = \frac{\sigma}{\sqrt{2}} (\frac{1}{\kappa} - \kappa)}{
       mu = sigma / sqrt(2) (1 / kappa - kappa)}.
}
\value{
  \dQuote{dasl} and \dQuote{dasla} give the densities,
  \dQuote{pasl} and \dQuote{pasla} give the distribution functions,
  \dQuote{qasl} and \dQuote{qasla} give the quantile functions, and
  \dQuote{rasl} and \dQuote{rasls} give the random numbers.

  \code{asl.optim} returns the MLE of data \code{x} including
  \code{theta}, \code{mu}, \code{kappa}, and \code{sigma}.
}
\references{
  Kotz S, Kozubowski TJ, Podgorski K. (2001) ``The Laplace distribution
    and generalizations: a revisit with applications to
    communications, economics, engineering, and finance.''
    Boston: Birkhauser.
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}.
}
\examples{
\dontrun{
suppressMessages(library(cubfits, quietly = TRUE))
set.seed(1234)

dasl(-2:2)
dasla(-2:2)
pasl(-2:2)
pasla(-2:2)
qasl(seq(0, 1, length = 5))
qasla(seq(0, 1, length = 5))

dasl(-2:2, log = TRUE)
dasla(-2:2, log = TRUE)
pasl(-2:2, log.p = TRUE)
pasla(-2:2, log.p = TRUE)
qasl(log(seq(0, 1, length = 5)), log.p = TRUE)
qasla(log(seq(0, 1, length = 5)), log.p = TRUE)

set.seed(123)
rasl(5)
rasla(5)

asl.optim(rasl(5000))
}
}
\keyword{utility}
