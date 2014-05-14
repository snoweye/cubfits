\name{Print}
\alias{print.mixnormerr}
\title{Functions for Printing Objects According to Classes}
\description{
  A Class \code{mixnormerr} is declared in \pkg{cubfits}, and this is the function
  to print and summary objects.
}
\usage{
\method{print}{mixnormerr}(x, digits = max(4, getOption("digits") - 3), ...)
}
\arguments{
  \item{x}{an object with the class attributes.}
  \item{digits}{for printing out numbers.}
  \item{\dots}{other possible options.}
}
\details{
  This is an useful function for summarizing and debugging.
}
\value{
  The results will cat or print on the STDOUT by default.
}
\references{
  \url{https://github.com/snoweye/cubfits/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}.
}
\seealso{
  \code{\link{mixnormerr.optim}()}.
}
\examples{
\dontrun{
suppressMessages(library(cubfits, quietly = TRUE))

### Get individual of phi.Obs.
GM <- apply(yassour[, -1], 1, function(x) exp(mean(log(x[x != 0]))))
phi.Obs.all <- yassour[, -1] / sum(GM) * 15000
phi.Obs.all[phi.Obs.all == 0] <- NA

### Run optimization.
X <- log(as.matrix(phi.Obs.all))
param.init <- list(K = 2, prop = c(0.95, 0.05), mu = c(-0.59, 3.11),
                   sigma2 = c(1.40, 0.59), sigma2.e = 0.03)
ret <- mixnormerr.optim(X, K = 2, param = param.init)
print(ret)
}
}
\keyword{summary}

