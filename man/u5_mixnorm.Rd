\name{Mixed Normal Optimization}
\alias{mixnormerr.optim}
\alias{dmixnormerr}
\title{Mixed Normal Optimization}
\description{
  Constrained optimization for mixed normal in 1D and typically for 2
  components.
}
\usage{
  mixnormerr.optim(X, K = 2, param = NULL)
  dmixnormerr(x, param)
}
\arguments{
  \item{X}{a gene expression data matrix of dimension \code{N * R}
           which has \code{N} genes and \code{R} replicates.}
  \item{K}{number of components to fit.}

  \item{x}{vector of quantiles.}
  \item{param}{parameters of \code{mixnormerr}, typically the element
               \code{param} of the \code{mixnormerr.optim()} returning object.}
}
\details{
  The function \code{mixnormerr.optim()}
  maximizes likelihood using \code{constrOptim()} based on
  the gene expression data \code{X} (usually in log scale)
  for \code{N} genes and \code{R} replicates (\code{NA} is allowed).
  The likelihood of each gene expression
  is a \code{K = 2} component mixed normal distribution
  (\eqn{\sum_k p_k N(mu_k, \sigma_k^2 + \sigma_e^2)}{
        sum_k p_k N(mu_k, sigma_k^2 + sigma_e^2)})
  with measurement errors of the replicates
  (\eqn{N(0, \sigma_e^2)}{N(0, sigma_e^2)}).

  The \code{sigma_k^2} is as the error of random component and
  the \code{sigma_e^2} is as the error of fixed component. Both
  are within a mixture model of two normal distributions.

  The function \code{dmixnormerr()} computes the density of the mixed
  normal distribution.

  \code{param} is a parameter list and contains five elements:
  \code{K} for number of components,
  \code{prop} for proportions,
  \code{mu} for centers of components,
  \code{sigma2} for variance of components, and
  \code{sigma2.e} for variance of measurement errors. 
}
\value{
  \code{mixnormerr.optim()} returns a list containing three main elements
  \code{param} is the final results (MLEs), \code{param.start} is the starting
  parameters, and \code{optim.ret} is the original returns of
  \code{constrOptim()}.
}
\note{
  This function is limited for small \code{K}. An equivalent EM algorithm
  should be done in a more stable way for large \code{K}.
}
\references{
  \url{https://github.com/snoweye/cubfits/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}.
}
\seealso{
  \code{\link{print.mixnormerr}()},
  \code{\link{simu.mixnormerr}()}.
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
\keyword{utility}
