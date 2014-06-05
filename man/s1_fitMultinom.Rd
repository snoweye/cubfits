\name{Fit Multinomial}
\alias{fitMultinom}
\title{Fit Multinomial Model (Generic)}
\description{
  This generic function estimates \code{\link{b}}
  (mutation (log(mu)) and selection (Delta.t) parameters)
  depending on options set by \code{\link{init.function}()}.
}
\usage{
  fitMultinom(reu13.df, phi, y, n, phi.new = NULL, coefstart = NULL)
}
\arguments{
  \item{reu13.df}{an object of format \code{\link{reu13.df}}.}
  \item{phi}{an object of format \code{\link{phi.Obs}}.}
  \item{y}{an object of format \code{\link{y}}.}
  \item{n}{an object of format \code{\link{n}}.}
  \item{phi.new}{an object of format \code{\link{phi.Obs}} for MCMC only.}
  \item{coefstart}{initial value for \code{\link{b}}
                   (mutation (log(mu)) and selection (Delta.t) parameters)
                   only used in \code{vglm()}.}
}
\details{
  \code{fitMultinom()} fits a multinomial logistic regression via
  vector generalized linear model fitting, \code{vglm()}.
  By default, for each amino acids, the last codon (order by characters)
  is assumed as a based line, and other codons are compared to the based
  line relatively.

  In MCMC, \code{phi.new} are new proposed expression values and
  used to propose new \code{\link{b}}. The \code{coefstart} is used to avoid
  randomization of estimating \code{\link{b}} in \code{vglm()},
  and speed up computation.
}
\value{
  A list of format \code{\link{b}} is returned which are modified from
  the returns of \code{vglm()}. Mainly, it includes
  \code{b$coefficient} (parameters in \code{vector}),
  \code{b$coef.mat} (parameters in \code{matrix}), and
  \code{b$R} (covariance matrix of parameters, *R* matrix in QR decomposition).
}
\references{
  \url{https://github.com/snoweye/cubfits/}

  Shah P. and Gilchrist M.A.
    ``Explaining complex codon usage patterns
      with selection for translational efficiency, mutation bias, and
      genetic drift''
    Proc Natl Acad Sci USA (2011) 108:10231--10236.
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}.
}
\seealso{
  \code{\link{init.function}()} and \code{\link{estimatePhi}()}.
}
\examples{
\dontrun{
suppressMessages(library(cubfits, quietly = TRUE))
set.seed(1234)

# Convert data.
reu13.list <- convert.reu13.df.to.list(ex.test$reu13.df)
y.list <- convert.y.to.list(ex.test$y)
n.list <- convert.n.to.list(ex.test$n)

# Get phi.pred.Init
init.function(model = "roc")
fitlist <- fitMultinom(ex.train$reu13.df, ex.train$phi.Obs, ex.train$y, ex.train$n)
phi.pred.Init <- estimatePhi(fitlist, reu13.list, y.list, n.list,
                         E.Phi = median(ex.test$phi.Obs),
                         lower.optim = min(ex.test$phi.Obs) * 0.9,
                         upper.optim = max(ex.test$phi.Obs) * 1.1)
}
}
\keyword{tool}
