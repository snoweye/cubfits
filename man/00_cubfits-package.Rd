\name{cubfits-package}
\alias{cubfits-package}
\docType{package}
\title{Codon Bias Usage Fits}
\description{
  Estimating mutation and selection coefficients on synonymous
  codon bias usage based on models of ribosome overhead cost (ROC).
  Multinomial logistic regression and Markov Chain Monte Carlo are used to
  estimate and predict protein production rates with/without the presence of
  expressions and measurement errors.
}
\details{
\tabular{ll}{
  Package: \tab cubfits\cr
  Type: \tab Package\cr
  License: \tab Mozilla Public License 2.0\cr
  LazyLoad: \tab yes\cr
}

  The install command is simply as \cr
  \cr
  \code{> R CMD INSTALL cubfits_*.tar.gz } \cr
  \cr
  from a command mode or \cr
  \cr
  \code{R> install.packages("cubfits") } \cr
  \cr
  inside an R session.
}
\references{
  \url{https://github.com/snoweye/cubfits/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}, Russell Zaretzki,
  William Howell, Drew Schmidt, and Michael Gilchrist.
}
\seealso{
  \code{\link{init.function}()}, \code{\link{cubfits}()},
  \code{\link{cubpred}()}, and \code{\link{cubappr}()}.
}
\examples{
\dontrun{
suppressMessages(library(cubfits, quietly = TRUE))

demo(roc.train, 'cubfits', ask = F, echo = F)
demo(roc.pred, 'cubfits', ask = F, echo = F)
demo(roc.appr, 'cubfits', ask = F, echo = F)
}
}
\keyword{package}

