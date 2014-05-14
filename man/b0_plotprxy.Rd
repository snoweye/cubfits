\name{Plotprxy}
\alias{plotprxy}
\title{Predictive X-Y Plot} 
\description{
  This utility function provides a basic plot of production rates.
}
\usage{
  plotprxy(x, y, x.ci = NULL, y.ci = NULL,
           log10.x = TRUE, log10.y = TRUE,
           add.lm = TRUE, add.one.to.one = TRUE, weights = NULL,
           add.legend = TRUE,
           xlim = NULL, ylim = NULL,
           xlab = "Predicted Production Rate (log10)",
           ylab = "Observed Production Rate (log10)",
           main = NULL)
}
\arguments{
  \item{x}{expression values.}
  \item{y}{expression values, of the same length of \code{x}.}
  \item{x.ci}{confidence interval of \code{x}, of dimension
              \code{length{x} * 2}, for outliers labeling.}
  \item{y.ci}{confidence interval of \code{y}, of dimension
              \code{length{y} * 2}, for outliers labeling.}
  \item{log10.x}{\code{log10()} and mean transformation of x axis.}
  \item{log10.y}{\code{log10()} and mean transformation of y axis.}
  \item{add.lm}{if add \code{lm()} fit.}
  \item{add.one.to.one}{if add one-to-one line.}
  \item{weights}{weights to \code{lm()}.}
  \item{add.legend}{if add default legend.}
  \item{xlim}{limits of x-axis.}
  \item{ylim}{limits of y-axis.}
  \item{xlab}{an option passed to \code{plot()}.}
  \item{ylab}{an option passed to \code{plot()}.}
  \item{main}{an option passed to \code{plot()}.}
}
\details{
  As the usual X-Y plot where \code{x} and \code{y} are expression values.

  If \code{add.lm = TRUE} and \code{weights} are given, then both ordinary
  and weighted least squares results will be plotted.
}
\value{
  A scatter plot with a fitted \code{lm()} line and R squared value.
}
\references{
  \url{https://github.com/snoweye/cubfits/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}.
}
\seealso{
  \code{\link{plotbin}()} and \code{\link{plotmodel}()}.
}
\examples{
\dontrun{
suppressMessages(library(cubfits, quietly = TRUE))

y.scuo <- convert.y.to.scuo(ex.train$y)
SCUO <- calc_scuo_values(y.scuo)$SCUO
plotprxy(ex.train$phi.Obs, SCUO)
}
}
\keyword{plotting}
