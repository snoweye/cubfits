\name{Plotbin}
\alias{prop.bin.roc}
\alias{plotbin}
\title{Plot Binning Results}
\description{
  Plot binning results to visualize the effects of mutation and selection
  along with expression levels empirically.
}
\usage{
  prop.bin.roc(reu13.df, phi.Obs = NULL, nclass = 20)

  plotbin(ret.bin, ret.model = NULL, main = NULL,
          xlab = "Production Rate (log10)", ylab = "Proportion",
          xlim = NULL, lty = 1, x.log10 = TRUE, stderr = FALSE, ...)
}
\arguments{
  \item{reu13.df}{a \code{\link{reu13.df}} object.}
  \item{phi.Obs}{a \code{\link{phi.Obs}} object.}
  \item{nclass}{number of binning classes across the range of \code{phi.Obs}.}

  \item{ret.bin}{binning results from \code{\link{prop.bin.roc}()}.}
  \item{ret.model}{model results from \code{\link{prop.model.roc}()}.}
  \item{main}{an option passed to \code{plot()}.}
  \item{xlab}{an option passed to \code{plot()}.}
  \item{ylab}{an option passed to \code{plot()}.}
  \item{xlim}{range of X-axis.}
  \item{lty}{line type if \code{ret.model} is provided.}
  \item{x.log10}{\code{log10()} transformation of X-axis.}
  \item{stderr}{plot stand error instead of stand deviation.}
  \item{...}{options passed to \code{plot()}.}
}
\details{
  The function \code{plotbin()} plots the binning results \code{ret.bin}
  returned from \code{prop.bin.roc()}. Fitted curves may be added if
  \code{ret.model} is provided which can be obtained from
  \code{\link{prop.model.roc}()}.

  \code{\link{plotaddmodel}()} can append model later if \code{ret.model}
  is not provided to \code{plotbin()}.

  Currently, only ROC model is supported.
  Colors are controlled by \code{\link{.CF.PT}}.
}
\value{
  A binning plot is drawn.
}
\references{
  \url{https://github.com/snoweye/cubfits/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}.
}
\seealso{
  \code{\link{plotmodel}()} and \code{\link{prop.model.roc}()}.
}
\examples{
\dontrun{
demo(plotbin, 'cubfits', ask = F, echo = F)
}
}
\keyword{plotting}
