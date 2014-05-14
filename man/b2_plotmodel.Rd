\name{Plotmodel}
\alias{prop.model.roc}
\alias{plotmodel}
\alias{plotaddmodel}
\title{Plot Fitted Models} 
\description{
  Plot model results to visualize the effects of mutation and selection
  along with expression levels. 
  The model can be fitted by MCMC or multinomial logistic regression.
}
\usage{
  prop.model.roc(bInit, phi.Obs.lim = c(0.01, 10), phi.Obs.scale = 1,
                 nclass = 40, x.log10 = TRUE)

  plotmodel(ret.model, main = NULL,
            xlab = "Production Rate (log10)", ylab = "Proportion",
            xlim = NULL, lty = 1, x.log10 = TRUE, ...)

  plotaddmodel(ret.model, lty, u.codon = NULL, color = NULL,
               x.log10 = TRUE)
}
\arguments{
  \item{bInit}{a \code{\link{b}} object.}
  \item{phi.Obs.lim}{range of \code{phi.Obs}.}
  \item{phi.Obs.scale}{optional scaling factor.}
  \item{nclass}{number of binning classes across the range of \code{phi.Obs}.}
  \item{x.log10}{\code{log10()} transformation of X-axis.}

  \item{ret.model}{model results from \code{\link{prop.model.roc}()}.}
  \item{main}{an option passed to \code{plot()}.}
  \item{xlab}{an option passed to \code{plot()}.}
  \item{ylab}{an option passed to \code{plot()}.}
  \item{xlim}{range of X-axis.}
  \item{lty}{line type.}

  \item{u.codon}{unique synonymous codon names.}
  \item{color}{a color vector for unique codon, typically returns of
               the internal function \code{get.color()}.}

  \item{...}{options passed to \code{plot()}.}
}
\details{
  The function \code{plotmodel()} plots the fitted curves obtained from
  \code{\link{prop.model.roc}()}.

  The function \code{plotaddmodel()} can append model curves to a binning plot
  provided unique synonymous codons and colors are given. This function is
  nearly for an internal call within \code{plotmodel()}, but is exported and
  useful for workflow.

  Currently, only ROC model is supported.
  Colors are controlled by \code{\link{.CF.PT}}.
}
\value{
  A fitted curve plot is drawn.
}
\references{
  \url{https://github.com/snoweye/cubfits/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}.
}
\seealso{
  \code{\link{plotbin}()}, \code{\link{prop.bin.roc}()}, and
  \code{\link{prop.model.roc}()}.
}
\examples{
\dontrun{
demo(plotbin, 'cubfits', ask = F, echo = F)
}
}
\keyword{plotting}
