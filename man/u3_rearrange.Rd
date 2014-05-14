\name{Rearrangment Utility}
\alias{rearrange.reu13.df}
\alias{rearrange.y}
\alias{rearrange.n}
\alias{rearrange.phi.Obs}
\title{Rearrange Data Structure by ORF Names}
\description{
  These utility functions rearrange data in the order of ORF names.
}
\usage{
  rearrange.reu13.df(reu13.df)
  rearrange.y(y)
  rearrange.n(n)
  rearrange.phi.Obs(phi.Obs)
}
\arguments{
  \item{reu13.df}{a list of \code{\link{reu13.df}} data frames
                  divided by amino acids.}
  \item{y}{a list of \code{\link{y}} data frames divided by amino acids.}
  \item{n}{a list of \code{\link{n}} vectors divided by amino acids.}
  \item{phi.Obs}{a vector of \code{\link{phi.Obs}} format.}
}
\details{
  These utility functions take inputs and return ordered outputs.
  It is necessary to rearrange data in a right order of ORF names which
  avoids subsetting data frame within MCMC and improve performance.
}
\value{
  The outputs are in the same format of inputs except the order of data
  is sorted by ORF names.
}
\references{
  \url{https://github.com/snoweye/cubfits/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}.
}
\seealso{
  \link{AllDataFormats},
  \code{\link{convert.n.to.list}()},
  \code{\link{convert.reu13.df.to.list}()}, and
  \code{\link{convert.y.to.list}()}.
}
\examples{
\dontrun{
suppressMessages(library(cubfits, quietly = TRUE))

reu13.df <- rearrange.reu13.df(ex.train$reu13.df)
y <- rearrange.y(ex.train$y)
n <- rearrange.n(ex.train$n)
phi.Obs <- rearrange.phi.Obs(ex.train$phi.Obs)
}
}
\keyword{utility}
