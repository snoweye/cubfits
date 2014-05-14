\name{Datasets}
\alias{bInit}
\alias{ex.test}
\alias{ex.train}
\docType{data}
\title{Datasets for Demonstrations}
\description{
  Examples of toy data to test and demonstrate \pkg{cubfits}.
}
\usage{
  bInit
  ex.test
  ex.train
}
\format{
  All are in list formats.
}
\details{
  \code{bInit} contains two sets (\code{roc} and \code{rocnse}) of
  initial coefficients including mutation and selection parameters for
  3 amino acids 'A', 'C', and 'D' in \code{matrix} format.
  Both sets are in \code{\link{b}} format.

  \code{ex.train} contains a training set of 100 sequences including
  3 \code{\link{reu13.df}} (codon counts in \code{reu13} data frame format
  divided by amino acids),
  3 \code{\link{y}} (codon counts in simplified data frame format
  divided by amino acids),
  3 \code{\link{n}} (total amino acid counts in vector format
  divided by amino acids), and
  \code{\link{phi.Obs}} (observed phi values in vector format).

  \code{ex.test} contains a testing set of the other 100 sequences in the
  same format of \code{ex.train}.
}
\references{
  \url{https://github.com/snoweye/cubfits/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}.
}
\seealso{
  \code{\link{init.function}()}, \code{\link{cubfits}()},
  \code{\link{cubpred}()}, and \code{\link{cubappr}()}.
}
\examples{
\dontrun{
suppressMessages(library(cubfits, quietly = TRUE))

str(bInit)
str(ex.test)
str(ex.train)
}
}
\keyword{datasets}
