\name{SCUO Index}
\alias{calc_scuo_values}
\title{Function for Synonymous Codon Usage Order (SCUO) Index}
\description{
  Calculate the Synonymous Codon Usage Order (SCUO) index for each gene.
  Used as a substitute for expression in cases of without expression
  measurements.
}
\usage{
  calc_scuo_values(codon.counts)
}
\arguments{
  \item{codon.counts}{an object of format \code{\link{scuo}}.}
}
\details{
  This function computes SCUO index for each gene. Typically, this method
  is completely based on entropy and information theory to estimate
  expression values of sequences according to their codon information.
}
\value{
  \code{SCUO} indices are returned.
}
\references{
  \url{https://www.tandfonline.com/doi/abs/10.1080/03081070500502967}

  Wan X.-F., Zhou J., Xu D.
    ``CodonO: a new informatics method for measuring synonymous codon usage
      bias within and across genomes''
    International Journal of General Systems Vol. 35, Iss. 1, 2006.
}
\author{
  Drew Schmidt.
}
\seealso{
  \code{\link{scuo.random}()}, \code{\link{calc_cai_values}()},
  \code{\link{calc_scu_values}()}.
}
\examples{
\dontrun{
suppressMessages(library(cubfits, quietly = TRUE))

y.scuo <- convert.y.to.scuo(ex.train$y)
SCUO <- calc_scuo_values(y.scuo)$SCUO
plotprxy(ex.train$phi.Obs, SCUO, ylab = "SCUO (log10)")
}
}
\keyword{tool}
