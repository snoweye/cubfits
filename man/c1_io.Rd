\name{Cedric IO Utilities}
\alias{readGenome}
\alias{normalizeDataSet}
\title{Cedric IO Utilities}
\description{
  This utility function provides basic IO by Cedric.
}
\usage{
  readGenome(fn.genome, ex.sh.aa = 0, rm.first.aa = 0)
  
  normalizeDataSet(data)
}
\arguments{
  \item{fn.genome}{Fasta file with sequences}
  \item{ex.sh.aa}{Ignore sequences with a length less than ex.sh.aa. (After removal of the first rm.first.aa amino acids)}
  \item{rm.first.aa}{Remove the first rm.first.aa amino acids (after start codon)}
  \item{data}{Vector to be normalized. Means will be set to 1}
}
\details{
  under development
}
\value{
  under development
}
\references{
  \url{https://github.com/clandere/cubfits/}
}
\author{
  Cedric Landerer \email{cedric.landerer@gmail.com}.
}
\seealso{
  under development
}
\examples{
  \dontrun{
    library(cubfits)
    seq.string <- readGenome("my_genome.fasta", 150, 10)
  }
}
\keyword{cedric}
