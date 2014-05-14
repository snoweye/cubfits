\name{Input and Output Utility}
\alias{get.expath}
\alias{read.seq}
\alias{read.phi.df}
\alias{write.seq}
\alias{write.phi.df}
\alias{DataIO}
\title{Input and Output Utility}
\description{
  These utility functions read and write data of FASTA and phi.df formats.
}
\usage{
  read.seq(file.name, forceDNAtolower = FALSE, convertDNAtoupper = TRUE)
  write.seq(seq.data, file.name)

  read.phi.df(file.name, header = TRUE, sep = "\t", quote = "")
  write.phi.df(phi.df, file.name)

  get.expath(file.name, path.root = "./ex_data/", pkg = "cubfits")
}
\arguments{
  \item{file.name}{a file name to read or write.}
  \item{forceDNAtolower}{an option passed to \code{read.fasta()} of
                         \pkg{seqinr} package.}
  \item{convertDNAtoupper}{force everything in upper case.}
  \item{header}{an option passed to \code{read.table()}.}
  \item{sep}{an option passed to \code{read.table()}.}
  \item{quote}{an option passed to \code{read.table()}.}
  \item{seq.data}{a \code{\link{seq.data}} object.}
  \item{phi.df}{a \code{\link{phi.df}} object.}
  \item{path.root}{root path for the file name relatively to the pkg.}
  \item{pkg}{package name for the path of root.}
}
\details{
  \code{read.seq()} and \code{write.seq()} typically read and
  write FASTA files (DNA ORFs or sequences).

  \code{read.phi.df()} and \code{write.phi.df()} typically read and write
  phi.df files (expression values of ORFs or sequences).

  \code{get.expath()} is only for demonstration returning a full path
  to the file.
}
\value{
  \code{read.seq()} returns an object of \code{\link{seq.data}} format
  which can be converted to \code{\link{seq.string}} format later via
  \code{\link{convert.seq.data.to.string}()}.

  \code{read.phi.df()} returns an object of \code{\link{phi.df}} format
  which contains expression values.
}
\references{
  \url{https://github.com/snoweye/cubfits/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}.
}
\seealso{
  \code{\link{convert.seq.data.to.string}()}.
}
\examples{
\dontrun{
suppressMessages(library(cubfits, quietly = TRUE))

seq.data <- read.seq(get.expath("seq_200.fasta"))
phi.df <- read.phi.df(get.expath("phi_200.tsv"))
aa.names <- c("A", "C", "D")

# Read in from FASTA file.
seq.string <- convert.seq.data.to.string(seq.data)
}
}
\keyword{utility}
