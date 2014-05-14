\name{Coverting Utility}
\alias{convert.reu13.df.to.list}
\alias{convert.y.to.list}
\alias{convert.n.to.list}
\alias{convert.y.to.scuo}
\alias{convert.seq.data.to.string}
\alias{codon.low2up}
\alias{codon.up2low}
\alias{dna.low2up}
\alias{dna.up2low}
\alias{convert.b.to.bVec}
\alias{convert.bVec.to.b}
\alias{DataConverting}
\title{Convert Data Frame to Other Formats}
\description{
  These utility functions convert data of format divided by amino acids into
  list of format divided by ORFs, or convert data to other formats.
}
\usage{
  convert.reu13.df.to.list(reu13.df)
  convert.y.to.list(y)
  convert.n.to.list(n)

  convert.y.to.scuo(y)
  convert.seq.data.to.string(seq.data)

  codon.low2up(x)
  codon.up2low(x)

  dna.low2up(x)
  dna.up2low(x)

  convert.b.to.bVec(b)
  convert.bVec.to.b(bVec, aa.names, model = .CF.CT$model[1])
}
\arguments{
  \item{reu13.df}{a list of \code{\link{reu13.df}} data frames
                  divided by amino acids.}
  \item{y}{a list of \code{\link{y}} data frames divided by amino acids.}
  \item{n}{a list of \code{\link{n}} vectors divided by amino acids.}
  \item{seq.data}{a vector of \code{\link{seq.data}} format.}

  \item{x}{a codon or dna string, such "ACG", "acg", or "A", "a".}

  \item{b}{a \code{\link{b}} object.}
  \item{bVec}{a \code{\link{bVec}} object.}
  \item{aa.names}{a vector contains amino acid names for analysis.}
  \item{model}{model fitted.}
}
\details{
  \code{convert.reu13.df.to.list()}, \code{convert.y.to.list()}, and
  \code{convert.n.to.list()}:
  these utility functions take the inputs divided by amino acids
  and return the outputs divided by ORFs.

  \code{convert.y.scuo()} converts \code{\link{y}} into \code{\link{scuo}}
  format.

  \code{convert.seq.data.to.string()} converts \code{\link{seq.data}} into
  \code{\link{seq.string}} format.

  \code{codon.low2up()} and \code{codon.up2low()} convert codon strings
  between lower or upper cases.

  \code{convert.bVec.to.b()} and \code{convert.b.to.bVec()} convert
  objects \code{\link{b}} and \code{\link{bVec}}.
}
\value{
  All functions return the corresponding formats.
}
\references{
  \url{https://github.com/snoweye/cubfits/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}.
}
\seealso{
  \link{AllDataFormats},
  \code{\link{rearrange.n}()},
  \code{\link{rearrange.reu13.df}()},
  \code{\link{rearrange.y}()}, and
  \code{\link{read.seq}()}.
}
\examples{
\dontrun{
suppressMessages(library(cubfits, quietly = TRUE))

reu13.list <- convert.reu13.df.to.list(ex.train$reu13.df)
y.list <- convert.y.to.list(ex.train$y)
n.list <- convert.n.to.list(ex.train$n)

scuo <- convert.y.to.scuo(ex.train$y)

seq.data <- read.seq(get.expath("seq_200.fasta"))
seq.string <- convert.seq.data.to.string(seq.data)

codon.low2up("acg")
codon.up2low("ACG")

dna.low2up(c("a", "c", "g"))
dna.up2low(c("A", "C", "G"))
}
}
\keyword{utility}
