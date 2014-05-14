\name{Generating Utility}
\alias{gen.reu13.df}
\alias{gen.reu13.list}
\alias{gen.y}
\alias{gen.n}
\alias{gen.phi.Obs}
\alias{gen.scuo}
\alias{DataGenerating}
\title{Generating Data Structure}
\description{
  These utility functions generate and summarize sequence strings into several
  useful formats such as \code{\link{reu13.df}}, \code{\link{y}}, and
  \code{\link{n}}, etc.
}
\usage{
  gen.reu13.df(seq.string, phi.df = NULL, aa.names = .CF.GV$amino.acid,
               split.S = TRUE, drop.X = TRUE, drop.MW = TRUE,
               drop.1st.codon = TRUE)
  gen.y(seq.string, aa.names = .CF.GV$amino.acid,
        split.S = TRUE, drop.X = TRUE, drop.MW = TRUE)
  gen.n(seq.string, aa.names = .CF.GV$amino.acid,
        split.S = TRUE, drop.X = TRUE, drop.MW = TRUE)

  gen.reu13.list(seq.string, aa.names = .CF.GV$amino.acid,
                 split.S = TRUE, drop.X = TRUE, drop.MW = TRUE,
                 drop.1st.codon = TRUE)
  gen.phi.Obs(phi.df)
  gen.scuo(seq.string, aa.names = .CF.GV$amino.acid,
           split.S = TRUE, drop.X = TRUE, drop.MW = TRUE)
}
\arguments{
  \item{seq.string}{a list of sequence strings.}
  \item{phi.df}{a \code{\link{phi.df}} object returned from \code{\link{read.phi.df}()}.}
  \item{aa.names}{a vector contains amino acid names for analysis.}
  \item{split.S}{split amino acid 'S' if any.}
  \item{drop.X}{drop amino acid 'X' if any.}
  \item{drop.MW}{drop amino acid 'M' and 'W' if any.}
  \item{drop.1st.codon}{if drop the first codon.}
}
\details{
  These functions mainly take inputs of sequence strings
  \code{\link{seq.string}} or \code{\link{phi.df}} and turn them
  into corresponding format.
}
\value{
  The outputs are data structure in corresponding formats.
  See \link{AllDataFormats} for details.
}
\references{
  \url{https://github.com/snoweye/cubfits/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}.
}
\seealso{
  \link{AllDataFormats},
  \code{\link{read.seq}()}, \code{\link{read.phi.df}()}, and
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
reu13.df <- gen.reu13.df(seq.string, phi.df, aa.names)
reu13.list.new <- gen.reu13.list(seq.string, aa.names)
y <- gen.y(seq.string, aa.names)
n <- gen.n(seq.string, aa.names)
scuo <- gen.scuo(seq.string, aa.names)

# Convert to list format.
reu13.list <- convert.reu13.df.to.list(reu13.df)
y.list <- convert.y.to.list(y)
n.list <- convert.n.to.list(n)
}
}
\keyword{utility}
