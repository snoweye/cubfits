\name{Simulation Tool}
\alias{simu.orf}
\alias{simu.phi.Obs}
\alias{simu.mixnormerr}
\title{Simulate ORFs and Expression Data}
\description{
  These utility functions generate data for simulation studies
  including fake ORFs and expression values.
}
\usage{
  simu.orf(n, bInit, phi.Obs = NULL, AA.prob = NULL, orf.length = NULL,
           orf.names = NULL, model = .CF.CT$model)
  simu.phi.Obs(Phi, sigmaW.lim = 1)
  simu.mixnormerr(n, param)
}
\arguments{
  \item{n}{number of ORFs or sequences.}
  \item{bInit}{parameters of mutation and selection of format
               \code{\link{b}}.}
  \item{phi.Obs}{an object of format \code{\link{phi.Obs}}.}
  \item{AA.prob}{proportion of amino acids.}
  \item{orf.length}{lengths of ORFs.}
  \item{orf.names}{names of ORFs.}
  \item{model}{model to be simulated.}
  \item{Phi}{expression values (potentially true expression).}
  \item{sigmaW.lim}{std of measurement errors (between Phi and phi.Obs).}
  \item{param}{as in \code{\link{dmixnormerr}()}}
}
\details{
  \code{simu.orf()} generates ORFs or sequences based on the \code{bInit}
  and \code{phi.Obs}.

  If \code{phi.Obs} is omitted, then standard log normal random variables
  are instead).

  If \code{AA.prob} is omitted, then uniform proportion is assigned.

  If \code{orf.length} is omitted, then 10 to 20 codons are randomly
  assigned.

  If \code{orf.names} is omitted, then "ORF1" to "ORFn" are assigned.

  \code{simu.phi.Obs()} generates \code{\link{phi.Obs}} by adding normal random
  errors to \code{Phi}, and errors have mean 0 and standard deviation
  \code{sigmaW.lim}.

  \code{simu.mixnormerr()} generates \code{Phi} according to the \code{param},
  and adds normal random errors to \code{Phi}.
}
\value{
  \code{simu.orf()} returns a list of format \code{\link{seq.data}}.

  \code{simu.phi.Obs()} returns a vector of format \code{\link{phi.Obs}}.

  \code{simu.mixnormerr()} returns a list contains three vectors of length
  \code{n}: one for expected gene expression \code{Phi}, one for observed
  gene expression \code{phi.Obs}, and one for the component id \code{id.K}.
}
\references{
  \url{https://github.com/snoweye/cubfits/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}.
}
\seealso{
  \code{\link{read.seq}()}, \code{\link{read.phi.df}()},
  \code{\link{write.seq}()}, \code{\link{write.phi.df}()}, and
  \code{\link{mixnormerr.optim}()}.
}
\examples{
\dontrun{
suppressMessages(library(cubfits, quietly = TRUE))
set.seed(1234)

# Generate sequences.
da.roc <- simu.orf(length(ex.train$phi.Obs), bInit$roc,
                   phi.Obs = ex.train$phi.Obs, model = "roc")
names(da.roc) <- names(ex.train$phi.Obs)
write.fasta(da.roc, names(da.roc), "toy_roc.fasta")
}
}
\keyword{tool}
