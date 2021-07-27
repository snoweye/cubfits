\name{Yassour2009}
\alias{yassour}
\docType{data}
\title{Yassour 2009 Yeast Experiment Dataset}
\description{
  Experiments and data are obtained from Yassour et. al. (2009).
}
\usage{
  yassour
}
\format{
  A \code{data.frame} contains 6303 rows and 5 columns:
  \code{ORF} is for gene names in character, and
  \code{YPD0.1}, \code{YPD0.2}, \code{YPD15.1}, and \code{YPD15.2} are
  gene expressions in positive double corresponding to 4 controlled
  Yeast experiments.
}
\details{
  The original data are available as the URL of the section of Source next.
  As the section of Examples next, data are selected from \code{SD3.xls} and
  reordered by \code{ORF}.

  For further analysis, the Examples section also provides how to
  convert them to \code{phi.Obs} values either in geometric means or individually.
}
\source{
  \url{https://www.pnas.org/content/early/2009/02/10/0812841106}

  \url{https://www.pnas.org/highwire/filestream/598612/field_highwire_adjunct_files/3/SD3.xls}

  Yassour M, Kaplan T, Fraser HB, Levin JZ, Pfiffner J, Adiconis X,
  Schroth G, Luo S, Khrebtukova I, Gnirke A, Nusbaum C, Thompson DA,
  Friedman N, Regev A. (2009) ``Ab initio construction of a eukaryotic
  transcriptome by massively parallel mRNA sequencing.''
  Proc Natl Acad Sci USA 106(9):3264-9. [PMID:19208812]
}
\references{
  Wallace E.W.J., Airoldi E.M., and Drummond D.A.
    ``Estimating Selection on Synonymous Codon Usage from Noisy Experimental
      Data''
    Mol Biol Evol (2013) 30(6):1438--1453.
}
\examples{
\dontrun{
### SD3.xls is available from the URL provided in the References.
da <- read.table("SD3.xls", header = TRUE, sep = "\t", quote = "",
                 stringsAsFactors = FALSE)

### Select ORF, YPD0.1, YPD0.2, YPD15.1, YPD15.2.
da <- da[, c(1, 8, 9, 10, 11)]
colnames(da) <- c("ORF", "YPD0.1", "YPD0.2", "YPD15.1", "YPD15.2")

### Drop inappropriate values (NaN, NA, Inf, -Inf, and 0).
tmp <- da[, 2:5]
id.tmp <- rowSums(is.finite(as.matrix(tmp)) & tmp != 0) >= 3
tmp <- da[id.tmp, 1:5]
yassour <- tmp[order(tmp$ORF),]    # cubfits::yassour

### Get geometric mean of phi.Obs and scaling similar to Wallace (2013).
GM <- apply(yassour[, -1], 1, function(x) exp(mean(log(x[x != 0]))))
phi.Obs <- GM / sum(GM) * 15000

### Get individual of phi.Obs.
GM <- apply(yassour[, -1], 1, function(x) exp(mean(log(x[x != 0]))))
phi.Obs.all <- yassour[, -1] / sum(GM) * 15000
phi.Obs.all[phi.Obs.all == 0] <- NA
}
}
\keyword{datasets}
