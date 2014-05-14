\name{Data Formats}
\alias{b}
\alias{bVec}
\alias{n}
\alias{n.list}
\alias{phi.df}
\alias{reu13.df}
\alias{reu13.list}
\alias{scuo}
\alias{seq.string}
\alias{seq.data}
\alias{phi.Obs}
\alias{y}
\alias{y.list}
\alias{AllDataFormats}
\docType{data}
\title{Data Formats}
\description{
  Data formats used in \pkg{cubfits}.
}
\format{
  All are in simple formats as S3 default lists or data frames.
}
\details{
  \itemize{
  \item
    Format \code{b}: \cr
    A named list \code{A} contains amino acids.
    Each element of the list \code{A[[i]]} is a list of elements
    \code{coefficients} (coefficients of log(mu) and Delta.t),
    \code{coef.mat} (matrix format of \code{coefficients}), and
    \code{R} (covariance matrix of \code{coefficients}).
    Note that \code{coefficients} and \code{R} are typically as in the output
    of \code{vglm()} of \pkg{VGAM} package.
    Also, \code{coef.mat} and \code{R} may miss in some cases.
    \cr
    e.g. \code{A[[i]]$coef.mat} is the regression beta matrix of \code{i}-th
         amino acid.

  \item
    Format \code{bVec}: \cr
    A vector simply contains all coefficients of a \code{b} object \code{A}.
    Note that this is probably only used inside MCMC or the output of
    \code{vglm()} of \pkg{VGAM} package.
    \cr
    e.g. \code{do.call("c", lapply(A, function(x) x$coefficients))}.

  \item
    Format \code{n}: \cr
    A named list \code{A} contains amino acids.
    Each element of the list \code{A[[i]]} is a vector containing total
    codon counts.
    \cr
    e.g. \code{A[[i]][j]} is for \code{j}-th ORF of \code{i}-th amino acid
         \code{names(A)[i]}.

  \item
    Format \code{n.list}: \cr
    A named list \code{A} contains ORFs.
    Each element of the list \code{A[[i]]} is a named list of amino acid
    containing total count.
    \cr
    e.g. \code{A[[i]][[j]]} contains total count of
         \code{j}-th amino acid in \code{i}-th ORF.

  \item
    Format \code{phi.df}: \cr
    A data frame \code{A} contains two columns \code{ORF} and \code{phi.value}.
    \cr
    e.g. \code{A[i,]} is for \code{i}-th ORF.

  \item
    Format \code{reu13.df}: \cr
    A named list \code{A} contains amino acids.
    Each element is a data frame summarizing ORF and expression.
    The data frame has four to five columns including
    \code{ORF}, \code{phi} (expression), \code{Pos} (amino acid position),
    \code{Codon} (synonymous codon), and
    \code{Codon.id} (synonymous codon id, for computing only).
    Note that \code{Codon.id} may miss in some cases.
    \cr
    e.g. \code{A[[i]][17,]} is the 17-th recode of \code{i}-th amino acid.

  \item
    Format \code{reu13.list}: \cr
    A named list \code{A} contains ORFs.
    Each element is a named list \code{A[[i]]} contains amino acids.
    Each element of nested list \code{A[[i]][[j]]} is a position vector
    of synonymous codon.
    \cr
    e.g. \code{A[[i]][[j]][k]} is the \code{k}-th synonymous codon position of
         \code{j}-th amino acid in the \code{i}-th ORF.

  \item
    Format \code{scuo}: \cr
    A data frame of 8 named columns includes
    \code{AA} (amino acid), \code{ORF}, \code{C1}, \dots, \code{C6}
    where \code{C*}'s are for codon counts.

  \item
    Format \code{seq.string}: \cr
    Default outputs of \code{read.fasta()} of \pkg{seqinr} package.
    A named list \code{A} contains ORFs.
    Each element of the list is a long string of a ORF.
    \cr
    e.g. \code{A[[i]][1]} or \code{A[[i]]} is the sequence of
         \code{i}-th ORF.

  \item
    Format \code{seq.data}: \cr
    Converted from \code{seq.string} format.
    A named list \code{A} contains ORFs.
    Each element of the list \code{A[[i]]} is a string vector.
    Each element of the vector is a codon string.
    \cr
    e.g. \code{A[[i]][j]} is \code{i}-th ORF and \code{j}-th codon.

  \item
    Format \code{phi.Obs}: \cr
    A named vector \code{A} of observed expression values and possibly
    with measurement errors.
    \cr
    e.g. \code{A[i]} is the observed phi value of \code{i}-th ORF.

  \item
    Format \code{y}: \cr
    A named list \code{A} contains amino acids.
    Each element of the list \code{A[[i]]} is a matrix
    where ORFs are in row and synonymous codons are in column.
    The element of the matrix contains codon counts.
    \cr
    e.g. \code{A[[i]][j, k]} is the count for \code{i}-th amino acid,
         \code{j}-th ORF, and \code{k}-th synonymous codon.

  \item
    Format \code{y.list}: \cr
    A named list \code{A} contains ORFs.
    Each element of the list \code{A[[i]]} is a named list \code{A[[i]][[j]]}
    contains amino acids.
    The element of amino acids list is a codon count vector.
    \cr
    e.g. \code{A[[i]][[j]][k]} is the count for \code{i}-th ORF,
         \code{j}-th amino acid, and \code{k}-th synonymous codon.
  }
}
\references{
  \url{https://github.com/snoweye/cubfits/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}.
}
\keyword{dataformats}
