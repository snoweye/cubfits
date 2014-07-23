\name{Controls}
\alias{.cubfitsEnv}
\alias{.CF.CT}
\alias{.CF.CONF}
\alias{.CF.DP}
\alias{.CF.GV}
\alias{.CF.OP}
\alias{.CF.AC}
\alias{.CF.PT}
\alias{.CF.PARAM}
\alias{.CO.CT}
\docType{data}
\title{Default Controlling Options}
\description{
  Default controls of \pkg{cubfits} include for models, optimizations, MCMC,
  plotting, global variables, etc.
}
\usage{
  .cubfitsEnv
  .CF.CT
  .CF.CONF
  .CF.GV
  .CF.DP
  .CF.OP
  .CF.AC
  .CF.PT
  .CF.PARAM
  .CO.CT
}
\format{
  All are in lists and contain several controlling options.
}
\details{
  See \code{\link{init.function}()} for use cases of these objects.

  \itemize{
  \item
    \code{.cubfitEnv} is a default environment to dynamically save functions
    and objects.

  \item
    \code{.CF.CT} is main controls of models. It currently includes
    \tabular{ll}{
      \code{model}     \tab main models \cr
      \code{type.p}    \tab proposal for hyper-parameters \cr
      \code{type.Phi}  \tab proposal for Phi \cr
      \code{model.Phi} \tab prior of Phi \cr
      \code{init.Phi}  \tab initial methods for Phi \cr
      \code{init.fit}  \tab how is coefficient proposed \cr
      \code{parallel}  \tab parallel functions \cr
      \code{adaptive}  \tab method for adaptive MCMC
    }

  \item
    \code{.CF.CONF} controls the initial and draw scaling. 
    It currently includes
    \tabular{ll}{
      \code{scale.phi.Obs}           \tab if phi were scaled to mean 1 \cr
      \code{init.b.Scale}            \tab initial b scale \cr
      \code{init.phi.Scale}          \tab initial phi scale \cr
      \code{p.nclass}                \tab number of classes if mixture phi \cr
      \code{b.DrawScale}             \tab drawing scale for b if random walk \cr
      \code{p.DrawScale}             \tab drawing scale for p if random walk \cr
      \code{phi.DrawScale}           \tab random walk scale for phi \cr
      \code{phi.pred.DrawScale}      \tab random walk scale for phi.pred \cr
      \code{sigma.Phi.DrawScale}     \tab random walk scale for sigma.Phi \cr
      \code{bias.Phi.DrawScale}      \tab random walk scale for bias.Phi \cr
      \code{estimate.bias.Phi}       \tab if estimate bias of phi during MCMC \cr
      \code{compute.logL}            \tab if compute logL in each iteration
    }

  \item
    \code{.CF.GV} contains global variables for amino acids and codons.
    It currently includes
    \tabular{ll}{
      \code{amino.acid}             \tab amino acids \cr
      \code{amino.acid.3}           \tab amino acids \cr
      \code{synonymous.codon}       \tab synonymous codons of amino acids \cr
      \code{amino.acid.split}       \tab amino acid 'S' is split \cr
      \code{amino.acid.split.3}     \tab amino acid 'S' is split \cr
      \code{synonymous.codon.split} \tab synonymous codons of split amino acid
    }

  \item
    \code{.CF.OP} controls optimizations. It currently includes
    \tabular{ll}{
      \code{optim.method}    \tab method for \code{optim}() \cr
      \code{stable.min.exp}  \tab minimum exponent \cr
      \code{stable.max.exp}  \tab maximum exponent \cr
      \code{E.Phi}           \tab expected Phi \cr
      \code{lower.optim}     \tab lower of derivative of logL(x) \cr
      \code{upper.optim}     \tab upper of derivative of logL(x) \cr
      \code{lower.integrate} \tab lower of integration of L(x) \cr
      \code{upper.integrate} \tab upper of integration of L(x)
    }

  \item
    \code{.CF.DP} is for dumping MCMC iterations. It currently includes
    \tabular{ll}{
      \code{dump}             \tab if dumping within MCMC \cr
      \code{iter}             \tab iterations per dumping \cr
      \code{prefix.dump}      \tab path and file names of dumping \cr
      \code{verbose}          \tab if verbose \cr
      \code{iterThin}         \tab iterations to thin chain \cr
      \code{report}           \tab iterations to report \cr
      \code{report.proc}      \tab iterations to report \code{proc.time}()
    }

  \item
    \code{.CF.AC} controls adaptive MCMC. It currently includes
    \tabular{ll}{
      \code{renew.iter}          \tab per renewing iterations \cr
      \code{target.accept.lower} \tab target acceptant rate lower bound \cr
      \code{target.accept.upper} \tab target acceptant rate upper bound \cr
      \code{scale.increase}      \tab increase scale size \cr
      \code{scale.decrease}      \tab decrease scale size \cr
      \code{sigma.lower}         \tab lower bound of relative scale size \cr
      \code{sigma.upper}         \tab upper bound of relative scale size
    }

  \item
    \code{.CF.PT} controls the plotting format. It currently includes
    \tabular{ll}{
      \code{color} \tab color for codons.
    }

  \item
    \code{.CF.PARAM} controls the parameters and hyperparameters of priors.
    It currently includes
    \tabular{ll}{
      \code{phi.meanlog}      \tab mean of phi in loca scale \cr
      \code{phi.sdlog}        \tab standard deviation of phi in loca scale
    }

  \item
    \code{.CO.CT} controls the constrained optimization function. It currently
    includes
    \tabular{ll}{
      \code{debug} \tab message printing level of debugging.
    }
  }
}
\references{
  \url{https://github.com/snoweye/cubfits/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}.
}
\seealso{
  \code{\link{init.function}()}, \code{\link{cubfits}()},
  \code{\link{cubpred}()}, \code{\link{cubappr}()}, and
  \code{\link{mixnormerr.optim}()}.
}
\examples{
\dontrun{
suppressMessages(library(cubfits, quietly = TRUE))

.CF.CT
.CF.CONF
.CF.DP
.CF.GV
.CF.OP
.CF.AC
.CF.PT
.CF.PARAM
.CO.CT

ls(.cubfitsEnv)
init.function()
ls(.cubfitsEnv)
}
}
\keyword{datasets}
