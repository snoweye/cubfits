\name{Initial Generic Functions}
\alias{init.function}
\title{Initial Generic Functions of Codon Usage Bias Fits} 
\description{
  Initial generic functions for model fitting/approximation/prediction
  of \pkg{cubfits}.
}
\usage{
  init.function(model = .CF.CT$model[1],
                type.p = .CF.CT$type.p[1],
                type.Phi = .CF.CT$type.Phi[1],
                model.Phi = .CF.CT$model.Phi[1],
                init.Phi = .CF.CT$init.Phi[1],
                init.fit = .CF.CT$init.fit[1],
                parallel = .CF.CT$parallel[1],
                adaptive = .CF.CT$adaptive[1])
}
\arguments{
  \item{model}{main fitted model.}
  \item{type.p}{proposal method for hyper-parameters.}
  \item{type.Phi}{proposal method for Phi (true expression values).}
  \item{model.Phi}{prior of Phi.}
  \item{init.Phi}{initial methods for Phi.}
  \item{init.fit}{how is coefficient initialed in \code{vglm()} of \pkg{VGAM}.}
  \item{parallel}{parallel functions.}
  \item{adaptive}{method for adaptive MCMC.}
}
\details{
  This function mainly takes the options, find the according generic
  functions, and assign those functions to \code{\link{.cubfitsEnv}}.
  Those generic functions can be executed accordingly later within functions
  for MCMC or multinomial logistic regression such as \code{\link{cubfits}()},
  \code{\link{cubappr}()}, and \code{\link{cubpred}()}.
  By default, those options are provided by \code{\link{.CF.CT}} which also
  leaves rooms for extensions of more complicated models and further
  optimizations.

  It is supposed to call this function before running any MCMC or
  multinomial logistic regression. This function may affect
  \code{\link{cubfits}()}, \code{\link{cubpred}()}, \code{\link{cubappr}()},
  \code{\link{estimatePhi}()}, and \code{\link{fitMultinom}()}.

  \itemize{
    \item
    \code{model} is the main fitting model, currently only \code{roc} is
    fully supported.

    \item
    \code{type.p} is for proposing hyper-parameters in Gibb sampler. Currently,
    \code{lognormal_fix} is suggested where mean 1 is fixed for log normal
    distribution. Conjugated prior and flat prior exist and are easily available
    in this step

    \item
    \code{type.Phi} is for proposing Phi (expression values) in the random walk
    chain updates. Only, \code{RW_Norm} is supported. Usually, the acceptance
    ratio can be adapted within 25\% and 50\% controlled by
    \code{\link{.CF.AC}} if \code{adaptive = simple}.

    \item
    \code{model.Phi} is for the distribution of Phi. Typically, log normal
    distribution \code{lognormal} is assumed.

    \item
    \code{init.Phi} is a way to initial Phi. Posterior mean \code{PM}
    is recommended which avoid boundary values.

    \item
    \code{init.fit} is a way of initial coefficients to fit mutation and
    selection coefficients (\eqn{\log\mu}{log.mu} and
    \eqn{\Delta t}{Delta.t} or \eqn{\omega}{omega})
    in \code{vglm()}. Option \code{current} means the \code{\link{b}}
    (log(mu) and Delta.t) of current MCMC iteration is the initial values, while
    \code{random} means \code{vglm()} provides the initial values.

    \item
    \code{parallel} is a way of parallel methods to speed up code.
    \code{lapply} means \code{lapply()} is used and no parallel;
    \code{mclapply} means \code{mclapply()} of \pkg{parallel} is used and
    good for shared memory machines;
    \code{task.pull} means \code{task.pull()} of \pkg{pbdMPI} is used and
    good for heterogeneous machines;
    \code{pbdLapply} means \code{pbdLapply()} of \pkg{pbdMPI} is used and
    good for homogeneous machines.
    Among those, \code{task.pull} is tested thoroughly and is the most reliable
    and efficient method.

    \item
    \code{adaptive} is a way for adaptive MCMC that propose better mixing
    distributions for random walks of Phi. The \code{simple} method is
    suggested and only the proposal distribution of Phi
    (\code{type.Phi = RW_Norm}) is adjusted gradually.
  }
}
\value{
  Return an invisible object which is a list contain all
  generic functions according to the input options.
  All functions are also assigned in the \code{\link{.cubfitsEnv}}
  for later evaluations called by MCMC or multinomial logistic regression.
}
\note{
  Note that all options are taken default values from the global control
  object \code{\link{.CF.CT}}, so one can utilize/alter the object's values
  to adjust those affected functions.

  Note that \code{phi.Obs} should be scaled to mean 1 before
  applying to MCMC.
}
\references{
  \url{https://github.com/snoweye/cubfits/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}.
}
\seealso{
  \code{\link{.CF.CT}}, \code{\link{.CF.CT}}, \code{\link{cubfits}()},
  \code{\link{cubpred}()}, and \code{\link{cubappr}()}. 
}
\examples{
\dontrun{
suppressMessages(library(cubfits, quietly = TRUE))
set.seed(1234)

# Convert data.
reu13.list <- convert.reu13.df.to.list(ex.test$reu13.df)
y.list <- convert.y.to.list(ex.test$y)
n.list <- convert.n.to.list(ex.test$n)

# Get phi.pred.Init
init.function(model = "roc")
fitlist <- fitMultinom(ex.train$reu13.df, ex.train$phi.Obs, ex.train$y,
                       ex.train$n)
phi.pred.Init <- estimatePhi(fitlist, reu13.list, y.list, n.list,
                         E.Phi = median(ex.test$phi.Obs),
                         lower.optim = min(ex.test$phi.Obs) * 0.9,
                         upper.optim = max(ex.test$phi.Obs) * 1.1)
}
}
\keyword{tool}
