\name{checkdegeneracy}
\docType{methods}
\alias{checkdegeneracy}
\alias{check.degeneracy}
\alias{degeneracy}
\alias{print.degeneracy}
\alias{checkdegeneracy,btergm-method}
\alias{checkdegeneracy,mtergm-method}
\alias{checkdegeneracy.btergm}
\alias{checkdegeneracy.mtergm}
\title{Degeneracy check for btergm and mtergm objects}
\description{
Assess degeneracy of btergm and mtergm models.
}
\details{
The methods for the generic \code{degeneracy} function implement a degeneracy check for btergm and mtergm objects. For btergm, this works by comparing the global statistics of simulated networks to those of the observed networks at each observed time step. If the global statistics differ significantly, this is indicated by small p-values. If there are many significant results, this indicates degeneracy. For mtergm, the \code{mcmc.diagnostics} function from the \pkg{ergm} package is used.}
\usage{
\S4method{checkdegeneracy}{mtergm}(object, ...)

\S4method{checkdegeneracy}{btergm}(object, nsim = 1000, 
    MCMC.interval = 1000, MCMC.burnin = 10000, verbose = FALSE)

\method{print}{degeneracy}(x, ...)
}
\arguments{
\item{object}{ A \code{btergm} or \code{mtergm} object, as estimated using the \code{btergm} or \code{mtergm} function. }
\item{nsim}{ The number of networks to be simulated at each time step. This number should be sufficiently large for a meaningful comparison. If possible, much more than 1,000 simulations. }
\item{MCMC.burnin}{ Internally, this package uses the simulation facilities of the \pkg{ergm} package to create new networks against which to compare the original network(s) for goodness-of-fit assessment. This argument sets the MCMC burnin to be passed over to the simulation command. The default value is \code{10000}. There is no general rule of thumb on the selection of this parameter, but if the results look suspicious (e.g., when the model fit is perfect), increasing this value may be helpful. }
\item{MCMC.interval}{ Internally, this package uses the simulation facilities of the \pkg{ergm} package to create new networks against which to compare the original network(s) for goodness-of-fit assessment. This argument sets the MCMC interval to be passed over to the simulation command. The default value is \code{1000}, which means that every 1000th simulation outcome from the MCMC sequence is used. There is no general rule of thumb on the selection of this parameter, but if the results look suspicious (e.g., when the model fit is perfect), increasing this value may be helpful. }
\item{verbose}{ Print details? }
\item{x}{ A \code{degeneracy} object created by the \code{checkdegeneracy} function. }
\item{...}{ Arbitrary further arguments. }
}
\seealso{
\link{btergm-package} \link{btergm} \link{gof}
}
\author{
Philip Leifeld (\url{http://www.philipleifeld.com})
}
