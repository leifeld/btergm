\name{simulate.btergm}
\alias{simulate.btergm}
\title{Simulate new networks from btergm objects}
\description{Simulate new networks from btergm objects.}
\usage{
\method{simulate}{btergm}(object, nsim = 1, seed = NULL, 
    index = NULL, formula = getformula(object), 
    coef = object@coef, verbose = TRUE, ...)
}
\details{
The \code{simulate.btergm} function is a wrapper for the \code{simulate.formula} function in the \pkg{ergm} package (see \code{help("\link[ergm]{simulate.formula}")}). It can be used to simulate new networks from a \code{btergm} object. The \code{index} argument specifies from which of the original networks the new network(s) should be simulated. For example, if \code{object} is an estimation based on cosponsorship networks from the 99th to the 107th Congress (as in Desmarais and Cranmer 2012), and the cosponsorship network in the 108th Congress should be predicted using the \code{simulate.btergm} function, then the argument \code{index = 9} should be passed to the function because the network should be based on the 9th network in the list (that is, the latest network, which is the cosponsorship network for the 107th Congress). Note that all relevant objects (the networks and the covariates) must be present in the workspace (as was the case during the estimation of the model).
}
\arguments{
\item{object}{ A \code{btergm} object, resulting from a call of the \code{\link{btergm}} function. }
\item{nsim}{ The number of networks to be simulated. Note that for values greater than one, a \code{network.list} object is returned, which can be indexed just like a \code{list} object, for example \code{mynetworks[[1]]} for the first simulated network in the object \code{mynetworks}. }
\item{seed}{ Random number integer seed. See \link[base]{set.seed}. }
\item{formula}{ A model formula from which the new network(s) should be simulated. By default, the formula is taken from the btergm \code{object}. }
\item{index}{ Index of the network from which the new network(s) should be simulated. The index refers to the list of response networks on the left-hand side of the model formula. Note that more recent networks are located at the end of the list. By default, the last (= most recent) network is used. }
\item{coef}{ A vector of parameter estimates. By default, the coefficients are extracted from the given btergm \code{object}. }
\item{verbose}{ Print additional details while running the simulations? }
\item{...}{ Arbitrary further arguments are handed over to the \code{\link[ergm]{simulate.formula}} function. For details, refer to the help page of the \code{\link[ergm]{simulate.formula}} function. }
}
\references{
Desmarais, Bruce A. and Skyler J. Cranmer (2012):
 Statistical Mechanics of Networks: Estimation and Uncertainty.
 \emph{Physica A} 391: 1865--1876.

Leifeld, Philip, Skyler J. Cranmer and Bruce A. Desmarais (2017): 
Temporal Exponential Random Graph Models with btergm: Estimation 
and Bootstrap Confidence Intervals. \emph{Journal of Statistical 
Software} 83(6): 1-36. \url{http://dx.doi.org/10.18637/jss.v083.i06}.
}
\seealso{
\link{btergm-package} \link{btergm} \link[btergm]{gof}
}
\examples{
# fit a TERGM to some toy data
library("statnet")
set.seed(5)
networks <- list()
for(i in 1:10){            # create 10 random networks with 10 actors
  mat <- matrix(rbinom(100, 1, .25), nrow = 10, ncol = 10)
  diag(mat) <- 0           # loops are excluded
  nw <- network(mat)       # create network object
  networks[[i]] <- nw      # add network to the list
}
covariates <- list()
for (i in 1:10) {          # create 10 matrices as covariate
  mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
  covariates[[i]] <- mat   # add matrix to the list
}
fit <- btergm(networks ~ edges + istar(2) +
    edgecov(covariates), R = 100)

# simulate 12 new networks from the last (= 10th) time step
sim1 <- simulate(fit, nsim = 12)

# simulate 1 network from the first time step
sim2 <- simulate(fit, index = 1)

# simulate network from t = 5 with larger covariate coefficient
coefs <- coef(fit)
coefs["edgecov.covariates[[i]]"] <- 0.5
sim3 <- simulate(fit, index = 5, coef = coefs)
}
