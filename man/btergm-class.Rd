\name{btergm-class}
\Rdversion{1.1}
\docType{class}
\alias{btergm-class}
\alias{mtergm-class}
\alias{coef,btergm-method}
\alias{coef,mtergm-method}
\alias{confint,btergm-method}
\alias{nobs,btergm-method}
\alias{nobs,mtergm-method}
\alias{show,btergm-method}
\alias{show,mtergm-method}
\alias{summary,btergm-method}
\alias{summary,mtergm-method}
\alias{timesteps.btergm}
\alias{timesteps.mtergm}
\alias{btergm.se}

\title{Classes \code{"btergm"} and \code{"mtergm"}}
\description{
\code{btergm} objects result from the estimation of a bootstrapped TERGM via the \code{btergm} function in the \pkg{xergm} package. \code{btergm} objects contain the coefficients, the bootstrapping samples of the coefficients, the number of replications, the number of observations, the number of time steps, the original formula, and the response, effects and weights objects that were fed into the \code{glm} call for estimating the model. \code{mtergm} objects result from MCMC-MLE-based estimation of a TERGM via the \code{mtergm} function. They contain the coefficients, standard errors, and p-values, among other details.
}
\section{Slots}{
  \describe{
    \item{\code{coef}:}{Object of class \code{"numeric"}. The coefficients. }
    \item{\code{bootsamp}:}{Object of class \code{"matrix"}. The bootstrapping sample. }
    \item{\code{R}:}{Object of class \code{"numeric"}. Number of replications. }
    \item{\code{nobs}:}{Object of class \code{"numeric"}. Number of observations. }
    \item{\code{time.steps}:}{Object of class \code{"numeric"}. Number of time steps. }
    \item{\code{formula}:}{Object of class \code{"formula"}. The original model formula (without indices for the time steps). }
    \item{\code{response}:}{Object of class \code{"integer"}. The response variable. }
    \item{\code{effects}:}{Object of class \code{"data.frame"}. The effects that went into the \code{glm} call. }
    \item{\code{weights}:}{Object of class \code{"integer"}. The weights of the observations. }
    \item{\code{auto.adjust}:}{Object of class \code{"logical"}. Indicates whether automatic adjustment of dimensions was done before estimation. }
    \item{\code{offset}:}{Object of class \code{"logical"}. Indicates whether an offset matrix with structural zeros was used. }
    \item{\code{directed}:}{Object of class \code{"logical"}. Are the dependent networks directed? }
    \item{\code{bipartite}:}{Object of class \code{"logical"}. Are the dependent networks bipartite? }
    \item{se}{Standard errors.}
    \item{pval}{p-values.}
    \item{estimate}{Estimate: either MCMC MLE or MPLE.}
    \item{loglik}{The log likelihood.}
    \item{aic}{Akaike's Information Criterion (AIC).}
    \item{bic}{The Bayesian Information Criterion (BIC).}
  }
}
\usage{
\S4method{summary}{btergm}(object, level = 0.95, type = "perc", 
    invlogit = FALSE, ...)

\S4method{summary}{mtergm}(object, ...)

\S4method{show}{btergm}(object)

\S4method{show}{mtergm}(object)

\S4method{nobs}{btergm}(object)

\S4method{nobs}{mtergm}(object)

\S4method{coef}{btergm}(object, invlogit = FALSE, ...)

\S4method{coef}{mtergm}(object, invlogit = FALSE, ...)

\S4method{confint}{btergm}(object, parm, level = 0.95, type = "perc", 
    invlogit = FALSE, ...)

btergm.se(object, print = FALSE)

timesteps.btergm(object)

timesteps.mtergm(object)
}
\details{
Various generic methods are available for \code{btergm} objects: The \code{coef} and \code{show} methods return the coefficients; the \code{summary} method gives a model summary. The \code{nobs} method returns the number of observations. The {confint} method returns confidence intervals from the bootstrap replications of \code{btergm} objects, and the user can specify the confidence level. The method returns a matrix with three columns: the estimate, the lower bound, and the upper bound of the confidence interval for each model term.

The \code{btergm.se} function computes standard errors and p values for \code{btergm} objects. It returns a matrix with four columns: the estimate, the standard error, the z value, and the p value for each model term. If the argument \code{print = TRUE} is used, the matrix is printed to the R console as a formatted coefficient matrix with significance stars instead. Note that confidence intervals are the preferred way of interpretation for bootstrapped TERGMs; standard errors are only accurate if the bootstrapped data are normally distributed, which is not always the case. Various methods for checking for normality for each model term are available, for example quantile-quantile plots (e.g., \code{qqnorm(x@boot$t[, 1])} for the first model term in the \code{btergm} object called \code{x}).

The \code{timesteps.btergm} function extracts the number of time steps from a \code{btergm} object. The number of time steps is the number of networks being modeled on the left-hand side of the model formula.

Some of these functions or methods are also available for \code{mtergm} objects.
}
\arguments{
\item{object}{ A \code{btergm} object. }
\item{level}{ The significance level for computation of the confidence intervals. The default is \code{0.95} (that is, an alpha value of 0.05). Other common values include \code{0.999}, \code{0.99}, \code{0.9}, and \code{0.5}. }
\item{parm}{ Parameters (specified by integer position or character string). }
\item{type}{ Type of confidence interval, e.g., basic bootstrap interval (\code{type = "basic"}), percentile-based interval (\code{type = "perc"}, which is the default option), or bias-adjusted and accelerated confidence interval (\code{type = "bca"}). All options from the \code{type} argument of the \link[boot]{boot.ci} function in the boot package can be used to generate confidence intervals. }
\item{invlogit}{ Apply inverse logit transformation to the estimates and/or confidence intervals? That is, \code{1 / (1 + exp(-x))}, where \code{x} is the respective value. }
\item{print}{ Should the formatted coefficient table be printed to the R console along with significance stars (\code{print = TRUE}), or should the plain coefficient matrix be returned (\code{print = FALSE})? }
\item{...}{ Further arguments to be handed over to subroutines. }
}
\references{
Cranmer, Skyler J., Tobias Heinrich and Bruce A. Desmarais (2014):
Reciprocity and the Structural Determinants of the International Sanctions
Network. \emph{Social Networks} 36(1): 5--22.
\url{http://dx.doi.org/10.1016/j.socnet.2013.01.001}.

Desmarais, Bruce A. and Skyler J. Cranmer (2012):
Statistical Mechanics of Networks: Estimation and Uncertainty.
\emph{Physica A} 391: 1865--1876.
\url{http://dx.doi.org/10.1016/j.physa.2011.10.018}.

Desmarais, Bruce A. and Skyler J. Cranmer (2010):
Consistent Confidence Intervals for Maximum Pseudolikelihood Estimators.
\emph{Neural Information Processing Systems 2010 Workshop on Computational 
Social Science and the Wisdom of Crowds}.

Leifeld, Philip, Skyler J. Cranmer and Bruce A. Desmarais (2017): 
Temporal Exponential Random Graph Models with btergm: Estimation 
and Bootstrap Confidence Intervals. \emph{Journal of Statistical 
Software}.
}
\seealso{
\link{btergm-package} \link{btergm} \link{simulate.btergm} \link[btergm]{gof} \link{knecht} \link{getformula} \link{interpret} \link{mtergm}
}
\keyword{classes}
