% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mtergm.R
\docType{class}
\name{mtergm-class}
\alias{mtergm-class}
\alias{show,mtergm-method}
\alias{coef,mtergm-method}
\alias{nobs,mtergm-method}
\alias{timesteps.mtergm}
\alias{summary,mtergm-method}
\title{An S4 Class to represent a fitted TERGM by MCMC-MLE}
\usage{
\S4method{show}{mtergm}(object)

\S4method{coef}{mtergm}(object, invlogit = FALSE, ...)

\S4method{nobs}{mtergm}(object)

timesteps.mtergm(object)

\S4method{summary}{mtergm}(object, ...)
}
\arguments{
\item{object}{An \code{mtergm} object.}

\item{invlogit}{Apply inverse logit transformation to the estimates and/or
confidence intervals? That is, \eqn{\frac{1}{1 + \exp(-x)}}, where \eqn{x}
is the respective value.}

\item{...}{Currently not in use.}
}
\description{
An S4 class to represent a fitted TERGM by MCMC-MLE.
}
\details{
\code{mtergm} objects result from MCMC-MLE-based estimation of a TERGM via
the \code{\link{mtergm}} function. They contain the coefficients, standard
errors, and p-values, among other details.
}
\section{Functions}{
\itemize{
\item \code{show(mtergm)}: Show the coefficients of an \code{mtergm} object.

\item \code{coef(mtergm)}: Return the coefficients of an \code{mtergm} object.

\item \code{nobs(mtergm)}: Return the coefficients of an \code{mtergm} object.

\item \code{timesteps.mtergm()}: Return the number of time steps saved in an
\code{mtergm} object.

\item \code{summary(mtergm)}: Return the coefficients of an \code{mtergm} object.

}}
\section{Slots}{

\describe{
\item{\code{coef}}{Object of class \code{"numeric"}. The coefficients.}

\item{\code{se}}{Object of class \code{"numeric"}. The standard errors.}

\item{\code{pval}}{Object of class \code{"numeric"}. The p-values.}

\item{\code{nobs}}{Object of class \code{"numeric"}. Number of observations.}

\item{\code{time.steps}}{Object of class \code{"numeric"}. Number of time steps.}

\item{\code{formula}}{Object of class \code{"formula"}. The original model formula
(without indices for the time steps).}

\item{\code{formula2}}{The revised formula with the object references after applying
the \code{\link{tergmprepare}} function.}

\item{\code{auto.adjust}}{Object of class \code{"logical"}. Indicates whether
automatic adjustment of dimensions was done before estimation.}

\item{\code{offset}}{Object of class \code{"logical"}. Indicates whether an offset
matrix with structural zeros was used.}

\item{\code{directed}}{Object of class \code{"logical"}. Are the dependent networks
directed?}

\item{\code{bipartite}}{Object of class \code{"logical"}. Are the dependent networks
bipartite?}

\item{\code{estimate}}{Estimate: either MLE or MPLE.}

\item{\code{loglik}}{Log likelihood of the MLE.}

\item{\code{aic}}{Akaike's Information Criterion.}

\item{\code{bic}}{Bayesian Information Criterion.}

\item{\code{ergm}}{The original \code{ergm} object as estimated by the
\code{\link[ergm]{ergm}} function in the \pkg{ergm} package.}

\item{\code{nvertices}}{Number of vertices.}

\item{\code{data}}{The data after processing by the \code{\link{tergmprepare}}
function.}
}}

\seealso{
Other tergm-classes: 
\code{\link{btergm-class}},
\code{\link{createBtergm}()},
\code{\link{createMtergm}()},
\code{\link{createTbergm}()},
\code{\link{tbergm-class}}
}
\author{
Philip Leifeld
}
\concept{tergm-classes}
