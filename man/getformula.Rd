\name{getformula}
\docType{methods}
\alias{getformula}
\alias{getformula-methods}
\alias{getformula,ergm-method}
\alias{getformula,btergm-method}
\alias{getformula,mtergm-method}
\title{Extract the formula from a model.}
\description{Extract the model formula from ergm or btergm objects.}
\usage{
\S4method{getformula}{ergm}(x)

\S4method{getformula}{btergm}(x)

\S4method{getformula}{mtergm}(x)
}
\details{
Extract the model formula from ergm or btergm objects.
}
\arguments{
\item{x}{ A model object, for example a \code{btergm} or an \code{ergm} object. }
}
\seealso{
\link[btergm]{gofmethods}
}
