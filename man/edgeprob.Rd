\name{edgeprob}
\alias{edgeprob}
\title{
Compute all dyadic edge probabilities for an ERGM or TERGM.
}
\description{
\code{edgeprob} is a convenience function that creates a data frame with
all dyads in the ERGM or TERGM along with their edge probabilities and
their predictor values (i.e., change statistics). This is useful for
creating marginal effects plots or contrasting multiple groups of
dyads. This function works faster than the \link{interpret} function.
See also the \link{interpret} help page.
}
\usage{
edgeprob(object, verbose = FALSE)
}
\arguments{
\item{object}{ An \code{ergm}, \code{btergm}, or \code{mtergm} object. }
\item{verbose}{ Print details? }
}
\value{
The first variable in the resulting data frame contains the edge value (i.e., the dependent variable, which is usually binary). The next variables contain all the predictors from the ERGM or TERGM (i.e., the change statistics). The next three variables contain the indices of the sender (i), the receiver (j), and the time step (t). These three indices serve to identify the dyad. The last variable contains the computed edge probabilities.
}
\author{
Philip Leifeld
}
\seealso{
\link{interpret} \link{btergm-package} \link{btergm}
}