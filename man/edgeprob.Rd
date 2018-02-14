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
\references{
Czarna, Anna Z., Philip Leifeld, Magdalena Smieja, Michael Dufner 
and Peter Salovey (2016): Do Narcissism and Emotional Intelligence 
Win Us Friends? Modeling Dynamics of Peer Popularity Using 
Inferential Network Analysis. \emph{Personality and Social 
Psychology Bulletin} 42(11): 1588--1599.

Desmarais, Bruce A. and Skyler J. Cranmer (2012):
Micro-Level Interpretation of Exponential Random Graph Models with
Application to Estuary Networks.
\emph{The Policy Studies Journal} 40(3): 402--434.

Leifeld, Philip, Skyler J. Cranmer and Bruce A. Desmarais (2017): 
Temporal Exponential Random Graph Models with btergm: Estimation 
and Bootstrap Confidence Intervals. \emph{Journal of Statistical 
Software} 83(6): 1-36. \url{http://dx.doi.org/10.18637/jss.v083.i06}.
}
\seealso{
\link{interpret} \link{btergm-package} \link{btergm}
}
