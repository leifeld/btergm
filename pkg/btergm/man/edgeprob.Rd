\name{edgeprob}
\alias{edgeprob}
\title{
Compute all dyadic edge probabilities for an ERGM or TERGM.
}
\description{
\code{edgeprob} is a convenience function that creates a data frame with all dyads in the ERGM or TERGM along with their edge probabilities and their predictor values (i.e., change statistics). This is useful for creating marginal effects plots or contrasting multiple groups of dyads. Parallel processing is available for this task. See also the \link{interpret} help page. The \link{interpret} function is used internally to compute the probabilities.
}
\usage{
edgeprob(object, parallel = c("no", "multicore", "snow"), 
    ncpus = 1, cl = NULL, verbose = TRUE)
}
\arguments{
\item{object}{ An \code{ergm}, \code{btergm}, or \code{mtergm} object. }
\item{parallel}{ Parallel processing. Valid options are "no" (for no parallel computing), "multicore" (for forking), and "snow" (for a PSOCK cluster). See the \link{btergm} and \link{gof-methods} help pages for details on these options. }
\item{ncpus}{ Number of cores to use for parallel processing. See the \link{btergm} and \link{gof-methods} help pages for details. }
\item{cl}{ Cluster object for parallel processing with the "snow" option (i.e., a PSOCK or MPI cluster). If \code{NULL}, a cluster object will be created internally. See the \link{btergm} and \link{gof-methods} help pages for details. }
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