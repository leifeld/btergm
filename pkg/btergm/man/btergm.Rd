\name{btergm}
\alias{btergm}
\alias{mtergm}
\alias{tergm}
\title{TERGM by bootstrapped pseudolikelihood or MCMC MLE}
\description{TERGM by bootstrapped pseudolikelihood or MCMC MLE.}
\usage{
btergm(formula, R = 500, offset = FALSE, parallel = c("no", 
    "multicore", "snow"), ncpus = 1, cl = NULL, 
    verbose = TRUE, ...)

mtergm(formula, offset = FALSE, constraints = ~ ., 
    estimate = c("MLE", "MPLE"), verbose = TRUE, ...)
}
\details{
The \code{btergm} function computes temporal exponential random graph models (TERGM) by bootstrapped pseudolikelihood, as described in Desmarais and Cranmer (2012).

The \code{mtergm} function computes TERGMs by MCMC MLE (or MPLE with uncorrected standard errors) via blockdiagonal matrices and structural zeros. The \code{btergm} function is faster than the \code{mtergm} function.
}
\arguments{
\item{formula}{ Formula for the TERGM. Model construction works like in the \pkg{ergm} package with the same model terms etc. (for a list of terms, see \code{help("\link[ergm]{ergm-terms}")}). The networks to be modeled on the left-hand side of the equation must be given either as a list of network objects with more recent networks last (i.e., chronological order) or as a list of matrices with more recent matrices at the end. \code{dyadcov} and \code{edgecov} terms accept time-independent covariates (as \code{network} or \code{matrix} objects) or time-varying covariates (as a list of networks or matrices with the same length as the list of networks to be modeled). }
\item{R}{ Number of bootstrap replications. The higher the number of replications, the more accurate but also the slower is the estimation. }
\item{offset}{ If \code{offset = TRUE} is set, a list of offset matrices (one for each time step) with structural zeros is handed over to the pseudolikelihood routine. The offset matrices contain structural zeros where either the dependent networks or any of the covariates have missing nodes (if \code{auto.adjust = TRUE} is used). All matrices and network objects are inflated to the dimensions of the largest object, and the offset matrices inform the estimation routine which dyads are constrained to be absent (that is, fixed at \code{-Inf}). If \code{offset = FALSE} is set (the default behavior) and \code{auto.adjust} is switched on, all nodes that are not present across all covariates and networks within a time step are removed completely from the respective object(s) before estimation begins. This is computationally more efficient than using an offset matrix. }
\item{parallel}{ Use multiple cores in a computer or nodes in a cluster to speed up bootstrapping computations. The default value \code{"no"} means parallel computing is switched off. If \code{"multicore"} is used, the \code{mclapply} function from the \pkg{parallel} package (formerly in the \pkg{multicore} package) is used for parallelization. This should run on any kind of system except MS Windows because it is based on forking. It is usually the fastest type of parallelization. If \code{"snow"} is used, the \code{parLapply} function from the \pkg{parallel} package (formerly in the \pkg{snow} package) is used for parallelization. This should run on any kind of system including cluster systems and including MS Windows. It is slightly slower than the former alternative if the same number of cores is used. However, \code{"snow"} provides support for MPI clusters with a large amount of cores, which \pkg{multicore} does not offer (see also the \code{cl} argument). The backend for the bootstrapping procedure is the \pkg{boot} package. }
\item{ncpus}{ The number of CPU cores used for parallel computing (only if \code{parallel} is activated). If the number of cores should be detected automatically on the machine where the code is executed, one can set \code{ncpus = detectCores()} after loading the \pkg{parallel} package. On some HPC clusters, the number of available cores is saved as an environment variable; for example, if MOAB is used, the number of available cores can sometimes be accessed using \code{Sys.getenv("MOAB_PROCCOUNT")}, depending on the implementation. }
\item{cl}{ An optional \pkg{parallel} or \pkg{snow} cluster for use if \code{parallel = "snow"}. If not supplied, a PSOCK cluster is created temporarily on the local machine. }
\item{constraints}{ Constraints of the ERGM. See \code{\link[ergm]{ergm}} for details. }
\item{estimate}{ Estimation procedure of the ERGM. MCMC MLE by default, but MPLE with uncorrected standard errors is possible. See \code{\link[ergm]{ergm}} for details. }
\item{verbose}{ Print details about data preprocessing and estimation settings. }
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
}
\seealso{
\link{btergm-package} \link{simulate.btergm} \link[btergm]{gof} \link{knecht} \link{btergm-class} \link{preprocess} \link{timecov}
}
\examples{
# A simple toy example:

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

summary(fit)               # show estimation results

\dontrun{
# The same example using MCMC MLE:

fit2 <- mtergm(networks ~ edges + istar(2) + 
    edgecov(covariates))

summary(fit2)
}

# For an example with real data, see help("knecht").


# Examples for parallel processing:

# Some preliminaries: 
# - "Forking" means running the code on multiple cores in the same 
#   computer. It's fast but consumes a lot of memory because all 
#   objects are copied for each node. It's also restricted to 
#   cores within a physical computer, i.e. no distribution over a 
#   network or cluster. Forking does not work on Windows systems.
# - "MPI" is a protocol for distributing computations over many 
#   cores, often across multiple physical computers/nodes. MPI 
#   is fast and can distribute the work across hundreds of nodes 
#   (but remember that R can handle a maximum of 128 connections, 
#   which includes file access and parallel connections). However, 
#   it requires that the Rmpi package is installed and that an MPI 
#   server is running (e.g., OpenMPI).
# - "PSOCK" is a TCP-based protocol. It can also distribute the 
#   work to many cores across nodes (like MPI). The advantage of 
#   PSOCK is that it can as well make use of multiple nodes within 
#   the same node or desktop computer (as with forking) but without 
#   consuming too much additional memory. However, the drawback is 
#   that it is not as fast as MPI or forking.
# The following code provides examples for these three scenarios.

# btergm works with clusters via the parallel package. That is, the 
# user can create a cluster object (of type "PSOCK", "MPI", or 
# "FORK") and supply it to the 'cl' argument of the 'btergm' 
# function. If no cluster object is provided, btergm will try to 
# create a temporary PSOCK cluster (if parallel = "snow") or it 
# will use forking (if parallel = "multicore").

\dontrun{

# To use a PSOCK cluster without providing an explicit cluster 
# object:
require("parallel")
fit <- btergm(networks ~ edges + istar(2) + edgecov(covariates), 
    R = 100, parallel = "snow", ncpus = 25)

# Equivalently, a PSOCK cluster can be provided as follows:
require("parallel")
cores <- 25
cl <- makeCluster(cores, type = "PSOCK")
fit <- btergm(networks ~ edges + istar(2) + edgecov(covariates), 
    R = 100, parallel = "snow", ncpus = cores, cl = cl)
stopCluster(cl)

# Forking (without supplying a cluster object) can be used as 
# follows.
require("parallel")
cores <- 25
fit <- btergm(networks ~ edges + istar(2) + edgecov(covariates), 
    R = 100, parallel = "multicore", ncpus = cores)
stopCluster(cl)

# Forking (by providing a cluster object) works as follows:
require("parallel")
cores <- 25
cl <- makeCluster(cores, type = "FORK")
fit <- btergm(networks ~ edges + istar(2) + edgecov(covariates), 
    R = 100, parallel = "snow", ncpus = cores, cl = cl)
stopCluster(cl)

# To use MPI, a cluster object MUST be created beforehand. In 
# this example, a MOAB HPC server is used. It stores the number of 
# available cores as a system option:
require("parallel")
cores <- as.numeric(Sys.getenv("MOAB_PROCCOUNT"))
cl <- makeCluster(cores, type = "MPI")
fit <- btergm(networks ~ edges + istar(2) + edgecov(covariates), 
    R = 100, parallel = "snow", ncpus = cores, cl = cl)
stopCluster(cl)

# In the following example, the Rmpi package is used to create a 
# cluster. This may not work on all systems; consult your local 
# support staff or the help files on your HPC server to find out how 
# to create a cluster object on your system.

# snow/Rmpi start-up
if (!is.loaded("mpi_initialize")) {
    library("Rmpi")
}
library(snow);

mpirank <- mpi.comm.rank (0)
if (mpirank == 0) {
   invisible(makeMPIcluster())
} else {
  sink (file="/dev/null")
  invisible(slaveLoop (makeMPImaster()))
  mpi.finalize()
  q()
}
# End snow/Rmpi start-up

cl <- getMPIcluster()

fit <- btergm(networks ~ edges + istar(2) + edgecov(covariates), 
    R = 100, parallel = "snow", ncpus = 25, cl = cl)

}

}
