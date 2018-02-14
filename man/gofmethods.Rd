\name{gof-methods}
\alias{gof-methods}
\alias{gofmethods}
\alias{gof.methods}
\alias{gof}
\alias{gof,btergm-method}
\alias{gof,mtergm-method}
\alias{gof,ergm-method}
\alias{gof,network-method}
\alias{gof,matrix-method}
\alias{gof,sienaFit-method}
\alias{gof.btergm}
\alias{gof.mtergm}
\alias{gof.sienaFit}
\alias{gof.network}
\alias{gof.matrix}
\alias{gof.ergm}
\alias{createGOF}
\docType{methods}
\title{Conduct Goodness-of-Fit Diagnostics on ERGMs, TERGMs, SAOMs, and logit models}
\description{
Assess goodness of fit of btergm and other network models.
}
\details{
The generic \code{gof} function provides goodness-of-fit measures and degeneracy checks for \code{btergm}, \code{mtergm}, \code{ergm}, \code{sienaFit}, and custom dyadic-independent models. The user can provide a list of network statistics for comparing simulated networks based on the estimated model with the observed network(s). See \link{gof-statistics}. The objects created by these methods can be displayed using various plot and print methods (see \link{gof-plot}).

In-sample GOF assessment is the default, which means that the same time steps are used for creating simulations and for comparison with the observed network(s). It is possible to do out-of-sample prediction by specifying a (list of) target network(s) using the \code{target} argument. If a formula is provided, the simulations are based on the networks and covariates specified in the formula. This is helpful in situations where complex out-of-sample predictions have to be evaluated. A usage scenario could be to simulate from a network at time \code{t} (provided through the \code{formula} argument) and compare to an observed network at time \code{t + 1} (the \code{target} argument). This can be done, for example, to assess predictive performance between time steps of the original networks, or to check whether the model performs well with regard to a newly measured network given the old data from the previous time step.

Predictive fit can also be assessed for stochastic actor-oriented models (SAOM) as implemented in the \pkg{RSiena} package. After compiling the usual objects (model, data, effects), one of the time steps can be predicted based on the previous time step and the SAOM using the \code{sienaFit} method of the \code{gof} function. By default, however, within-sample fit is used for SAOMs, just like for (T)ERGMs.

The \code{gof} methods for networks and matrices serve to assess the goodness of fit of a dyadic-independence model. To do this, the method requires a vector of coefficients (one coefficient for the intercept or \code{edges} term and one coefficient for each covariate), a list of covariates (in matrix or network shape), and a dependent network or matrix. This is useful for assessing the goodness of fit of QAP-adjusted logistic regression models (as implemented in the \code{netlogit} function in the \pkg{sna} package) or other dyadic-independence models, such as models fitted using \code{glm}. Note that this method only works with cross-sectional models and does not accept lists of networks as input data.

The \code{createGOF} function is used internally by the \code{gof} function in order to create a \code{gof} object from a list of simulated networks and a list of target networks to compare against. It can also be used directly by the end user if the user wants to supply lists of simulated and target networks from other sources.
}
\usage{
\S4method{gof}{btergm}(object, target = NULL, formula = getformula(object), 
    nsim = 100, MCMC.interval = 1000, MCMC.burnin = 10000, 
    parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL, 
    statistics = c(dsp, esp, deg, ideg, geodesic, rocpr, 
    walktrap.modularity), verbose = TRUE, ...)

\S4method{gof}{mtergm}(object, target = NULL, formula = getformula(object), 
    nsim = 100, MCMC.interval = 1000, MCMC.burnin = 10000, 
    parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL, 
    statistics = c(dsp, esp, deg, ideg, geodesic, rocpr, 
    walktrap.modularity), verbose = TRUE, ...)

\S4method{gof}{ergm}(object, target = NULL, formula = getformula(object), 
    nsim = 100, MCMC.interval = 1000, MCMC.burnin = 10000, 
    parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL, 
    statistics = c(dsp, esp, deg, ideg, geodesic, rocpr, 
    walktrap.modularity), verbose = TRUE, ...)

\S4method{gof}{matrix}(object, covariates, coef, target = NULL, nsim = 100, 
    mcmc = FALSE, MCMC.interval = 1000, MCMC.burnin = 10000, 
    parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL, 
    statistics = c(dsp, esp, deg, ideg, geodesic, rocpr, 
    walktrap.modularity), verbose = TRUE, ...)

\S4method{gof}{network}(object, covariates, coef, target = NULL, 
    nsim = 100, mcmc = FALSE, MCMC.interval = 1000, 
    MCMC.burnin = 10000, parallel = c("no", "multicore", "snow"), 
    ncpus = 1, cl = NULL, statistics = c(dsp, esp, deg, ideg, 
    geodesic, rocpr, walktrap.modularity), verbose = TRUE, ...)

\S4method{gof}{sienaFit}(object, period = NULL, parallel = c("no", 
    "multicore", "snow"), ncpus = 1, cl = NULL, structzero = 10, 
    statistics = c(esp, deg, ideg, geodesic, rocpr, walktrap.modularity), 
    groupName = object$f$groupNames[[1]], varName = NULL, 
    outofsample = FALSE, sienaData = NULL, sienaEffects = NULL, 
    nsim = NULL, verbose = TRUE, ...)

createGOF(simulations, target, statistics = c(dsp, esp, deg, 
    ideg, geodesic, rocpr, walktrap.modularity), parallel = "no", 
    ncpus = 1, cl = NULL, verbose = TRUE, ...)
}
\arguments{
\item{cl}{ An optional \pkg{parallel} or \pkg{snow} cluster for use if \code{parallel = "snow"}. If not supplied, a cluster on the local machine is created temporarily. }
\item{coef}{ A vector of coefficients. }
\item{covariates}{ A list of matrices or network objects that serve as covariates for the dependent network. The covariates in this list are automatically added to the formula as \code{edgecov} terms. }
\item{formula}{ A model formula from which networks are simulated for comparison. By default, the formula from the \code{btergm} object \code{x} is used. It is possible to hand over a formula with only a single response network and/or dyad or edge covariates or with lists of response networks and/or covariates. It is also possible to use indices like \code{networks[[4]]} or \code{networks[3:5]} inside the formula. }
\item{groupName}{ The group name used in the Siena model. }
\item{mcmc}{ Should statnet's MCMC methods be used for simulating new networks? If \code{mcmc = FALSE}, new networks are simulated based on predicted tie probabilities of the regression equation. }
\item{MCMC.burnin}{ Internally, this package uses the simulation facilities of the \pkg{ergm} package to create new networks against which to compare the original network(s) for goodness-of-fit assessment. This argument sets the MCMC burnin to be passed over to the simulation command. The default value is \code{10000}. There is no general rule of thumb on the selection of this parameter, but if the results look suspicious (e.g., when the model fit is perfect), increasing this value may be helpful. }
\item{MCMC.interval}{ Internally, this package uses the simulation facilities of the \pkg{ergm} package to create new networks against which to compare the original network(s) for goodness-of-fit assessment. This argument sets the MCMC interval to be passed over to the simulation command. The default value is \code{1000}, which means that every 1000th simulation outcome from the MCMC sequence is used. There is no general rule of thumb on the selection of this parameter, but if the results look suspicious (e.g., when the model fit is perfect), increasing this value may be helpful. }
\item{ncpus}{ The number of CPU cores used for parallel GOF assessment (only if \code{parallel} is activated). If the number of cores should be detected automatically on the machine where the code is executed, one can try the \code{detectCores()} function from the \pkg{parallel} package. On some HPC clusters, the number of available cores is saved as an environment variable; for example, if MOAB is used, the number of available cores can sometimes be accessed using \code{Sys.getenv("MOAB_PROCCOUNT")}, depending on the implementation. Note that the maximum number of connections in a single R session (i.e., to other cores or for opening files etc.) is 128, so fewer than 128 cores should be used at a time. }
\item{nsim}{ The number of networks to be simulated at each time step. Example: If there are six time steps in the \code{formula} and \code{nsim = 100}, a total of 600 new networks is simulated. The comparison between simulated and observed networks is only done within time steps. For example, the first 100 simulations are compared with the first observed network, simulations 101-200 with the second observed network etc. }
\item{object}{ A \code{btergm}, \code{ergm}, or \code{sienaFit} object (for the \code{btergm}, \code{ergm}, and \code{sienaFit} methods, respectively). Or a network object or matrix (for the \code{network} and \code{matrix} methods, respectively). }
\item{outofsample}{ Should out-of-sample prediction be attempted? If so, some additional arguments must be provided: \code{sienaData}, \code{sienaEffects}, and \code{nsim}. The \code{sienaData} object must contain a base and a target network for out-of-sample prediction. The \code{sienaEffects} must contain the effects to be used for the simulations. The estimates will be taken from the estimated \code{object}, and they will be injected into a new SAOM and fixed during the sampling procedure. \code{nsim} determines how many simulations are used for the out-of-sample comparison. }
\item{parallel}{ Use multiple cores in a computer or nodes in a cluster to speed up the simulations. The default value \code{"no"} means parallel computing is switched off. If \code{"multicore"} is used (only available for \code{sienaAlgorithm} and \code{sienaModel} objects), the \code{mclapply} function from the \pkg{parallel} package (formerly in the \pkg{multicore} package) is used for parallelization. This should run on any kind of system except MS Windows because it is based on forking. It is usually the fastest type of parallelization. If \code{"snow"} is used, the \code{parLapply} function from the \pkg{parallel} package (formerly in the \pkg{snow} package) is used for parallelization. This should run on any kind of system including cluster systems and including MS Windows. It is slightly slower than the former alternative if the same number of cores is used. However, \code{"snow"} provides support for MPI clusters with a large amount of cores, which \pkg{multicore} does not offer (see also the \code{cl} argument). Note that \code{"multicore"} will only work if all cores are on the same node. For example, if there are three nodes with eight cores each, a maximum of eight CPUs can be used. Parallel computing is described in more detail on the help page of \link{btergm}. }
\item{period}{ Which transition between time periods should be used for GOF assessment? By default, all transitions between all time periods are used. For example, if there are three consecutive networks, this will extract simulations from the transitions between 1 and 2 and between 2 and 3, respectively, and these simulations will be compared to the networks at time steps 2 and 3, respectively. The time period can be provided as a numeric, e.g., \code{period = 4} for extracting the simulations between time steps 4 and 5 (= the fourth transition) and predicting the fifth network. Values lower than 1 or larger than the number of consecutive networks minus 1 are therefore not permitted. This argument is only used if out-of-sample prediction is switched off. }
\item{sienaData}{ An object of the class \code{siena}, which is usually created using the \code{sienaDataCreate} function in the \pkg{RSiena} package. This argument is only used for out-of-sample prediction. The object must be based on a \code{sienaDependent} object that contains two networks: the base network from which to simulate forward, and the target network which you want to predict out-of-sample. The object can contain further objects for storing covariates etc. that are necessary for estimating new networks. The best practice is to create an object that is identical to the \code{siena} object used for estimating the model, except that it contains the base and the target network instead of the dependent variable/networks. }
\item{sienaEffects}{ An object of the class \code{sienaEffects}, which is usually created using the \code{getEffects()} and the \code{includeEffects()} functions in the \code{RSiena} package. The best practice is to provide a \code{sienaEffects} object that is identical to the object used to create the original model (that is, it should contain the same effects), except that it should be based on the \code{siena} object provided through the \code{sienaData} argument. In other words, the \code{sienaEffects} object should be based on the base and target network used for out-of-sample prediction, and it should contain the same effects as those used for the original estimation. This argument is used only for out-of-sample prediction. }
\item{simulations}{ A list of \code{network} objects or sparse matrices (generated using the \pkg{Matrix} package) representing simulated networks. }
\item{statistics}{ A list of functions used for comparison of observed and simulated networks. Note that the list should contain the actual functions, not a character representation of them. See \link{gof-statistics} for details. }
\item{target}{ In the \code{gof} function: A network or list of networks to which the simulations are compared. If left empty, the original networks from the \code{btergm} object \code{x} are used as observed networks. In the \code{createGOF} function: a list of sparse matrices (generated using the \pkg{Matrix} package) or a list of \code{network} objects (generated using the \pkg{network} package). The simulations are compared against these target networks. }
\item{structzero}{ Which value was used for structural zeros (usually nodes that have dropped out of the network or have not yet joined the network) in the dependent variable/network? These nodes are removed from the observed network and the simulations before comparison. Usually, the value \code{10} is used for structural zeros in Siena. }
\item{varName}{ The variable name that denotes the dependent networks in the Siena model. }
\item{verbose}{ Print details? }
\item{...}{ Arbitrary further arguments to be passed on to the statistics. See also the help page for the \link{gof-statistics}. }
}
\references{
Leifeld, Philip, Skyler J. Cranmer and Bruce A. Desmarais (2017): 
Temporal Exponential Random Graph Models with btergm: Estimation 
and Bootstrap Confidence Intervals. \emph{Journal of Statistical 
Software} 83(6): 1-36. \url{http://dx.doi.org/10.18637/jss.v083.i06}.

Leifeld, Philip and Skyler J. Cranmer (2014): A Theoretical and 
Empirical Comparison of the Temporal Exponential Random Graph 
Model and the Stochastic Actor-Oriented Model. Paper presented 
at the 7th Political Networks Conference, McGill University, 
Montreal, Canada, May 30. \url{http://arxiv.org/abs/1506.06696}.
}
\seealso{
\link{btergm-package} \link{btergm} \link{simulate.btergm} \link[ergm]{simulate.formula} \link{gof} \link{gof-statistics} \link{gof-plot}
}
\examples{
\dontrun{
# First, create data and fit a TERGM...
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

# Then assess the goodness of fit:
g <- gof(fit, statistics = c(triad.directed, esp, maxmod.modularity, 
    rocpr), nsim = 50)
g
plot(g)  # see ?"gof-plot" for details

# createGOF can also be used with user-supplied simulations:
library("statnet")
data(florentine)
gest <- ergm(flomarriage ~ edges + absdiff("wealth"))
sim <- simulate(gest, nsim = 50)
g <- createGOF(sim, list(flomarriage), statistics = c(esp, ideg), roc = FALSE)
g
plot(g)

# The help page for the Knecht dataset (?knecht) contains another example.
}
}
\author{
Philip Leifeld (\url{http://www.philipleifeld.com})
}
\keyword{methods}
\keyword{gof}
