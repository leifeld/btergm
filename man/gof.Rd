\name{gof-methods}
\docType{methods}
\alias{gof-methods}
\alias{gof,btergm-method}
\alias{gof,ergm-method}
\alias{gof,network-method}
\alias{gof,matrix-method}
\alias{gof,sienaAlgorithm-method}
\alias{gof,sienaModel-method}
\alias{gof}
\alias{gof.btergm}
\alias{gof.sienaAlgorithm}
\alias{gof.sienaModel}
\alias{gof.network}
\alias{gof.matrix}
\alias{gof.ergm}
\title{Conduct Goodness-of-Fit Diagnostics on ERGMs, TERGMs, SAOMs, and logit models}
\description{
Assess goodness of fit and degeneracy of btergm and other network models.
}
\usage{
\S4method{gof}{btergm}(object, target = NULL, 
    formula = getformula(object), nsim = 100, MCMC.interval = 1000, 
    MCMC.burnin = 10000, parallel = c("no", "MPI", "SOCK"), 
    ncpus = 1, cl = NULL, classicgof = TRUE, rocprgof = TRUE, 
    checkdegeneracy = TRUE, statistics = c("dsp", "esp", "geodist", 
    "degree", "idegree", "odegree", "kstar", "istar", "ostar"), 
    pr.impute = "poly4", verbose = TRUE, ...)

\S4method{gof}{ergm}(object, target = NULL, 
    formula = getformula(object), nsim = 100, MCMC.interval = 1000, 
    MCMC.burnin = 10000, parallel = c("no", "MPI", "SOCK"), 
    ncpus = 1, cl = NULL, classicgof = TRUE, rocprgof = TRUE, 
    checkdegeneracy = TRUE, statistics = c("dsp", "esp", "geodist", 
    "degree", "idegree", "odegree", "kstar", "istar", "ostar"), 
    pr.impute = "poly4", verbose = TRUE, ...)

\S4method{gof}{sienaAlgorithm}(object, siena.data, siena.effects, 
    predict.period = NULL, nsim = 50, parallel = c("no", "multicore", 
    "snow"), ncpus = 1, cl = NULL, target.na = NA, 
    target.na.method = "remove", target.structzero = 10, 
    classicgof = TRUE, rocprgof = TRUE, statistics = c("dsp", "esp", 
    "geodist", "degree", "idegree", "odegree", "kstar", "istar", 
    "ostar"), pr.impute = "poly4", ...)

\S4method{gof}{sienaModel}(object, siena.data, siena.effects, 
    predict.period = NULL, nsim = 50, parallel = c("no", "multicore", 
    "snow"), ncpus = 1, cl = NULL, target.na = NA, 
    target.na.method = "remove", target.structzero = 10, 
    classicgof = TRUE, rocprgof = TRUE, statistics = c("dsp", "esp", 
    "geodist", "degree", "idegree", "odegree", "kstar", "istar", 
    "ostar"), pr.impute = "poly4", ...)

\S4method{gof}{network}(object, covariates, coef, target = NULL, 
    nsim = 100, mcmc = FALSE, MCMC.interval = 1000, 
    MCMC.burnin = 10000, parallel = c("no", "MPI", "SOCK"), 
    ncpus = 1, cl = NULL, classicgof = TRUE, rocprgof = TRUE, 
    statistics = c("dsp", "esp", "geodist", "degree", "idegree", 
    "odegree", "kstar", "istar", "ostar"), pr.impute = "poly4", 
    verbose = TRUE, ...)

\S4method{gof}{matrix}(object, covariates, coef, target = NULL, 
    nsim = 100, mcmc = FALSE, MCMC.interval = 1000, 
    MCMC.burnin = 10000, parallel = c("no", "MPI", "SOCK"), 
    ncpus = 1, cl = NULL, classicgof = TRUE, rocprgof = TRUE, 
    statistics = c("dsp", "esp", "geodist", "degree", "idegree", 
    "odegree", "kstar", "istar", "ostar"), pr.impute = "poly4", 
    verbose = TRUE, ...)
}
\details{
The generic \code{gof} function provides goodness-of-fit measures and degeneracy checks for \code{btergm}, \code{ergm}, SAOM, and custom dyadic-independent models. Three different types of GOF/degeneracy assessment are possible with this function:

(1) Classic statnet-type GOF assessment by comparing summary statistics of observed and simulated networks. The \code{gof} function has six built-in statistics: dyad-wise shared partners (dsp), edge-wise shared partners (esp), degree (for undirected networks only), indegree (for directed networks only), outdegree (for directed networks only), and geodesic distances. The comparison can be plotted using boxplots for the simulations and lines for the observed network(s) or printed using t-tests (testing whether simulated and observed networks are significantly different for all values in the distributions of the summary statistics).

(2) An assessment of the classification performance using receiver operating characteristics (ROC) and precision-recall (PR) curves as well as the area under the curve (AUC) for the ROC curve.

(3) For bootstrapped TERGMs: A degeneracy check by comparing the global statistics of simulated networks to those of the observed networks at each observed time step. If the global statistics differ significantly, this is indicated by small p values. If there are many significant results, this indicates degeneracy.

For all three types of GOF assessment, by default, in-sample predictive performance is assessed by comparing all observed networks to all simulations from the same networks (just like in the \pkg{ergm} package, but aggregated over several time steps). If an observed network or a list of observed networks is provided as the \code{target} argument, the simulations are compared to these networks instead. This is useful for out-of-sample prediction. If a formula is provided, the simulations are based on the networks and covariates specified in the formula. This is helpful in situations where complex out-of-sample predictions have to be evaluated. A usage scenario could be to simulate from a network at time \code{t} (provided through the \code{formula} argument) and compare to an observed network at time \code{t + 1} (the \code{target} argument). This can be done, for example, to assess predictive performance between time steps of the original networks, or to check whether the model performs well with regard to a newly measured network given the old data from the previous time step.

Predictive fit can also be assessed for stochastic actor-oriented models (SAOM) as implemented in the \pkg{RSiena} package. After compiling the usual objects (model, data, effects), one of the time steps can be predicted based on the previous time step and the SAOM using the \code{sienaAlgorithm} (for \pkg{RSiena} >= 1.1-227) or \code{sienaModel} (for \pkg{RSiena} < 1.1-227) method of the \code{gof} function.

The \code{gof} methods for networks and matrices serve to assess the goodness of fit of a dyadic-independence model. To do this, the method requires a vector of coefficients (one coefficient for the intercept or \code{edges} term and one coefficient for each covariate), a list of covariates (in matrix or network shape), and a dependent network or matrix. This is useful for assessing the goodness of fit of QAP-adjusted logistic regression models (as implemented in the \code{netlogit} function in the \pkg{sna} package) or other dyadic-independence models, such as models fitted using \code{glm}. Note that this method only works with cross-sectional models and does not accept lists of networks as input data.

See also the \code{\link{plot.btergmgof}} help page for details on the plotting and printing options for GOF assessment.
}
\arguments{
\item{object}{ A \code{btergm}, \code{ergm}, \code{sienaAlgorithm}, or \code{sienaModel} object (for the \code{btergm}, \code{ergm}, \code{sienaAlgorithm}, and \code{sienaModel} methods, respectively). Or a network object or matrix (for the \code{network} and \code{matrix} methods, respectively). }
\item{siena.data}{ An object of the class \code{siena}, which is usually created using the \code{sienaDataCreate} function in the \code{RSiena} package. }
\item{siena.effects}{ An object of the class \code{sienaEffects}, which is usually created using the \code{getEffects()} and the \code{includeEffects()} function in the \code{RSiena} package. }
\item{predict.period}{ Which time period should be predicted? By default, the last time period is predicted based on the last simulation of the second-last time period. The time period can be provided as a numeric, e.g., \code{predict.period = 4} for predicting the fourth network. }
\item{target}{ A network or list of networks to which the simulations are compared. If left empty, the original networks from the \code{btergm} object \code{x} are used as observed networks. }
\item{formula}{ A model formula from which networks are simulated for comparison. By default, the formula from the \code{btergm} object \code{x} is used. It is possible to hand over a formula with only a single response network and/or dyad or edge covariates or with lists of response networks and/or covariates. It is also possible to use indices like \code{networks[[4]]} or \code{networks[3:5]} inside the formula. }
\item{nsim}{ The number of networks to be simulated at each time step. Example: If there are six time steps in the \code{formula} and \code{nsim = 100}, a total of 600 new networks is simulated. }
\item{MCMC.interval}{ Internally, this package uses the simulation facilities of the \pkg{ergm} package to create new networks against which to compare the original network(s) for goodness-of-fit assessment. This argument sets the MCMC interval to be passed over to the simulation command. The default value is \code{1000}, which means that every 1000th simulation outcome from the MCMC sequence is used. There is no general rule of thumb on the selection of this parameter, but if the results look suspicious (e.g., when the model fit is perfect), increasing this value may be helpful. }
\item{MCMC.burnin}{ Internally, this package uses the simulation facilities of the \pkg{ergm} package to create new networks against which to compare the original network(s) for goodness-of-fit assessment. This argument sets the MCMC burnin to be passed over to the simulation command. The default value is \code{10000}. There is no general rule of thumb on the selection of this parameter, but if the results look suspicious (e.g., when the model fit is perfect), increasing this value may be helpful. }
\item{parallel}{ Use multiple cores in a computer or nodes in a cluster to speed up the simulations. The default value \code{"no"} means parallel computing is switched off. If \code{"multicore"} is used (only available for \code{sienaAlgorithm} and \code{sienaModel} objects), the \code{mclapply} function from the \pkg{parallel} package (formerly in the \pkg{multicore} package) is used for parallelization. This should run on any kind of system except MS Windows because it is based on forking. It is usually the fastest type of parallelization. If \code{"snow"} is used (only available for \code{sienaAlgorithm} and \code{sienaModel} objects), the \code{parLapply} function from the \pkg{parallel} package (formerly in the \pkg{snow} package) is used for parallelization. This should run on any kind of system including cluster systems and including MS Windows. It is slightly slower than the former alternative if the same number of cores is used. However, \code{"snow"} provides support for MPI clusters with a large amount of cores, which \pkg{multicore} does not offer (see also the \code{cl} argument). If \code{"MPI"} is used (only available for \code{btergm} and \code{ergm} objects), MPI parallelization as implemented in the \pkg{ergm} package is used. And if \code{"SOCK"} is used (only available for \code{btergm} and \code{ergm} objects), a SOCK cluster as implemented in the \pkg{ergm} package (via the deprecated \pkg{snow} package) is used for parallelization. Note that \code{"multicore"} and \code{"SOCK"} will only work if all cores are on the same node. For example, if there are three nodes with eight cores each, a maximum of eight CPUs can be used. }
\item{ncpus}{ The number of CPU cores used for parallel simulations (only if \code{parallel} is activated). If the number of cores should be detected automatically on the machine where the code is executed, one can try the \code{detectCores()} function from the \pkg{parallel} package. On some HPC clusters, the number of available cores is saved as an environment variable; for example, if MOAB is used, the number of available cores can sometimes be accessed using \code{Sys.getenv("MOAB_PROCCOUNT")}, depending on the implementation. Note that the maximum number of connections in a single R session (i.e., to other cores or for opening files etc.) is 128, so fewer than 128 cores should be used at a time. }
\item{cl}{ An optional \pkg{parallel} or \pkg{snow} cluster for use if \code{parallel = "snow"}. If not supplied, a cluster on the local machine is created temporarily. }
\item{target.na}{ Which value was used for missing data in the dependent variable? }
\item{target.na.method}{ How should missing data be handled when comparing the simulations to the empirical (= observed) network? Two options are possible: \code{remove} drops nodes with missing ties both from the simulations (after running the simulations) and from the observed network before the comparison. \code{fillmode} replaces missing values by the mode of the network matrix (usually \code{0}). }
\item{target.structzero}{ Which value was used for structural zeroes (usually nodes which have dropped out of the network or have not yet joined the network) in the dependent variable? These nodes are removed from the observed network and the simulations before comparison. }
\item{classicgof}{ If \code{classicgof = TRUE} is set, the classic statnet-style goodness-of-fit comparison is conducted. This means that shared-partner statistics, the geodesic distance distribution and the degree distribution are compared between observed and simulated networks. The results can be plotted as boxplots or printed as tables. Note that the \code{classicgof}, \code{rocprgof} and \code{checkdegeneracy} arguments can be used together. In that case, the resulting \code{btergmgof} object will contain all three types of GOF/degeneracy assessment. }
\item{rocprgof}{ If \code{rocprgof = TRUE} is set, the coordinates of ROC and PR curves as well as the AUC measure are stored in the resulting \code{btergmgof} object. The results can be plotted as curves or printed as tables. Note that the \code{classicgof}, \code{rocprgof} and \code{checkdegeneracy} arguments can be used together. In that case, the resulting \code{btergmgof} object will contain all three types of GOF/degeneracy assessment. }
\item{checkdegeneracy}{ If \code{checkdegeneracy = TRUE} is set, the global statistics of the observed and simulated networks are compared for each observed time step separately. Frequent significant deviations indicate degeneracy. The results can be printed as tables. Note that the \code{classicgof}, \code{rocprgof} and \code{checkdegeneracy} arguments can be used together. In that case, the resulting \code{btergmgof} object will contain all three types of GOF/degeneracy assessment. }
\item{statistics}{ A character vector of auxiliary statistics used for comparison of observed and simulated networks. Valid values are \code{"dsp"}, \code{"esp"}, \code{"geodist"}, \code{"degree"}, \code{"idegree"}, \code{"odegree"}, \code{"kstar"}, \code{"istar"}, and \code{"ostar"}. By default, all of these statistics are used. }
\item{pr.impute}{ In some cases, the first precision value of the precision-recall curve is undefined. The \code{pr.impute} argument serves to impute this missing value to ensure that the AUC-PR value is not severely biased.

Possible values are \code{"no"} for no imputation, \code{"one"} for using a value of \code{1.0}, \code{"second"} for using the next (= adjacent) precision value, \code{"poly1"} for fitting a straight line through the remaining curve to predict the first value, \code{"poly2"} for fitting a second-order polynomial curve etc. until \code{"poly9"}

Warning: this is a pragmatic solution. Please double-check whether the imputation makes sense. This can be checked by plotting the resulting \code{btergmgof} object and using the \code{pr.poly} argument to plot the predicted curve on top of the actual PR curve. }
\item{covariates}{ A list of matrices or network objects that serve as covariates for the dependent network. The covariates in this list are automatically added to the formula as \code{edgecov} terms. }
\item{coef}{ A vector of coefficients. }
\item{mcmc}{ Should statnet's MCMC methods be used for simulating new networks? If \code{mcmc = FALSE}, new networks are simulated based on predicted tie probabilities of the regression equation. }
\item{verbose}{ Print details? }
\item{...}{ Arbitrary further arguments are handed over to the \code{\link[ergm]{simulate.formula}} function or the \code{\link[RSiena]{siena07}} function. For details, refer to the help page of these functions. }
}
\seealso{
\link{btergm-package} \link{btergm} \link{simulate.btergm} \link[ergm]{simulate.formula} \link{plot.btergmgof}
}
\keyword{methods}
\keyword{gof}
