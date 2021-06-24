#' Goodness-of-fit diagnostics for ERGMs, TERGMs, SAOMs, and logit models
#' 
#' Assess goodness of fit of \code{btergm} and other network models.
#' 
#' The generic \code{gof} function provides goodness-of-fit measures and
#' degeneracy checks for \code{btergm}, \code{mtergm}, \code{tbergm},
#' \code{ergm}, \code{sienaFit}, and custom dyadic-independent models. The user
#' can provide a list of network statistics for comparing simulated networks
#' based on the estimated model with the observed network(s). See
#' \code{\link{gof-statistics}}. The objects created by these methods can be
#' displayed using various plot and print methods (see \code{\link{gof-plot}}).
#' 
#' In-sample GOF assessment is the default, which means that the same time steps
#' are used for creating simulations and for comparison with the observed
#' network(s). It is possible to do out-of-sample prediction by specifying a
#' (list of) target network(s) using the \code{target} argument. If a formula is
#' provided, the simulations are based on the networks and covariates specified
#' in the formula. This is helpful in situations where complex out-of-sample
#' predictions have to be evaluated. A usage scenario could be to simulate from
#' a network at time \code{t} (provided through the \code{formula} argument) and
#' compare to an observed network at time \code{t + 1} (the \code{target}
#' argument). This can be done, for example, to assess predictive performance
#' between time steps of the original networks, or to check whether the model
#' performs well with regard to a newly measured network given the old data from
#' the previous time step.
#' 
#' Predictive fit can also be assessed for stochastic actor-oriented models
#' (SAOM) as implemented in the \pkg{RSiena} package. After compiling the usual
#' objects (model, data, effects), one of the time steps can be predicted based
#' on the previous time step and the SAOM using the \code{sienaFit} method of
#' the \code{gof} function. By default, however, within-sample fit is used for
#' SAOMs, just like for (T)ERGMs.
#' 
#' The \code{gof} methods for networks and matrices serve to assess the goodness
#' of fit of a dyadic-independence model. To do this, the method requires a
#' vector of coefficients (one coefficient for the intercept or \code{edges}
#' term and one coefficient for each covariate), a list of covariates (in matrix
#' or network shape), and a dependent network or matrix. This is useful for
#' assessing the goodness of fit of QAP-adjusted logistic regression models (as
#' implemented in the \code{netlogit} function in the \pkg{sna} package) or
#' other dyadic-independence models, such as models fitted using \code{glm}.
#' Note that this method only works with cross-sectional models and does not
#' accept lists of networks as input data.
#' 
#' The \code{createGOF} function is used internally by the \code{gof} function
#' in order to create a \code{gof} object from a list of simulated networks and
#' a list of target networks to compare against. It can also be used directly by
#' the end user if the user wants to supply lists of simulated and target
#' networks from other sources.
#' 
#' @param cl An optional \pkg{parallel} or \pkg{snow} cluster for use if
#'   \code{parallel = "snow"}. If not supplied, a cluster on the local machine
#'   is created temporarily.
#' @param coef A vector of coefficients.
#' @param covariates A list of matrices or network objects that serve as
#'   covariates for the dependent network. The covariates in this list are
#'   automatically added to the formula as \code{edgecov} terms.
#' @param formula A model formula from which networks are simulated for
#'   comparison. By default, the formula from the \code{btergm} object \code{x}
#'   is used. It is possible to hand over a formula with only a single response
#'   network and/or dyad or edge covariates or with lists of response networks
#'   and/or covariates. It is also possible to use indices like
#'   \code{networks[[4]]} or \code{networks[3:5]} inside the formula.
#' @param groupName The group name used in the Siena model.
#' @param mcmc Should statnet's MCMC methods be used for simulating new
#'   networks? If \code{mcmc = FALSE}, new networks are simulated based on
#'   predicted tie probabilities of the regression equation.
#' @param MCMC.burnin Internally, this package uses the simulation facilities of
#'   the \pkg{ergm} package to create new networks against which to compare the
#'   original network(s) for goodness-of-fit assessment. This argument sets the
#'   MCMC burnin to be passed over to the simulation command. The default value
#'   is \code{10000}. There is no general rule of thumb on the selection of this
#'   parameter, but if the results look suspicious (e.g., when the model fit is
#'   perfect), increasing this value may be helpful.
#' @param MCMC.interval Internally, this package uses the simulation facilities
#'   of the \pkg{ergm} package to create new networks against which to compare
#'   the original network(s) for goodness-of-fit assessment. This argument sets
#'   the MCMC interval to be passed over to the simulation command. The default
#'   value is \code{1000}, which means that every 1000th simulation outcome from
#'   the MCMC sequence is used. There is no general rule of thumb on the
#'   selection of this parameter, but if the results look suspicious (e.g., when
#'   the model fit is perfect), increasing this value may be helpful.
#' @param ncpus The number of CPU cores used for parallel GOF assessment (only
#'   if \code{parallel} is activated). If the number of cores should be detected
#'   automatically on the machine where the code is executed, one can try the
#'   \code{detectCores()} function from the \pkg{parallel} package. On some HPC
#'   clusters, the number of available cores is saved as an environment
#'   variable; for example, if MOAB is used, the number of available cores can
#'   sometimes be accessed using \code{Sys.getenv("MOAB_PROCCOUNT")}, depending
#'   on the implementation. Note that the maximum number of connections in a
#'   single R session (i.e., to other cores or for opening files etc.) is 128,
#'   so fewer than 128 cores should be used at a time.
#' @param nsim The number of networks to be simulated at each time step.
#'   Example: If there are six time steps in the \code{formula} and
#'   \code{nsim = 100}, a total of 600 new networks is simulated. The
#'   comparison between simulated and observed networks is only done within time
#'   steps. For example, the first 100 simulations are compared with the first
#'   observed network, simulations 101-200 with the second observed network etc.
#' @param object A \code{btergm}, \code{ergm}, or \code{sienaFit} object (for
#'   the \code{btergm}, \code{ergm}, and \code{sienaFit} methods, respectively).
#'   Or a network object or matrix (for the \code{network} and \code{matrix}
#'   methods, respectively).
#' @param outofsample Should out-of-sample prediction be attempted? If so, some
#'   additional arguments must be provided: \code{sienaData},
#'   \code{sienaEffects}, and \code{nsim}. The \code{sienaData} object must
#'   contain a base and a target network for out-of-sample prediction. The
#'   \code{sienaEffects} must contain the effects to be used for the
#'   simulations. The estimates will be taken from the estimated \code{object},
#'   and they will be injected into a new SAOM and fixed during the sampling
#'   procedure. \code{nsim} determines how many simulations are used for the
#'   out-of-sample comparison.
#' @param parallel Use multiple cores in a computer or nodes in a cluster to
#'   speed up the simulations. The default value \code{"no"} means parallel
#'   computing is switched off. If \code{"multicore"} is used (only available
#'   for \code{sienaAlgorithm} and \code{sienaModel} objects), the
#'   \code{mclapply} function from the \pkg{parallel} package (formerly in the
#'   \pkg{multicore} package) is used for parallelization. This should run on
#'   any kind of system except MS Windows because it is based on forking. It is
#'   usually the fastest type of parallelization. If \code{"snow"} is used, the
#'   \code{parLapply} function from the \pkg{parallel} package (formerly in the
#'   \pkg{snow} package) is used for parallelization. This should run on any
#'   kind of system including cluster systems and including MS Windows. It is
#'   slightly slower than the former alternative if the same number of cores is
#'   used. However, \code{"snow"} provides support for MPI clusters with a large
#'   amount of cores, which \pkg{multicore} does not offer (see also the
#'   \code{cl} argument). Note that \code{"multicore"} will only work if all
#'   cores are on the same node. For example, if there are three nodes with
#'   eight cores each, a maximum of eight CPUs can be used. Parallel computing
#'   is described in more detail on the help page of \link{btergm}.
#' @param period Which transition between time periods should be used for GOF
#'   assessment? By default, all transitions between all time periods are used.
#'   For example, if there are three consecutive networks, this will extract
#'   simulations from the transitions between 1 and 2 and between 2 and 3,
#'   respectively, and these simulations will be compared to the networks at
#'   time steps 2 and 3, respectively. The time period can be provided as a
#'   numeric, e.g., \code{period = 4} for extracting the simulations between
#'   time steps 4 and 5 (= the fourth transition) and predicting the fifth
#'   network. Values lower than 1 or larger than the number of consecutive
#'   networks minus 1 are therefore not permitted. This argument is only used if
#'   out-of-sample prediction is switched off.
#' @param sienaData An object of the class \code{siena}, which is usually
#'   created using the \code{sienaDataCreate} function in the \pkg{RSiena}
#'   package. This argument is only used for out-of-sample prediction. The
#'   object must be based on a \code{sienaDependent} object that contains two
#'   networks: the base network from which to simulate forward, and the target
#'   network which you want to predict out-of-sample. The object can contain
#'   further objects for storing covariates etc. that are necessary for
#'   estimating new networks. The best practice is to create an object that is
#'   identical to the \code{siena} object used for estimating the model, except
#'   that it contains the base and the target network instead of the dependent
#'   variable/networks.
#' @param sienaEffects An object of the class \code{sienaEffects}, which is
#'   usually created using the \code{getEffects()} and the
#'   \code{includeEffects()} functions in the \code{RSiena} package. The best
#'   practice is to provide a \code{sienaEffects} object that is identical to
#'   the object used to create the original model (that is, it should contain
#'   the same effects), except that it should be based on the \code{siena}
#'   object provided through the \code{sienaData} argument. In other words, the
#'   \code{sienaEffects} object should be based on the base and target network
#'   used for out-of-sample prediction, and it should contain the same effects
#'   as those used for the original estimation. This argument is used only for
#'   out-of-sample prediction.
#' @param simulations A list of \code{network} objects or sparse matrices
#'   (generated using the \pkg{Matrix} package) representing simulated networks.
#' @param statistics A list of functions used for comparison of observed and
#'   simulated networks. Note that the list should contain the actual functions,
#'   not a character representation of them. See \link{gof-statistics} for
#'   details.
#' @param target In the \code{gof} function: A network or list of networks to
#'   which the simulations are compared. If left empty, the original networks
#'   from the \code{btergm} object \code{x} are used as observed networks. In
#'   the \code{createGOF} function: a list of sparse matrices (generated using
#'   the \pkg{Matrix} package) or a list of \code{network} objects (generated
#'   using the \pkg{network} package). The simulations are compared against
#'   these target networks.
#' @param structzero Which value was used for structural zeros (usually nodes
#'   that have dropped out of the network or have not yet joined the network) in
#'   the dependent variable/network? These nodes are removed from the observed
#'   network and the simulations before comparison. Usually, the value \code{10}
#'   is used for structural zeros in Siena.
#' @param varName The variable name that denotes the dependent networks in the
#'   Siena model.
#' @param verbose Print details?
#' @param ... Arbitrary further arguments to be passed on to the statistics. See
#'   also the help page for the \link{gof-statistics}.
#'
#' @references
#' Leifeld, Philip, Skyler J. Cranmer and Bruce A. Desmarais (2018): Temporal
#' Exponential Random Graph Models with btergm: Estimation and Bootstrap
#' Confidence Intervals. \emph{Journal of Statistical Software} 83(6): 1--36.
#' \doi{10.18637/jss.v083.i06}.
#' 
#' Leifeld, Philip and Skyler J. Cranmer (2019): A Theoretical and Empirical
#' Comparison of the Temporal Exponential Random Graph Model and the Stochastic
#' Actor-Oriented Model. Network Science 7(1): 20--51.
#' \doi{10.1017/nws.2018.26}.
#' 
#' @docType methods
#' @aliases gof-methods gofmethods
#' @importFrom ergm simulate_formula control.simulate.formula
#' @export
setGeneric("gof", function(object, ...) standardGeneric("gof"),
           package = "btergm")

#' Reduce a statistic x nsim matrix and compute summary statistics
#' 
#' Reduce a statistic x nsim matrix and compute summary statistics.
#' 
#' A function which reduces a statistic x nsim matrix and computes summary
#' statistics with the goal of getting rid of empty rows at the end, e.g., where
#' \code{dsp(37)} or so is usually not observed.
#' 
#' @param sim A simulated matrix of a certain type of statistics.
#' @param obs An observed matrix of a certain type of statistics.
#' 
#' @return Object containing the summary statistics.
#' 
#' @importFrom stats pnorm t.test
#' 
#' @noRd
reduce.matrix <- function(sim, obs) {
  
  numsim <- ncol(as.matrix(sim))
  numobs <- ncol(as.matrix(obs))
  xlabels <- rownames(obs)
  
  # if geodist statistic: put aside the last 'Inf' row
  if (is.null(rownames(sim)) || rownames(sim)[nrow(sim)] == "Inf") {
    geo <- TRUE
    inf.sim <- sim[nrow(sim), ]  # put aside this row for now and reuse later
    sim <- sim[-nrow(sim), ]
    if (is.matrix(obs)) {
      inf.obs <- obs[nrow(obs), ]
      obs <- matrix(obs[-nrow(obs), ], ncol = numobs)
    } else {
      inf.obs <- obs[length(obs)]
      obs <- matrix(obs[-length(obs)], ncol = numobs)
    }
  } else {
    geo <- FALSE
  }
  
  # find first empty row for simulation matrix
  sim.rs <- rowSums(sim)
  sim.remove <- length(sim.rs)  # at which row index can removal start?
  for (i in (length(sim.rs) - 1):1) {
    if (sim.rs[i] == 0 && sim.remove == (i + 1)) {
      sim.remove <- i  # remember which is the first empty row
    }
  }
  
  if (!is.matrix(obs)) {  # one network is compared
    obs.remove <- length(obs)
    for (i in (length(obs) - 1):1) {
      if (obs[i] == 0 && obs.remove == (i + 1)) {
        obs.remove <- i
      }
    }
  } else {  # several networks are compared
    obs.rs <- rowSums(obs)
    obs.remove <- length(obs.rs)
    for (i in (length(obs.rs) - 1):1) {
      if (obs.rs[i] == 0 && obs.remove == (i + 1)) {
        obs.remove <- i
      }
    }
  }
  rem <- max(c(obs.remove, sim.remove), na.rm = TRUE)  # which one is longer?
  
  # remove unnecessary rows
  if (!is.matrix(obs)) {  # get rid of empty observations or rows of obs
    obs <- matrix(obs[-(rem:length(obs))], ncol = numobs)
  } else {
    obs <- matrix(obs[-(rem:nrow(obs)), ], ncol = numobs)
  }
  sim <- matrix(sim[-(rem:nrow(sim)), ], ncol = numsim)  # same for sim stats
  
  if (nrow(obs) < rem) {
    for (i in (nrow(obs) + 1):rem) {
      obs <- rbind(obs, rep(0, ncol(obs)))
    }
  }
  if (nrow(sim) < rem) {
    for (i in (nrow(sim) + 1):rem) {
      sim <- rbind(sim, rep(0, ncol(sim)))
    }
  }
  
  # for geodist, add Inf row again
  if (geo == TRUE) {
    sim <- rbind(sim, inf.sim)
    obs <- rbind(obs, inf.obs)
    rownames(sim) <- c(1:(nrow(sim) - 1), "Inf")
  }
  
  # create final object which will contain the raw simulations + the comparison
  reducedobject <- list()
  reducedobject$sim <- as.data.frame(sim)
  
  rownames(sim) <- NULL
  rownames(obs) <- NULL
  
  # compute means, p values, etc. and put them in a data frame
  x <- matrix()
  if (ncol(obs) == 1 || ncol(sim) == 1) {  # compute all the summary statistics
    x.obs <- obs
    x.mean <- apply(sim, 1, mean)
    x.min <- apply(sim, 1, min)
    x.max <- apply(sim, 1, max)
    x.median <- apply(sim, 1, median)
    zval <- (x.mean - x.obs) / sd(x.mean)
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
    x.pval <- pval
    x <- data.frame(x.obs, x.mean, x.median, x.min, x.max, x.pval)
    colnames(x) <- c("obs", "sim: mean", "median", "min", "max", "Pr(>z)")
  } else {  # for several target networks, compute distribution
    x.obs.mean <- apply(obs, 1, mean)
    x.obs.min <- apply(obs, 1, min)
    x.obs.max <- apply(obs, 1, max)
    x.obs.median <- apply(obs, 1, median)
    x.mean <- apply(sim, 1, mean)
    x.min <- apply(sim, 1, min)
    x.max <- apply(sim, 1, max)
    x.median <- apply(sim, 1, median)
    x.pval <- numeric()
    for (i in 1:nrow(sim)) {
      tryCatch(
        expr = {
          x.pval[i] <- t.test(obs[i, ], sim[i, ])$p.value  # compare group means
        }, 
        error = function(e) {
          x.pval[i] <- 1  # if both are 0, a t.test cannot be computed...
        }, 
        finally = {}
      )
      
      if (is.nan(x.pval[i])) {  # geodist contains "Inf"
        x.pval[i] <- 1
      }
    }
    x <- data.frame(x.obs.mean, x.obs.median, x.obs.min, x.obs.max, x.mean, 
                    x.median, x.min, x.max, x.pval)
    colnames(x) <- c("obs: mean", "median", "min", "max", 
                     "sim: mean", "median", "min", "max", "Pr(>z)")
  }
  if (geo == TRUE) {
    rownames(x) <- c(xlabels[c(1:(nrow(x) - 1))], "Inf")
  } else {
    rownames(x) <- xlabels[1:nrow(x)]
  }
  rownames(reducedobject$sim) <- rownames(x)
  
  reducedobject$comparison <- as.data.frame(x)
  
  return(reducedobject)
}

#' @rdname gof
#' @importFrom parallel makeCluster mclapply clusterEvalQ parSapply stopCluster
#' @export
createGOF <- function(simulations,
                      target,
                      statistics = c(dsp,
                                     esp,
                                     deg,
                                     ideg,
                                     geodesic,
                                     rocpr,
                                     walktrap.modularity),
                      parallel = "no", 
                      ncpus = 1,
                      cl = NULL,
                      verbose = TRUE,
                      ...) {
  
  # prepare parallel processing
  stop.parallel <- FALSE
  if (parallel[1] == "snow" && is.null(cl) && ncpus > 1) {
    cl <- makeCluster(ncpus)
    stop.parallel <- TRUE
  }
  
  # make sure target list is a list of sparse matrices
  if (!"list" %in% class(target)) {
    stop(paste("The target network(s) must be provided as a list of network", 
               "objects or sparse matrix objects (using the Matrix package)."))
  }
  target <- lapply(target, function(x) {
    if (!class(x)[1] %in% c("dgCMatrix", "dgTMatrix", "dsCMatrix", "dsTMatrix", 
                            "dgeMatrix")) {
      mat <- as.matrix(x)
      mat <- Matrix(mat)
      return(mat)
    } else {
      return(x)
    }
  })
  
  # make sure simulation list is a list of network objects
  if ("network.list" %in% class(simulations)) {
    simulations <- lapply(simulations, function(x) Matrix(as.matrix(x)))
  } else if (!"list" %in% class(simulations)) {
    stop("The simulations must be provided as a list of network objects.")
  }
  
  # go through statistics functions and compute results
  goflist <- list()
  for (z in 1:length(statistics)) {
    args <- names(formals(statistics[[z]]))  # determine number of arguments
    if ("sim" %in% args && "obs" %in% args) {
      args <- 2
    } else {
      args <- 1
    }
    if (args == 1) {  # either sim or obs
      tryCatch(
        expr = {
          label <- suppressMessages(attributes(statistics[[z]](
            simulations[[1]])))$label
          if (verbose == TRUE) {
            message(paste("Processing statistic:", label))
          }
          
          if (parallel[1] == "no") {
            simulated <- suppressMessages(sapply(simulations, statistics[[z]], 
                                                 ...))
            observed <- suppressMessages(sapply(target, statistics[[z]], ...))
          } else if (parallel[1] == "multicore") {
            test <- suppressMessages(statistics[[z]](simulations[[1]]))
            if (is.numeric(test) && length(test) == 1) {
              simulated <- suppressMessages(unlist(mclapply(simulations, 
                                                            statistics[[z]], ..., mc.cores = ncpus)))
              observed <- suppressMessages(unlist(mclapply(target, 
                                                           statistics[[z]], ..., mc.cores = ncpus)))
            } else {  # no mcsapply available because different length vectors
              simulated <- suppressMessages(mclapply(simulations, 
                                                     statistics[[z]], ..., mc.cores = ncpus))
              observed <- suppressMessages(mclapply(target, statistics[[z]], 
                                                    ..., mc.cores = ncpus))
              max.length.sim <- max(sapply(simulated, length), na.rm = TRUE)
              max.length.obs <- max(sapply(observed, length), na.rm = TRUE)
              max.length <- max(max.length.sim, max.length.obs, na.rm = TRUE)
              simulated <- sapply(simulated, function(x) {
                c(x, rep(0, max.length - length(x)))
              })
              observed <- sapply(observed, function(x) {
                c(x, rep(0, max.length - length(x)))
              })
            }
          } else {
            clusterEvalQ(cl, library("ergm"))
            simulated <- suppressMessages(parSapply(cl = cl, simulations, 
                                                    statistics[[z]], ...))
            observed <- suppressMessages(parSapply(cl = cl, target, 
                                                   statistics[[z]], ...))
          }
          
          # if simulations have different dimensions, convert list to matrix
          if ("list" %in% class(simulated)) {
            lengths <- sapply(simulated, length)
            l <- max(lengths)
            index <- which(lengths == l)[1]
            rn <- names(simulated[[index]])
            simulated <- sapply(simulated, function(x) {
              c(x, rep(0, l - length(x)))
            })
            rownames(simulated) <- rn
          }
          if ("list" %in% class(observed)) {
            lengths <- sapply(observed, length)
            l <- max(lengths)
            index <- which(lengths == l)[1]
            rn <- names(observed[[index]])
            observed <- sapply(observed, function(x) {
              c(x, rep(0, l - length(x)))
            })
            rownames(observed) <- rn
          }
          
          gofobject <- list()
          gofobject$label <- label
          if (is.matrix(simulated)) {  # boxplot-type GOF
            reduced <- reduce.matrix(simulated, observed)
            gofobject$type <- "boxplot"
            gofobject$stats <- reduced$comparison
            gofobject$raw <- Matrix::Matrix(as.matrix(reduced$sim))
            class(gofobject) <- "boxplot"
          } else if (is.numeric(simulated)) {  # density-type GOF
            gofobject$type <- "univariate"
            gofobject$obs <- observed
            gofobject$sim <- simulated
            class(gofobject) <- "univariate"
          }
          goflist[[length(goflist) + 1]] <- gofobject
          names(goflist)[length(goflist)] <- label
        }, 
        error = function(e) {
          if (verbose == TRUE) {
            cat(paste("  Skipping statistic for the following reason:", e))
          }
        }, 
        finally = {}
      )
    } else if (args == 2) {  # sim and obs; ROCPR-type GOF
      tryCatch(
        expr = {
          label <- "Tie prediction"
          if (verbose == TRUE) {
            message(paste("Processing statistic:", label))
          }
          gofobject <- statistics[[z]](sim = simulations, obs = target, ... = ...)
          goflist[[length(goflist) + 1]] <- gofobject
          names(goflist)[length(goflist)] <- label
        }, 
        error = function(e) {
          cat(paste("  Skipping this statistic for the following reason:", 
                    e))
        }, 
        finally = {}
      )
    }
  }
  class(goflist) <- "gof"
  if (stop.parallel == TRUE) {
    stopCluster(cl)
  }
  return(goflist)
}

#' @noRd
gof.btergm <- function(object, target = NULL, formula = getformula(object), 
                       nsim = 100, MCMC.interval = 1000, MCMC.burnin = 10000,
                       parallel = c("no", "multicore", "snow"), ncpus = 1,
                       cl = NULL, statistics = c(dsp, esp, deg, ideg, geodesic,
                                                 rocpr, walktrap.modularity),
                       verbose = TRUE, ...) {
  
  if (nsim < 2) {
    stop("The 'nsim' argument must be greater than 1.")
  }
  
  if (is.function(statistics)) {
    statistics <- c(statistics)
  }
  
  # prepare parallel processing; translate options into statnet arguments
  if (is.null(ncpus) || ncpus == 0) {
    ncpus <- 1
  }
  if (!parallel[1] %in% c("no", "multicore", "snow")) {
    parallel <- "no"
    warning("'parallel' argument not recognized. Using 'no' instead.")
  }
  
  if (verbose == TRUE) {
    if (parallel[1] == "no") {
      parallel.msg <- "on a single computing core."
    } else if (parallel[1] == "multicore") {
      parallel.msg <- paste("using multicore forking on", ncpus, "cores.")
    } else if (parallel[1] == "snow") {
      parallel.msg <- paste("using parallel processing on", ncpus, "cores.")
    }
    message("\nStarting GOF assessment ", parallel.msg, "...")
  }
  
  # call tergmprepare and integrate results as a child environment in the chain
  if ("btergm" %in% class(object)) {
    l <- tergmprepare(formula = formula, offset = object@offset, 
                      verbose = verbose)
    offset <- object@offset
    form <- as.formula(l$form, env = environment())
  } else {
    l <- tergmprepare(formula = formula, offset = FALSE, verbose = FALSE)
    offset <- TRUE
    form <- as.formula(paste(l$form, "+ offset(edgecov(offsmat[[i]]))"), 
                       env = environment())
  }
  for (i in 1:length(l$covnames)) {
    assign(l$covnames[i], l[[l$covnames[i]]])
  }
  assign("offsmat", l$offsmat)
  
  # check and rearrange target network(s)
  if (is.null(target)) {
    if (verbose == TRUE) {
      message(paste("\nNo 'target' network(s) provided. Using networks on the",
                    "left-hand side of the model formula as observed networks.\n"))
    }
    target <- l$networks
  } else if (is.network(target) || is.matrix(target)) {
    target <- list(target)
    if (verbose == TRUE) {
      message("\nOne observed ('target') network was provided.\n")
    }
  } else if ("list" %in% class(target)) {
    if (verbose == TRUE) {
      message(paste("\n", length(target), "observed ('target') networks were",
                    "provided.\n"))
    }
  } else {
    stop("'target' must be a network, matrix, or list of matrices or networks.")
  }
  
  # extract coefficients from object
  if ("tbergm" %in% class(object)) {
    coefs <- object@bergm$Theta[sample(1:nrow(object@bergm$Theta), nsim), ,
                                drop = FALSE]
    if (isFALSE(offset)) {
      coefs <- coefs[, -ncol(coefs)]
    }
  } else {
    if (offset == TRUE) {
      coefs <- c(coef(object), -Inf)  # -Inf for offset matrix
    } else {
      coefs <- coef(object)
    }
  }
  
  # adjust formula at each step, and simulate networks
  sim <- list()
  degen <- list()
  for (index in 1:l$time.steps) {
    i <- index  # index 'i' is used in formula construction!
    # simulations for statnet-style and rocpr GOF
    if (verbose == TRUE) {
      if ("btergm" %in% class(object) || "mtergm" %in% class(object) || "tbergm" %in% class(object)) {
        f.i <- gsub("\\[\\[i\\]\\]", paste0("[[", index, "]]"), 
                    paste(deparse(form), collapse = ""))
        f.i <- gsub("\\s+", " ", f.i)
        if ("btergm" %in% class(object)) {
          f.i <- gsub("^networks", l$lhs.original, f.i)
        }
      } else if ("ergm" %in% class(object)) {
        f.i <- paste(deparse(formula), collapse = "")
        f.i <- gsub("\\s+", " ", f.i)
      } else {
        stop(paste("Unknown object type:", paste(class(object), collapse = ", ")))
      }
      message(paste("Simulating", nsim, 
                    "networks from the following formula:\n", f.i, "\n"))
    }
    
    # check for mismatch in coefficients and data
    if (!"tbergm" %in% class(object) && length(coefs[!is.infinite(coefs)]) > ncol(ergm::ergmMPLE(form)$predictor)) {
      mismatch <- names(coefs)[which(!names(coefs) %in% colnames(ergm::ergmMPLE(form)$predictor))]
      mismatch <- mismatch[mismatch != ""]
      msg <- paste0("At t = ", i, ", at least one of the covariates has missing ",
                    "levels for which coefficients but no data were available. ",
                    "This can happen, for example, when a categorical node ",
                    "variable has levels 0 to 2, but at a specific time step ",
                    "only levels 0 and 2 are found while level 1 is absent from ",
                    "the network. In this case, the coefficient from the TERGM ",
                    "for that absent level is present, causing a mismatch ",
                    "between the number of coefficients and the data structure. ",
                    "At the moment, this problem cannot be fixed because it ",
                    "would require changes in another package outside of the ",
                    "control of the btergm package authors. The only known fix ",
                    "at this point is to make sure that all variable levels are ",
                    "present at least once at each time point. Here is a list of ",
                    "model terms which cause the mismatch at t = ", i, ":")
      for (j in 1:length(mismatch)) {
        msg <- paste(msg, mismatch[j])
        if (j == length(mismatch)) {
          msg <- paste0(msg, ". Aborting now.")
        } else {
          msg <- paste0(msg, "; ")
        }
      }
      message(paste0(strwrap(msg), collapse = "\n"))
      stop("Coefficient mismatch error.")
    }
    
    tryCatch(
      expr = {
        if ("tbergm" %in% class(object)) {
          sim[[index]] <- apply(coefs, 1, function(x) {
            ergm::simulate_formula(form,
                                   nsim = 1,
                                   coef = x,
                                   constraints = ~ .,
                                   control = ergm::control.simulate.formula(MCMC.interval = MCMC.interval,
                                                                            MCMC.burnin = MCMC.burnin))
          })
        } else {
          sim[[index]] <- ergm::simulate_formula(form,
                                                 nsim = nsim,
                                                 coef = coefs,
                                                 constraints = ~ .,
                                                 control = ergm::control.simulate.formula(MCMC.interval = MCMC.interval,
                                                                                          MCMC.burnin = MCMC.burnin))
        }
      }, 
      error = function(e) {
        if (grep("elements, while the model requires", as.character(e))) {
          message("Error: The number of ERGM coefficients does not correspond to the number of parameters the simulation function is expecting. This can occur if you employ curved model terms, such as gwesp, gwdegree etc., without explicitly adding ', fixed = TRUE' to the model term, for example '+ gwesp(0.5, fixed = TRUE). The MPLE-based estimation strategy currently does not permit curved terms (i.e., the parameters must be fixed). Please add ', fixed = TRUE' to the respective model term(s) while estimating the model and then use the gof function again. The original error message from the ergm package will follow.")
          stop(e)
        } else {
          message("Error: Could not simulate any networks. Original error message follows.")
          stop(e)
        }
      }, 
      finally = {}
    )
  }
  
  # check basis network(s)
  if (verbose == TRUE) {
    if (l$time.steps == 1) {
      message("One network from which simulations are drawn was provided.\n")
    } else {
      message(paste(l$time.steps, "networks from which simulations are",
                    "drawn were provided.\n"))
    }
  }
  
  # unpack nested lists of simulations and target statistics
  simulations <- list()
  for (i in 1:length(sim)) {
    for (j in 1:length(sim[[i]])) {
      simulations[[length(simulations) + 1]] <- sim[[i]][[j]]
    }
  }
  rm(sim)
  
  # if NA in target networks, put them in the base network, too, and vice-versa
  if (length(l$networks) == length(target)) {
    for (i in 1:l$time.steps) {
      l$networks[[i]] <- as.matrix(l$networks[[i]])
      target[[i]] <- as.matrix(target[[i]])
      if (nrow(target[[i]]) != nrow(l$networks[[i]])) {
        stop(paste0("Dimensions of observed network and target do not match ", 
                    "at t = ", i, ": observed network has ", nrow(l$networks[[i]]), 
                    " rows while target has ", nrow(target[[i]]), " rows."))
      }
      if (ncol(target[[i]]) != ncol(l$networks[[i]])) {
        stop(paste0("Dimensions of observed network and target do not match ", 
                    "at t = ", i, ": observed network has ", ncol(l$networks[[i]]), 
                    " columns while target has ", ncol(target[[i]]), " columns."))
      }
      l$networks[[i]][is.na(as.matrix(target[[i]]))] <- NA
      l$networks[[i]] <- network::network(l$networks[[i]], 
                                          directed = l$directed, bipartite = l$bipartite)
      target[[i]][is.na(as.matrix(l$networks[[i]]))] <- NA
      target[[i]] <- network::network(target[[i]], directed = l$directed, 
                                      bipartite = l$bipartite)
    }
  }
  
  # data preparation
  sptypes <- c("dgCMatrix", "dgTMatrix", "dsCMatrix", "dsTMatrix", "dgeMatrix")
  
  directed <- logical()
  twomode <- logical()
  for (i in 1:length(target)) {
    if (is.network(target[[i]])) {
      directed[i] <- is.directed(target[[i]])
      twomode[i] <- is.bipartite(target[[i]])
      target[[i]] <- Matrix(as.matrix(target[[i]]))
    } else if (is.matrix(target[[i]])) {
      directed[i] <- is.mat.directed(target[[i]])
      twomode[i] <- !is.mat.onemode(target[[i]])
      target[[i]] <- Matrix(target[[i]])
    } else if (class(target[[i]])[1] %in% sptypes) {
      # OK
      directed[i] <- is.mat.directed(target[[i]])
      twomode[i] <- !is.mat.onemode(target[[i]])
    }
  }
  
  simulations <- lapply(simulations, function(x) Matrix(as.matrix(x)))
  goflist <- createGOF(simulations = simulations, target = target, 
                       statistics = statistics, parallel = parallel, ncpus = ncpus, cl = cl, 
                       verbose = verbose, ... = ...)
  return(goflist)
}

#' @rdname gof
#' @importFrom network network is.network is.bipartite
#' @export
setMethod("gof", signature = className("btergm", "btergm"), 
          definition = gof.btergm)

#' @rdname gof
#' @export
setMethod("gof", signature = className("ergm", "ergm"), 
          definition = gof.btergm)

#' @rdname gof
#' @export
setMethod("gof", signature = className("mtergm", "btergm"), 
          definition = gof.btergm)

#' @rdname gof
#' @export
setMethod("gof", signature = className("tbergm", "btergm"), 
          definition = gof.btergm)

#' @noRd
gof.sienaFit <- function(object,
                         period = NULL,
                         parallel = c("no", "multicore", "snow"),
                         ncpus = 1,
                         cl = NULL,
                         structzero = 10,
                         statistics = c(esp,
                                        deg,
                                        ideg,
                                        geodesic,
                                        rocpr,
                                        walktrap.modularity), 
                         groupName = object$f$groupNames[[1]],
                         varName = NULL, 
                         outofsample = FALSE,
                         sienaData = NULL,
                         sienaEffects = NULL, 
                         nsim = NULL,
                         verbose = TRUE,
                         ...) {
  
  if (!requireNamespace("RSiena", quietly = TRUE)) {
    stop("gof.sienaFit requires the 'RSiena' package to be installed.\n",
         "To do this, enter 'install.packages(\"RSiena\")'.")
  }
  
  # check correct specification of all data
  if (outofsample == FALSE) {
    if (is.null(varName)) {
      varName <- object$f$depNames[[1]]
    }
    if (is.null(object$sims) || any(sapply(object$sims, is.null))) {
      stop(paste("The object does not contain any simulations. Please", 
                 "re-estimate the model again with argument 'returnDeps = TRUE'."))
    }
    if (verbose == TRUE) {
      if (!is.null(sienaData)) {
        message(paste("'sienaData' argument is ignored because", 
                      "outofsample = FALSE."))
      }
      if (!is.null(sienaEffects)) {
        message(paste("'sienaEffects' argument is ignored because", 
                      "outofsample = FALSE."))
      }
      if (!is.null(nsim)) {
        message(paste("'nsim' argument is ignored because", 
                      "outofsample = FALSE."))
      }
    }
  } else {  # check if all data are there for out-of-sample prediction
    if (is.null(sienaData) || !"siena" %in% class(sienaData)) {
      stop(paste("For out-of-sample prediction, the 'sienaData'", 
                 "argument should provide a 'siena' object (as created by the", 
                 "'sienaDependent' function), and this object should contain two", 
                 "networks: the base network to simulate from and the target network", 
                 "to be predicted."))
    }
    if (sienaData$observations > 2) {
      stop(paste("For out-of-sample prediction, the 'sienaData' argument", 
                 "should contain only two networks: the base network to simulate", 
                 "from and the target network to be predicted."))
    }
    if (is.null(sienaEffects) || !"sienaEffects" %in% class(sienaEffects)) {
      stop(paste("For out-of-sample GOF assessment, a 'sienaEffects' object", 
                 "must be provided. This object must contain the same effects that", 
                 "were used in the original estimation, and it must be based on", 
                 "the 'sienaData' object."))
    }
    if (is.null(varName)) {
      varName <- names(sienaData$depvars)[1]
    }
    if (is.null(nsim)) {
      if (verbose == TRUE) {
        message("'nsim' not given. Using a default of 100 simulations.")
      }
      nsim <- 100
    }
    if (verbose == TRUE) {
      if (!is.null(period)) {
        message("The 'period' argument is ignored because outofsample = TRUE.")
      }
    }
  }
  if (verbose == TRUE) {
    message(paste("Variable name:", varName))
    if (outofsample == FALSE) {
      message(paste("Group name:", groupName))
    }
  }
  
  # interpret period argument
  if (outofsample == FALSE) {
    if (is.null(period)) {
      base <- 1:(attr(object$f[[1]]$depvars[[1]], "netdims")[3] - 1)
    } else {
      if (period < 2) {
        stop("The 'period' argument should have a minimum value of 2.")
      } else if (period > object$observations) {
        stop(paste("The 'period' argument should not be larger than the number", 
                   "of observed time periods."))
      }
      base <- period - 1
    }
  } else {
    base <- 1
  }
  
  # out of sample: extract, inject and fix coefficients and simulate networks
  if (outofsample == TRUE) {
    coefs <- object$theta[object$BasicRateFunction == FALSE]
    numtheta.new <- length(sienaEffects$initialValue[sienaEffects$include])
    sienaEffects$initialValue[sienaEffects$include][2:numtheta.new] <- coefs
    sim_model <- RSiena::sienaAlgorithmCreate(projname = "sim_model", 
                                              cond = FALSE,
                                              useStdInits = FALSE,
                                              nsub = 0,
                                              n3 = nsim,
                                              simOnly = TRUE)
    sim_ans <- RSiena::siena07(sim_model,
                               data = sienaData,
                               effects = sienaEffects,
                               batch = TRUE,
                               silent = TRUE,
                               returnDeps = TRUE)
    object <- sim_ans  # replace object by new simulations
    if (is.null(groupName)) {
      groupName <- names(sim_ans$f)[1]
      if (verbose == TRUE) {
        message(paste("Group name:", groupName))
      }
    }
  }
  
  # save the target object in a list and remove/handle structural zeros
  target <- list()
  for (i in base) {
    if (outofsample == TRUE) {
      dv.temp <- sienaData$depvars[[varName]][, , i + 1]
    } else {
      dv.temp <- object$f[[groupName]]$depvars[[varName]][, , i + 1]
    }
    rownames(dv.temp) <- 1:nrow(dv.temp)
    colnames(dv.temp) <- 1:ncol(dv.temp)
    dv.temp <- suppressMessages(handleMissings(dv.temp, na = structzero, 
                                               method = "remove", logical = FALSE))
    target[[length(target) + 1]] <- dv.temp
  }
  
  # define function for extracting all simulations and making them conformable
  simSiena <- function(i, myobject, mybase, mygroup, mydvname, target) {
    l <- list(length(target))
    count <- 0
    for (j in mybase) {
      count <- count + 1
      s <- RSiena::sparseMatrixExtraction(i = i, obsData = myobject$f, 
                                          sims = myobject$sims, period = j, groupName = mygroup, 
                                          varName = mydvname)
      rownames(s) <- 1:nrow(s)
      colnames(s) <- 1:ncol(s)
      l[[count]] <- Matrix::Matrix(adjust(as.matrix(s), target[[count]]))
    }
    return(l)
  }
  
  # extract the simulations, possibly in parallel
  if (is.null(ncpus) || ncpus == 0) {
    ncpus <- 1
  }
  if (!parallel[1] %in% c("no", "multicore", "snow")) {
    parallel <- "no"
    warning("'parallel' argument not recognized. Using 'no' instead.")
  }
  if (parallel[1] == "snow") {
    if (is.null(cl)) {
      cl <- parallel::makeCluster(ncpus)
    }
    if (verbose == TRUE) {
      message(paste("Using snow parallelization with", ncpus, "cores."))
    }
    parallel::clusterEvalQ(cl, library("Matrix"))
    sim <- parallel::parLapply(cl,
                               1:length(object$sims),
                               simSiena, 
                               myobject = object,
                               mybase = base,
                               mygroup = groupName,
                               mydvname = varName,
                               target = target)
  } else if (parallel[1] == "multicore") {
    if (verbose == TRUE) {
      message(paste("Using multicore parallelization with", ncpus, "cores."))
    }
    sim <- parallel::mclapply(1:length(object$sims),
                              simSiena,
                              myobject = object, 
                              mybase = base,
                              mygroup = groupName,
                              mydvname = varName, 
                              target = target,
                              mc.cores = ncpus)
  } else {
    if (verbose == TRUE) {
      message(paste0("Parallelization is switched off. Extracting ", 
                     "simulations sequentially."))
    }
    sim <- lapply(1:length(object$sims),
                  simSiena,
                  myobject = object, 
                  mybase = base,
                  mygroup = groupName,
                  mydvname = varName,
                  target = target)
  }
  
  # unpack nested lists of simulations
  simulations <- list(length(base) * length(object$sims))
  count <- 1
  for (t in 1:length(base)) {
    for (i in 1:length(object$sims)) {
      simulations[[count]] <- sim[[i]][[t]]
      count <- count + 1
    }
  }
  rm(sim)
  
  # apply auxiliary functions and return list of comparisons
  goflist <- createGOF(simulations = simulations,
                       target = target, 
                       statistics = statistics,
                       parallel = parallel,
                       ncpus = ncpus,
                       cl = cl, 
                       verbose = verbose,
                       ... = ...)
  return(goflist)
}

#' @rdname gof
#' @importFrom parallel makeCluster clusterEvalQ parLapply mclapply
#' @export
setMethod("gof", signature = className("sienaFit", "RSiena"), 
          definition = gof.sienaFit)

#' @noRd
gof.network <- function(object, covariates, coef, target = NULL, 
                        nsim = 100, mcmc = FALSE, MCMC.interval = 1000, MCMC.burnin = 10000, 
                        parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL, 
                        statistics = c(dsp, esp, deg, ideg, geodesic, rocpr, walktrap.modularity), 
                        verbose = TRUE, ...) {
  
  if (nsim < 2) {
    stop("The 'nsim' argument must be greater than 1.")
  }
  
  # check dependent network
  nw <- object
  if (is.network(nw)) {
    directed <- network::is.directed(nw)
    bipartite <- network::is.bipartite(nw)
  } else if (is.matrix(nw)) {
    directed <- !isSymmetric(nw)
    if (nrow(nw) == ncol(nw)) {
      bipartite <- FALSE
    } else {
      bipartite <- TRUE
    }
    nw <- network::network(nw, bipartite = bipartite, directed = directed)
  } else {
    stop("'object' must be a network object or a matrix.")
  }
  time.steps <- 1
  num.vertices <- nrow(as.matrix(nw))
  
  # check and rearrange target network(s)
  if (is.null(target)) {
    if (verbose == TRUE) {
      message(paste("\nNo 'target' network(s) provided. Using networks on the",
                    "left-hand side of the model formula as observed networks.\n"))
    }
    target <- nw
  } else if (is.network(target) || is.matrix(target)) {
    # do nothing
    if (verbose == TRUE) {
      message("\nOne observed ('target') network was provided.\n")
    }
  } else if ("list" %in% class(target)) {
    if (verbose == TRUE) {
      message(paste("\n", length(target), "observed ('target') networks were",
                    "provided. Using the first network\n"))
    }
    target <- target[[1]]
    if (!is.matrix(target) && !is.network(target)) {
      stop("First target network was not a matrix or network object.")
    }
  } else {
    stop("'target' must be a network, matrix, or list of matrices or networks.")
  }
  
  # check predictors and assemble formula
  if (!"list" %in% class(covariates)) {
    stop("Covariates must be provided as a list of matrices.")
  }
  numcov <- length(covariates)
  if (numcov + 1 != length(coef)) {
    stop(paste("The 'coef' vector ought to have a coefficient for edges", 
               "plus the same number of coefficients as there are covariates.", 
               "Right now, there are", length(coef), "coefficients and", numcov, 
               "covariates."))
  }
  rhs <- "edges"
  for (i in 1:numcov) {
    if (!is.network(covariates[[i]]) && !is.matrix(covariates[[i]])) {
      stop(paste("Covariate", i, "is not a matrix or network object."))
    }
    if (nrow(as.matrix(covariates[[i]])) != nrow(as.matrix(nw))) {
      stop(paste("Number of row nodes of covariate", i, "is not compatible."))
    }
    if (ncol(as.matrix(covariates[[i]])) != ncol(as.matrix(nw))) {
      stop(paste("Number of column nodes of covariate", i, 
                 "is not compatible."))
    }
    rhs <- paste(rhs, "+ edgecov(covariates[[", i, "]])")
  }
  form <- as.formula(paste("nw ~", rhs))
  
  # simulations for statnet-style and rocpr GOF
  message(paste("Simulating", nsim, 
                "networks from the following formula:\n", 
                gsub("\\s+", " ", paste(deparse(form), collapse = "")), "\n"))
  if (mcmc == TRUE) {  # TODO: IMPLEMENT PARALLEL PROCESSING HERE
    simulations <- ergm::simulate_formula(form,
                                          nsim = nsim,
                                          coef = coef,
                                          constraints = ~ .,
                                          control = ergm::control.simulate.formula(MCMC.interval = MCMC.interval,
                                                                                   MCMC.burnin = MCMC.burnin))
  } else {
    dat <- sapply(covariates, function(x) c(as.matrix(x)))
    dat <- cbind(rep(1, nrow(dat)), dat)
    prob <- plogis(coef %*% t(dat))
    simval <- t(sapply(prob, function(x) rbinom(nsim, 1, x)))
    simulations <- apply(simval, 2, function(x) Matrix(x, 
                                                       nrow = num.vertices, byrow = FALSE))
  }
  
  # if NA in target networks, put them in the base network, too, and vice-versa
  nw <- as.matrix(nw)
  nw[is.na(as.matrix(target))] <- NA
  nw <- network::network(nw, directed = directed, bipartite = bipartite)
  target <- as.matrix(target)
  target[is.na(as.matrix(nw))] <- NA
  target <- list(Matrix(target))
  
  # data preparation
  sptypes <- c("dgCMatrix", "dgTMatrix", "dsCMatrix", "dsTMatrix", "dgeMatrix")
  
  directed <- logical()
  twomode <- logical()
  for (i in 1:length(target)) {
    if (is.network(target[[i]])) {
      directed[i] <- is.directed(target[[i]])
      twomode[i] <- is.bipartite(target[[i]])
      target[[i]] <- Matrix(as.matrix(target[[i]]))
    } else if (is.matrix(target[[i]])) {
      directed[i] <- is.mat.directed(target[[i]])
      twomode[i] <- !is.mat.onemode(target[[i]])
      target[[i]] <- Matrix(target[[i]])
    } else if (class(target[[i]])[1] %in% sptypes) {
      # OK
      directed[i] <- is.mat.directed(target[[i]])
      twomode[i] <- !is.mat.onemode(target[[i]])
    }
  }
  goflist <- createGOF(simulations = simulations, target = target, 
                       statistics = statistics, parallel = parallel, ncpus = ncpus, cl = cl, 
                       verbose = verbose, ... = ...)
  return(goflist)
}

#' @rdname gof
#' @export
setMethod("gof", signature = className("network", "network"), 
          definition = gof.network)

#' @rdname gof
#' @export
setMethod("gof", signature = className("matrix", "base"), 
          definition = gof.network)