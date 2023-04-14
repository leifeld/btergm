#' Temporal Exponential Random Graph Models by Bootstrapped Pseudolikelihood
#'
#' Temporal Exponential Random Graph Models by Bootstrapped Pseudolikelihood.
#'
#' Temporal Exponential Random Graph Models (TERGM) estimated by maximum
#' pseudolikelihood with bootstrapped confidence intervals, Markov Chain Monte
#' Carlo maximum likelihood, or Bayesian estimation. Goodness of fit assessment
#' for ERGMs, TERGMs, and SAOMs. Micro-level interpretation of ERGMs and TERGMs.
#'
#' The \pkg{btergm} package implements TERGMs with MPLE and bootstrapped
#' confidence intervals (\code{\link{btergm}} function), MCMC MLE
#' (\code{\link{mtergm}} function), or Bayesian estimation (\code{\link{tbergm}}
#' function). Goodness of fit assessment for ERGMs, TERGMs, SAOMs, and dyadic
#' independence models is possible with the generic \code{\link{gof}} function
#' and its various methods defined here in the \pkg{btergm} package. New
#' networks can be simulated from TERGMs using the \code{\link{simulate.btergm}}
#' function. The package also implements micro-level interpretation for ERGMs
#' and TERGMs using the \code{\link{interpret}} function. Furthermore, the
#' package contains the \code{\link{chemnet}} and \code{\link{knecht}} (T)ERGM
#' datasets. To display citation information, type \code{citation("btergm")}.
#'
#' @author Philip Leifeld, Skyler J. Cranmer, Bruce A. Desmarais
#'
#' @references
#' Cranmer, Skyler J., Tobias Heinrich and Bruce A. Desmarais (2014):
#' Reciprocity and the Structural Determinants of the International Sanctions
#' Network. \emph{Social Networks} 36(1): 5-22.
#' \doi{10.1016/j.socnet.2013.01.001}.
#'
#' Desmarais, Bruce A. and Skyler J. Cranmer (2012): Statistical Mechanics of
#' Networks: Estimation and Uncertainty. \emph{Physica A} 391: 1865--1876.
#' \doi{10.1016/j.physa.2011.10.018}.
#'
#' Desmarais, Bruce A. and Skyler J. Cranmer (2010): Consistent Confidence
#' Intervals for Maximum Pseudolikelihood Estimators. \emph{Neural Information
#' Processing Systems 2010 Workshop on Computational Social Science and the
#' Wisdom of Crowds}.
#'
#' Leifeld, Philip, Skyler J. Cranmer and Bruce A. Desmarais (2018): Temporal
#' Exponential Random Graph Models with btergm: Estimation and Bootstrap
#' Confidence Intervals. \emph{Journal of Statistical Software} 83(6): 1--36.
#' \doi{10.18637/jss.v083.i06}.
#'
#' @name btergm-package
#' @docType package
NULL

#' Display version number and date when the package is loaded
#'
#' Display version number and date when the package is loaded.
#'
#' @importFrom utils packageDescription
#' @noRd
.onAttach <- function(libname, pkgname) {
  desc  <- packageDescription(pkgname, libname)
  packageStartupMessage(
    'Package:  btergm\n',
    'Version:  ', desc$Version, '\n',
    'Date:     ', desc$Date, '\n',
    'Authors:  Philip Leifeld (University of Essex)\n',
    '          Skyler J. Cranmer (The Ohio State University)\n',
    '          Bruce A. Desmarais (Pennsylvania State University)\n'
  )
}

#' Redefine S3 as S4 class for proper handling as part of the \code{btergm}
#' class.
#'
#' @noRd
setOldClass(c("boot", "boot"))

#' An S4 class to represent a fitted TERGM by bootstrapped MPLE
#'
#' An S4 class to represent a fitted TERGM by bootstrapped MPLE.
#'
#' \code{btergm} objects result from the estimation of a bootstrapped TERGM via
#' the \code{\link{btergm}} function. \code{btergm} objects contain the
#' coefficients, the bootstrapping samples of the coefficients, the number of
#' replications, the number of observations, the number of time steps, the
#' original formula, and the response, effects and weights objects that were fed
#' into the \code{glm} call for estimating the model.
#'
#' @slot coef Object of class \code{"numeric"}. The coefficients.
#' @slot boot Object of class \code{"matrix"}. The bootstrapping sample.
#' @slot R Object of class \code{"numeric"}. Number of replications.
#' @slot nobs Object of class \code{"numeric"}. Number of observations.
#' @slot time.steps Object of class \code{"numeric"}. Number of time steps.
#' @slot formula Object of class \code{"formula"}. The original model formula
#'   (without indices for the time steps).
#' @slot formula2 The revised formula with the object references after applying
#'   the \code{\link{tergmprepare}} function.
#' @slot response Object of class \code{"integer"}. The response variable.
#' @slot effects Object of class \code{"data.frame"}. The effects that went
#'   into the \code{glm} call.
#' @slot weights Object of class \code{"integer"}. The weights of the
#'   observations.
#' @slot auto.adjust Object of class \code{"logical"}. Indicates whether
#'   automatic adjustment of dimensions was done before estimation.
#' @slot offset Object of class \code{"logical"}. Indicates whether an offset
#'   matrix with structural zeros was used.
#' @slot directed Object of class \code{"logical"}. Are the dependent networks
#'   directed?
#' @slot bipartite Object of class \code{"logical"}. Are the dependent networks
#'   bipartite?
#' @slot nvertices Number of vertices.
#' @slot data The data after processing by the \code{\link{tergmprepare}}
#'   function.
#'
#' @family tergm-classes
#'
#' @export
setClass(Class = "btergm",
    representation = representation(
        coef = "numeric",
        boot = "boot",
        R = "numeric",
        nobs = "numeric",
        time.steps = "numeric",
        formula = "formula",
        formula2 = "character",
        response = "integer",
        effects = "data.frame",
        weights = "numeric",
        auto.adjust = "logical",
        offset = "logical",
        directed = "logical",
        bipartite = "logical",
        nvertices = "matrix",
        data = "list"
    ),
    validity = function(object) {
        if (!"numeric" %in% class(object@coef)) {
          stop("'coef' must be a 'numeric' object.")
        }
        if (!"boot" %in% class(object@boot)) {
          stop("'boot' must be a 'boot' object.")
        }
        if (!is.numeric(object@R)) {
          stop("'R' must be a numeric value of length 1.")
        }
        if (!is.numeric(object@nobs)) {
          stop("'nobs' must be a numeric value of length 1.")
        }
        if (!is.numeric(object@time.steps)) {
          stop("'time.steps' must be a numeric value of length 1.")
        }
        if (!"formula" %in% class(object@formula)) {
          stop("'formula' is not a 'formula' object.")
        }
        if (!is.integer(object@response)) {
          stop("'response' must consist of 'integer' values.")
        }
        if (!is.data.frame(object@effects)) {
          stop("'effects' must be a 'data.frame'.")
        }
        if (nrow(object@boot$t) != object@R) {
          error("The sample size does not correspond to the 'R' parameter.")
        }
        if (length(object@coef) != ncol(object@boot$t)) {
          stop("Number of terms differs between 'boot' and 'coef'")
        }
        if (length(object@response) != nrow(object@effects)) {
          stop("'response' and 'effects' must have the same length.")
        }
        if (!is.numeric(object@weights)) {
          stop("'weights' must consist of 'integer' or 'numeric' values.")
        }
        return(TRUE)
    }
)

#' Constructor for \linkS4class{btergm} objects
#'
#' Constructor for \linkS4class{btergm} objects.
#'
#' Create an S4 \linkS4class{btergm} object using this constructor function.
#'
#' @param coef Object of class \code{"numeric"}. The coefficients.
#' @param boot Object of class \code{"matrix"}. The bootstrapping sample.
#' @param R Object of class \code{"numeric"}. Number of replications.
#' @param nobs Object of class \code{"numeric"}. Number of observations.
#' @param time.steps Object of class \code{"numeric"}. Number of time steps.
#' @param formula Object of class \code{"formula"}. The original model formula
#'   (without indices for the time steps).
#' @param formula2 The revised formula with the object references after applying
#'   the \code{\link{tergmprepare}} function.
#' @param response Object of class \code{"integer"}. The response variable.
#' @param effects Object of class \code{"data.frame"}. The effects that went
#'   into the \code{glm} call.
#' @param weights Object of class \code{"integer"}. The weights of the
#'   observations.
#' @param auto.adjust Object of class \code{"logical"}. Indicates whether
#'   automatic adjustment of dimensions was done before estimation.
#' @param offset Object of class \code{"logical"}. Indicates whether an offset
#'   matrix with structural zeros was used.
#' @param directed Object of class \code{"logical"}. Are the dependent networks
#'   directed?
#' @param bipartite Object of class \code{"logical"}. Are the dependent networks
#'   bipartite?
#' @param nvertices Number of vertices.
#' @param data The data after processing by the \code{\link{tergmprepare}}
#'   function.
#'
#' @author Philip Leifeld
#'
#' @family tergm-classes
#'
#' @importFrom methods new
createBtergm <- function(coef, boot, R, nobs, time.steps, formula,
    formula2, response, effects, weights, auto.adjust, offset,
    directed, bipartite, nvertices, data) {
  new("btergm", coef = coef, boot = boot, R = R, nobs = nobs,
      time.steps = time.steps, formula = formula, formula2 = formula2,
      response = response, effects = effects, weights = weights,
      auto.adjust = auto.adjust, offset = offset, directed = directed,
      bipartite = bipartite, nvertices = nvertices, data = data)
}

#' Show the coefficients of a \code{btergm} object
#'
#' Show the coefficients of a \code{btergm} object.
#'
#' @param object A \code{btergm} object.
#'
#' @rdname btergm-class
#'
#' @importFrom methods show
#' @export
setMethod(f = "show", signature = "btergm", definition = function(object) {
    message("MLE Coefficients:")
    print(object@coef)
  }
)

#' @describeIn btergm-class Return the coefficients of a \code{btergm} object.
#'
#' @param object A \code{btergm} object.
#' @param invlogit Apply inverse logit transformation to the estimates and/or
#'   confidence intervals? That is, \eqn{\frac{1}{1 + \exp(-x)}}, where \eqn{x}
#'   is the respective value.
#'
#' @export
setMethod(f = "coef", signature = "btergm", definition = function(object,
      invlogit = FALSE, ...) {
    if (invlogit == FALSE) {
      return(object@coef)
    } else {
      return(1 / (1 + exp(-object@coef)))
    }
  }
)

#' @describeIn btergm-class Return the number of observations saved in a
#'   \code{btergm} object.
#'
#' @param object A \code{btergm} object.
#'
#' @export
setMethod(f = "nobs", signature = "btergm", definition = function(object) {
    n <- object@nobs
    t <- object@time.steps
    rep <- object@R
    return(c("Number of time steps" = t, "Number of dyads" = n,
        "Bootstrap replications" = rep))
  }
)

# function which can extract a coefficient matrix with SEs and p values
#' @describeIn btergm-class Create a coefficient table from a \code{btergm}
#'   object
#'
#' Create a coefficient matrix with standard errors and p-values.
#'
#' This function can create a coefficient matrix with coefficients, standard
#' errors, z-scores, and p-values, based on a fitted \code{btergm} object.
#' If the argument \code{print = TRUE} is used, the matrix is printed to the R
#' console as a formatted coefficient matrix with significance stars instead.
#' Note that confidence intervals are the preferred way of interpretation for
#' bootstrapped TERGMs; standard errors are only accurate if the bootstrapped
#' data are normally distributed, which is not always the case. Various methods
#' for checking for normality for each model term are available, for example
#' quantile-quantile plots (e.g., \code{qqnorm(x@boot$t[, 1])} for the first
#' model term in the \code{btergm} object called \code{x}).
#'
#' @param object A \code{btergm} object.
#' @param print Should the formatted coefficient table be printed to the R
#'   console along with significance stars (\code{print = TRUE}), or should the
#'   plain coefficient matrix be returned (\code{print = FALSE})?
#'
#' @export
btergm.se <- function(object, print = FALSE) {
  co <- object@coef
  sdev <- numeric()
  for (i in 1:ncol(object@boot$t)) {
    currentcol <- numeric()
    for (j in 1:nrow(object@boot$t)) {
      currentcol[j] <- (object@boot$t[j, i] - co[i])^2
    }
    sdev[i] <- sqrt(sum(currentcol) / length(currentcol))
  }
  zval <- (0 - apply(object@boot$t, 2, mean)) / sdev
  pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
  cmat <- cbind(co, sdev, zval, pval)
  colnames(cmat) <- c("Estimate", "Std.Err", "Z value", "Pr(>z)")
  warning(paste("Standard errors and p values may be misleading because the",
      "distribution of bootstrapped thetas may not be normal. Please rely on",
      "the confidence intervals instead or make sure the thetas are normally",
      "distributed (e.g., using qqnorm(object@boot$t[, 1]) etc."))
  if (print == TRUE) {
    printCoefmat(cmat)
  } else {
    return(cmat)
  }
}

#' @describeIn btergm-class Return the confidence intervals for estimates in a
#'   \code{btergm} object.
#'
#' @param object A \code{btergm} object.
#' @param parm Parameters (specified by integer position or character string).
#' @param level The significance level for computation of the confidence
#'   intervals. The default is \code{0.95} (that is, an alpha value of 0.05).
#'   Other common values include \code{0.999}, \code{0.99}, and \code{0.9}.
#' @param type Type of confidence interval, e.g., basic bootstrap interval
#'   (\code{type = "basic"}), percentile-based interval (\code{type = "perc"},
#'   which is the default option), or bias-adjusted and accelerated confidence
#'   interval (\code{type = "bca"}). All options from the \code{type} argument
#'   of the \link[boot]{boot.ci} function in the boot package can be used to
#'   generate confidence intervals.
#' @param invlogit Apply inverse logit transformation to the estimates and/or
#'   confidence intervals? That is, \eqn{\frac{1}{1 + \exp(-x)}}, where \eqn{x}
#'   is the respective value.
#' @param ... Further arguments for subroutines (not currently in use here).
#'
#' @export
setMethod(f = "confint", signature = "btergm", definition = function(object,
    parm, level = 0.95, type = "perc", invlogit = FALSE, ...) {
    cf <- coef(object, invlogit = invlogit)
    pnames <- names(cf)
    if (missing(parm)) {
      parm <- pnames
    } else if (is.numeric(parm)) {
      parm <- pnames[parm]
    }
    n.orig <- nrow(object@boot$t)
    object@boot$t <- object@boot$t[complete.cases(object@boot$t), ]
    n.ret <- nrow(object@boot$t)
    perc <- 100 * (n.orig - n.ret) / n.orig
    if (n.orig != n.ret) {
      warning(paste0("Too little variation in the model. ", n.orig - n.ret,
          " replications (", perc, "%) are dropped from CI estimation."))
    }
    if (invlogit == TRUE) {
      object@boot$t <- apply(object@boot$t, 1:2, function(x) 1 / (1 + exp(-x)))
      object@boot$t0 <- sapply(object@boot$t0, function(x) 1 / (1 + exp(-x)))
    }
    if (type == "perc") {
      type2 <- "percent"
    } else if (type == "norm") {
      type2 <- "normal"
    } else if (type == "basic") {
      type2 <- "basic"
    } else if (type == "stud") {
      type2 <- "student"
    } else if (type == "bca") {
      type2 <- "bca"
    } else {
      stop(paste("'type' not supported. Use 'perc', 'bca', 'norm', 'basic',",
          "or 'stud'."))
    }
    ci <- sapply(1:length(cf), function(x) {
        b <- boot::boot.ci(object@boot, conf = level, type = type, index = x)
        b[[type2]][4:5]
    })
    ci <- cbind(cf, colMeans(object@boot$t), t(ci))
    if ("numeric" %in% class(ci)) {
      ci.nam <- names(ci)
      ci <- matrix(ci, nrow = 1)
      colnames(ci) <- ci.nam
      rownames(ci) <- names(cf)
    }
    ci <- ci[parm, ]
    if (!is.matrix(ci)) {
      ci <- matrix(ci, ncol = 3)
      rownames(ci) <- parm
    }
    label1 <- paste0(100 * (1 - level) / 2, "%")
    label2 <- paste0(100 * (1 - (1 - level) / 2), "%")
    colnames(ci) <- c("Estimate", "Boot mean", label1, label2)
    return(ci)
  }
)

#' @describeIn btergm-class Return the number of time steps saved in a
#'   \code{btergm} object.
#'
#' @param object A \code{btergm} object.
#'
#' @export
timesteps.btergm <- function(object) {
  return(object@time.steps)
}

#' @describeIn btergm-class Summary of a fitted \code{btergm} object.
#'
#' @param object A \code{btergm} object.
#' @param level The significance level for computation of the confidence
#'   intervals. The default is \code{0.95} (that is, an alpha value of 0.05).
#'   Other common values include \code{0.999}, \code{0.99}, and \code{0.9}.
#' @param type Type of confidence interval, e.g., basic bootstrap interval
#'   (\code{type = "basic"}), percentile-based interval (\code{type = "perc"},
#'   which is the default option), or bias-adjusted and accelerated confidence
#'   interval (\code{type = "bca"}). All options from the \code{type} argument
#'   of the \link[boot]{boot.ci} function in the boot package can be used to
#'   generate confidence intervals.
#' @param invlogit Apply inverse logit transformation to the estimates and/or
#'   confidence intervals? That is, \eqn{\frac{1}{1 + \exp(-x)}}, where \eqn{x}
#'   is the respective value.
#' @param ... Further arguments to be passed through to the \code{confint}
#'   function.
#'
#' @export
setMethod(f = "summary", signature = "btergm", definition = function(object,
    level = 0.95, type = "perc", invlogit = FALSE, ...) {
    message(paste(rep("=", 26), collapse=""))
    message("Summary of model fit")
    message(paste(rep("=", 26), collapse=""))
    message(paste("\nFormula:  ", gsub("\\s+", " ",
        paste(deparse(object@formula), collapse = "")), "\n"))
    message(paste("Time steps:", object@time.steps, "\n"))
    message(paste("Bootstrapping sample size:", object@R, "\n"))

    message(paste0("Estimates and ", 100 * level, "% confidence intervals:"))
    cmat <- confint(object, level = level, type = type, invlogit = invlogit,
        ...)
    printCoefmat(cmat, cs.ind = 1, tst.ind = 3:4)
  }
)

#' Estimate a TERGM by MPLE with temporal bootstrapping
#'
#' Estimate a TERGM by MPLE with temporal bootstrapping.
#'
#' The \code{btergm} function computes temporal exponential random graph models
#' (TERGM) by bootstrapped pseudolikelihood, as described in Desmarais and
#' Cranmer (2012). It is faster than MCMC-MLE but only asymptotically unbiased
#' the longer the time series of networks because it uses temporal bootstrapping
#' to correct the standard errors.
#'
#' @param formula Formula for the TERGM. Model construction works like in the
#'   \pkg{ergm} package with the same model terms etc. (for a list of terms, see
#'   \code{help("\link[ergm]{ergm-terms}")}). The networks to be modeled on the
#'   left-hand side of the equation must be given either as a list of network
#'   objects with more recent networks last (i.e., chronological order) or as a
#'   list of matrices with more recent matrices at the end. \code{dyadcov} and
#'   \code{edgecov} terms accept time-independent covariates (as \code{network}
#'   or \code{matrix} objects) or time-varying covariates (as a list of networks
#'   or matrices with the same length as the list of networks to be modeled).
#' @param R Number of bootstrap replications. The higher the number of
#'   replications, the more accurate but also the slower is the estimation.
#' @param offset If \code{offset = TRUE} is set, a list of offset matrices
#'   (one for each time step) with structural zeros is handed over to the
#'   pseudolikelihood preparation routine. The offset matrices contain
#'   structural zeros where either the dependent networks or any of the
#'   covariates have missing nodes (if \code{auto.adjust = TRUE} is used). All
#'   matrices and network objects are inflated to the dimensions of the largest
#'   object, and the offset matrices inform the estimation preparation routine
#'   which dyads are constrained to be absent. After MPLE data preparation, the
#'   dyads with these structural zeros are removed before the GLM is estimated.
#'   If \code{offset = FALSE} is set (the default behavior), all nodes that are
#'   not present across all covariates and networks within a time step are
#'   removed completely from the respective object(s) before estimation begins.
#' @param returndata Return the processed input data instead of estimating and
#'   returning the model? In the \code{btergm} case, this will return a data
#'   frame with the dyads of the dependent variable/network and the change
#'   statistics for all covariates. In the \code{mtergm} case, this will return
#'   a list object with the blockdiagonal network object for the dependent
#'   variable and blockdiagonal matrices for all dyadic covariates and the
#'   offset matrix for the structural zeros.
#' @param parallel Use multiple cores in a computer or nodes in a cluster to
#'   speed up bootstrapping computations. The default value \code{"no"} means
#'   parallel computing is switched off. If \code{"multicore"} is used, the
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
#'   \code{cl} argument). The backend for the bootstrapping procedure is the
#'   \pkg{boot} package.
#' @param ncpus The number of CPU cores used for parallel computing (only if
#'   \code{parallel} is activated). If the number of cores should be detected
#'   automatically on the machine where the code is executed, one can set
#'   \code{ncpus = detectCores()} after loading the \pkg{parallel} package.
#'   On some HPC clusters, the number of available cores is saved as an
#'   environment variable; for example, if MOAB is used, the number of
#'   available cores can sometimes be accessed using
#'   \code{Sys.getenv("MOAB_PROCCOUNT")}, depending on the implementation.
#' @param cl An optional \pkg{parallel} or \pkg{snow} cluster for use if
#'   \code{parallel = "snow"}. If not supplied, a PSOCK cluster is created
#'   temporarily on the local machine.
#' @param control.ergm ergm controls for \code{\link[ergm]{ergmMPLE}} calls. See
#'   \code{\link[ergm]{control.ergm}} for details.
#' @param usefastglm Controls whether to use the \code{\link[fastglm]{fastglm}}
#'   estimation routine from the \pkg{fastglm} package with \code{method = 3}.
#'   Defaults to \code{FALSE} (and then uses
#'   \code{\link[speedglm:speedglm]{speedglm.wfit}} instead if available).
#' @param verbose Print details about data preprocessing and estimation
#'   settings.
#' @param ... Further arguments to be handed over to the
#'   \code{\link[boot]{boot}} function.
#'
#' @author Philip Leifeld, Skyler J. Cranmer, Bruce A. Desmarais
#'
#' @seealso \code{\link{mtergm}} \code{\link{tbergm}}
#'
#' @references Cranmer, Skyler J., Tobias Heinrich and Bruce A. Desmarais
#'   (2014): Reciprocity and the Structural Determinants of the International
#'   Sanctions Network. \emph{Social Networks} 36(1): 5-22.
#'   \doi{10.1016/j.socnet.2013.01.001}.
#'
#' Desmarais, Bruce A. and Skyler J. Cranmer (2012): Statistical Mechanics of
#'   Networks: Estimation and Uncertainty. \emph{Physica A} 391: 1865--1876.
#'   \doi{10.1016/j.physa.2011.10.018}.
#'
#' Desmarais, Bruce A. and Skyler J. Cranmer (2010): Consistent Confidence
#'   Intervals for Maximum Pseudolikelihood Estimators. \emph{Neural Information
#'   Processing Systems 2010 Workshop on Computational Social Science and the
#'   Wisdom of Crowds}.
#'
#' Leifeld, Philip, Skyler J. Cranmer and Bruce A. Desmarais (2017):
#'   Temporal Exponential Random Graph Models with btergm: Estimation and
#'   Bootstrap Confidence Intervals. \emph{Journal of Statistical Software}
#'   83(6): 1-36. \doi{10.18637/jss.v083.i06}.
#'
#' @examples
#' set.seed(5)
#'
#' networks <- list()
#' for (i in 1:10) {              # create 10 random networks with 10 actors
#'   mat <- matrix(rbinom(100, 1, .25), nrow = 10, ncol = 10)
#'   diag(mat) <- 0               # loops are excluded
#'   nw <- network::network(mat)  # create network object
#'   networks[[i]] <- nw          # add network to the list
#' }
#'
#' covariates <- list()
#' for (i in 1:10) {              # create 10 matrices as covariate
#'   mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
#'   covariates[[i]] <- mat       # add matrix to the list
#' }
#'
#' fit <- btergm(networks ~ edges + istar(2) + edgecov(covariates), R = 100)
#' summary(fit)                   # show estimation results
#'
#' # For examples with real data, see help("knecht") or help("alliances").
#'
#'
#' # Examples for parallel processing:
#'
#' # Some preliminaries:
#' # - "Forking" means running the code on multiple cores in the same
#' #   computer. It's fast but consumes a lot of memory because all
#' #   objects are copied for each node. It's also restricted to
#' #   cores within a physical computer, i.e. no distribution over a
#' #   network or cluster. Forking does not work on Windows systems.
#' # - "MPI" is a protocol for distributing computations over many
#' #   cores, often across multiple physical computers/nodes. MPI
#' #   is fast and can distribute the work across hundreds of nodes
#' #   (but remember that R can handle a maximum of 128 connections,
#' #   which includes file access and parallel connections). However,
#' #   it requires that the Rmpi package is installed and that an MPI
#' #   server is running (e.g., OpenMPI).
#' # - "PSOCK" is a TCP-based protocol. It can also distribute the
#' #   work to many cores across nodes (like MPI). The advantage of
#' #   PSOCK is that it can as well make use of multiple nodes within
#' #   the same node or desktop computer (as with forking) but without
#' #   consuming too much additional memory. However, the drawback is
#' #   that it is not as fast as MPI or forking.
#' # The following code provides examples for these three scenarios.
#'
#' # btergm works with clusters via the parallel package. That is, the
#' # user can create a cluster object (of type "PSOCK", "MPI", or
#' # "FORK") and supply it to the 'cl' argument of the 'btergm'
#' # function. If no cluster object is provided, btergm will try to
#' # create a temporary PSOCK cluster (if parallel = "snow") or it
#' # will use forking (if parallel = "multicore").
#'
#' \dontrun{
#' # To use a PSOCK cluster without providing an explicit cluster
#' # object:
#' require("parallel")
#' fit <- btergm(networks ~ edges + istar(2) + edgecov(covariates),
#'               R = 100, parallel = "snow", ncpus = 25)
#'
#' # Equivalently, a PSOCK cluster can be provided as follows:
#' require("parallel")
#' cores <- 25
#' cl <- makeCluster(cores, type = "PSOCK")
#' fit <- btergm(networks ~ edges + istar(2) + edgecov(covariates),
#'               R = 100, parallel = "snow", ncpus = cores, cl = cl)
#' stopCluster(cl)
#'
#' # Forking (without supplying a cluster object) can be used as
#' # follows.
#' require("parallel")
#' cores <- 25
#' fit <- btergm(networks ~ edges + istar(2) + edgecov(covariates),
#'               R = 100, parallel = "multicore", ncpus = cores)
#' stopCluster(cl)
#'
#' # Forking (by providing a cluster object) works as follows:
#' require("parallel")
#' cores <- 25
#' cl <- makeCluster(cores, type = "FORK")
#' fit <- btergm(networks ~ edges + istar(2) + edgecov(covariates),
#'               R = 100, parallel = "snow", ncpus = cores, cl = cl)
#' stopCluster(cl)
#'
#' # To use MPI, a cluster object MUST be created beforehand. In
#' # this example, a MOAB HPC server is used. It stores the number of
#' # available cores as a system option:
#' require("parallel")
#' cores <- as.numeric(Sys.getenv("MOAB_PROCCOUNT"))
#' cl <- makeCluster(cores, type = "MPI")
#' fit <- btergm(networks ~ edges + istar(2) + edgecov(covariates),
#'               R = 100, parallel = "snow", ncpus = cores, cl = cl)
#' stopCluster(cl)
#'
#' # In the following example, the Rmpi package is used to create a
#' # cluster. This may not work on all systems; consult your local
#' # support staff or the help files on your HPC server to find out how
#' # to create a cluster object on your system.
#'
#' # snow/Rmpi start-up
#' if (!is.loaded("mpi_initialize")) {
#'   library("Rmpi")
#' }
#' library(snow);
#'
#' mpirank <- mpi.comm.rank (0)
#' if (mpirank == 0) {
#'   invisible(makeMPIcluster())
#' } else {
#'   sink (file="/dev/null")
#'   invisible(slaveLoop (makeMPImaster()))
#'   mpi.finalize()
#'   q()
#' }
#' # End snow/Rmpi start-up
#'
#' cl <- getMPIcluster()
#'
#' fit <- btergm(networks ~ edges + istar(2) + edgecov(covariates),
#'               R = 100, parallel = "snow", ncpus = 25, cl = cl)
#' }
#'
#' @importFrom boot boot
#' @importFrom ergm control.ergm ergm.design ergm_model ergm.getnetwork ergm.pl ergmMPLE
#' @importFrom Matrix Matrix
#' @import stats
#' @export
btergm <- function(formula,
                   R = 500,
                   offset = FALSE,
                   returndata = FALSE,
                   parallel = c("no", "multicore", "snow"),
                   ncpus = 1,
                   cl = NULL,
                   control.ergm = NULL,
                   usefastglm = FALSE,
                   verbose = TRUE,
                   ...) {

  if (is.null(control.ergm)) {
    control.ergm <- ergm::control.ergm()
  }

  # call tergmprepare and integrate results in local environment
  l <- tergmprepare(formula = formula, offset = offset, verbose = verbose)
  for (i in 1:length(l$covnames)) {
    assign(l$covnames[i], l[[l$covnames[i]]])
  }
  assign("offsmat", l$offsmat)
  form <- as.formula(l$form)

  # check number of time steps
  if (l$time.steps == 1) {
    warning(paste("The confidence intervals and standard errors are",
        "meaningless because only one time step was provided."))
  }

  # verbose reporting
  if (verbose == TRUE && returndata == FALSE) {
    if (parallel[1] == "no") {
      parallel.msg <- "on a single computing core"
    } else if (parallel[1] == "multicore") {
      parallel.msg <- paste("using multicore forking on", ncpus, "cores")
    } else if (parallel[1] == "snow") {
      parallel.msg <- paste("using parallel processing on", ncpus, "cores")
    }
    if (offset == TRUE) {
      offset.msg <- "with offset matrix and "
    } else {
      offset.msg <- "with "
    }
    message("\nStarting pseudolikelihood estimation ", offset.msg,
          R, " bootstrapping replications ", parallel.msg, "...")
  } else if (verbose == TRUE && returndata == TRUE) {
    message("\nReturning data frame with change statistics.")
  }

  # create MPLE data structures
  if (offset == TRUE) {  # via structural zeros and an offset term
    # create the data for MPLE with structural zeros
    Y <- NULL  # dependent variable
    X <- NULL  # independent variables data frame
    W <- NULL  # weights for each observation
    O <- NULL  # offset term
    for (i in 1:length(l$networks)) {
      nw <- ergm::ergm.getnetwork(form)
      model <- ergm::ergm_model(form, nw, initialfit = TRUE)
      Clist.miss <- ergm::ergm.design(nw, verbose = FALSE)
      control.ergm$init <- c(rep(NA, length(l$rhs.terms) - 1), 1)
      pl <- ergm::ergm.pl(nw,
                          Clist.miss,
                          model,
                          theta.offset = c(rep(FALSE, length(l$rhs.terms) - 1), TRUE),
                          verbose = FALSE,
                          control = control.ergm)
      Y <- c(Y, pl$zy[pl$foffset == 0])
      X <- rbind(X, cbind(data.frame(pl$xmat[pl$foffset == 0, ], check.names = FALSE), i))
      W <- c(W, pl$wend[pl$foffset == 0])
      O <- c(O, pl$foffset[pl$foffset == 0])
    }
    term.names <- colnames(X)[-length(colnames(X))]
    term.names <- c(term.names, "time")
    colnames(X) <- term.names
  } else {  # by deleting structural zero observations per time step
    Y <- NULL
    X <- NULL
    W <- NULL
    O <- NULL  # will remain NULL and will be fed into GLM
    for (i in 1:length(l$networks)) {
      mpli <- ergm::ergmMPLE(form, control = control.ergm)
      Y <- c(Y, mpli$response)

      # fix different factor levels across time points
      if (i > 1 && ncol(X) != ncol(mpli$predictor) + 1) {
        cn.x <- colnames(X)[-ncol(X)]  # exclude last column "i"
        cn.i <- colnames(mpli$predictor)
        names.x <- cn.x[which(!cn.x %in% cn.i)]
        names.i <- cn.i[which(!cn.i %in% cn.x)]
        if (length(names.x) > 0) {
          for (nm in 1:length(names.x)) {
            mpli$predictor <- cbind(mpli$predictor, rep(0,
                nrow(mpli$predictor)))
            colnames(mpli$predictor)[ncol(mpli$predictor)] <- names.x[nm]
          }
        }
        if (length(names.i) > 0) {
          for (nm in 1:length(names.i)) {
            X <- cbind(X[, 1:(ncol(X) - 1)], rep(0, nrow(X)), X[, ncol(X)])
            colnames(X)[ncol(X) - 1] <- names.i[nm]
          }
        }
      }  # end of column fix

      X <- rbind(X, cbind(mpli$predictor, i))
      W <- c(W, mpli$weights)
    }
    term.names <- colnames(X)[-length(colnames(X))]
    term.names <- c(term.names, "time")
    X <- data.frame(X)
    colnames(X) <- term.names
  }

  # remove time variable for estimation
  unique.time.steps <- unique(X$time)
  x <- X[, -ncol(X)]
  x <- as.data.frame(x)  # in case there is only one column/model term
  time <- X$time
  rm(X)
  if (returndata == TRUE) {
    return(cbind(Y, x))
  }

  # create sparse matrix and compute start values for GLM
  if (ncol(x) == 1) {
    stop("At least two model terms must be provided to estimate a TERGM.")
  }

  if (isTRUE(usefastglm)) {
    if (requireNamespace("fastglm")) {
      xsparse <- NULL
      est <- fastglm::fastglm(y = Y,
                              x = as.matrix(x),
                              weights = W,
                              offset = O,
                              family = binomial(link = logit),
                              sparse = TRUE,
                              method = 3)

      startval <- est$coefficients

      # report NAs as warnings
      if (isTRUE(verbose) && any(is.na(startval))) {
        warning("The following coefficients yielded NA estimates: ",
                paste(names(startval)[which(is.na(startval))], collapse = ", "),
                ". Bootstrapping may not work with NA coefficients.")
      }

      estimate <- function(unique.time.steps,
                           bsi,
                           Yi = Y,
                           xsparsei = xsparse,
                           Wi = W,
                           Oi = O,
                           timei = time,
                           startvali = startval) {
        indic <- unlist(lapply(bsi, function(x) which(timei == x)))
        fastglm::fastglm(y = Yi[indic],
                         x = as.matrix(x)[indic, ],
                         weights = Wi[indic],
                         offset = Oi[indic],
                         family = binomial(link = logit),
                         method = 3)$coefficients
      }
    } else {
      stop("The 'fastglm' package was not found.")
    }
  } else {
    if (requireNamespace("speedglm", quietly = TRUE)) {
      xsparse <- Matrix(as.matrix(x), sparse = TRUE)
      est <- speedglm::speedglm.wfit(y = Y,
                                     X = xsparse,
                                     weights = W,
                                     offset = O,
                                     family = binomial(link = logit),
                                     sparse = TRUE)

      startval <- coef(est)

      # report NAs as warnings
      if (isTRUE(verbose) && any(is.na(startval))) {
        warning("The following coefficients yielded NA estimates: ",
                paste(names(startval)[which(is.na(startval))], collapse = ", "),
                ". Bootstrapping may not work with NA coefficients.")
      }

      # define function for bootstrapping and estimation
      estimate <- function(unique.time.steps, bsi, Yi = Y, xsparsei = xsparse,
                           Wi = W, Oi = O, timei = time, startvali = startval) {
        indic <- unlist(lapply(bsi, function(x) which(timei == x)))
        tryCatch(
          expr = {
            coef(speedglm::speedglm.wfit(y = Yi[indic],
                                         X = xsparsei[indic, ],
                                         weights = Wi[indic],
                                         offset = Oi[indic],
                                         family = binomial(link = logit),
                                         sparse = TRUE))
          },
          error = function(e) {
            # when fitted probabilities of 0 or 1 occur or when the algorithm does
            # not converge, use glm because it only throws a warning, not an error
            coef(glm.fit(y = Yi[indic], x = as.matrix(x)[indic, ],
                         weights = Wi[indic], offset = Oi[indic],
                         family = binomial(link = logit)))
          },
          warning = function(w) {
            warning(w)
          },
          finally = {}
        )
      }
    } else {
      warning("The 'speedglm' package was not found. Using the 'glm' function.")
      startval <- coef(glm.fit(y = Y,
                               x = as.matrix(x),
                               weights = W,
                               offset = O,
                               family = binomial(link = logit)))

      estimate <- function(unique.time.steps,
                           bsi,
                           Yi = Y,
                           xsparsei = xsparse,
                           Wi = W,
                           Oi = O,
                           timei = time,
                           startvali = startval) {
        indic <- unlist(lapply(bsi, function(x) which(timei == x)))
        coef(glm.fit(y = Yi[indic], x = as.matrix(x)[indic, ],
                     weights = Wi[indic], offset = Oi[indic],
                     family = binomial(link = logit)))
      }
    }
  }

  coefs <- boot(unique.time.steps, estimate, R = R, Yi = Y, xsparsei = xsparse,
                Wi = W, Oi = O, timei = time, startvali = startval,
                parallel = parallel, ncpus = ncpus, cl = cl, ...)

  if (coefs$t[1, 1] == "glm.fit: algorithm did not converge") {
    warning(paste("Algorithm did not converge. There might be a collinearity ",
        "between predictors and/or dependent networks at one or more time",
        "steps."))
  }
  if (sum(is.na(coefs$t)) > 0) {
    warning("NAs generated during bootstrap. This may be an indication of collinearity between predictors and/or dependent networks at one or more time steps.")
  }

  # create and return btergm object
  colnames(coefs$t) <- term.names[1:(length(term.names) - 1)]
  names(startval) <- colnames(coefs$t)
  data <- list()
  for (i in 1:length(l$covnames)) {
    data[[l$covnames[i]]] <- l[[l$covnames[i]]]
  }
  data$offsmat <- l$offsmat

  if (isTRUE(l$bipartite)) {
    nobs <- sum(sapply(data$offsmat, function(x) {
      length(x[x == 0])
    }))
  } else {
    nobs <- sum(sapply(data$offsmat, function(x) {
      diag(x) <- 1
      length(x[x == 0])
    }))
  }

  btergm.object <- createBtergm(startval, coefs, R, nobs, l$time.steps,
      formula, l$form, Y, x, W, l$auto.adjust, offset, l$directed, l$bipartite,
      nvertices = l$nvertices, data)
  if (verbose == TRUE) {
    message("Done.")
  }
  return(btergm.object)
}

#' Simulate Networks from a \code{btergm} Object
#'
#' Simulate networks from a \code{btergm} object using MCMC sampler.
#'
#' The \code{simulate.btergm} function is a wrapper for the
#' \code{\link[ergm]{simulate_formula}} function in the \pkg{ergm} package (see
#' \code{help("simulate.ergm")}). It can be used to simulate new networks from a
#' \code{btergm} object. The \code{index} argument specifies from which of the
#' original networks the new network(s) should be simulated. For example, if
#' \code{object} is an estimation based on cosponsorship networks from the 99th
#' to the 107th Congress (as in Desmarais and Cranmer 2012), and the
#' cosponsorship network in the 108th Congress should be predicted using the
#' \code{simulate.btergm} function, then the argument \code{index = 9} should be
#' passed to the function because the network should be based on the 9th network
#' in the list (that is, the latest network, which is the cosponsorship network
#' for the 107th Congress). Note that all relevant objects (the networks and the
#' covariates) must be present in the workspace (as was the case during the
#' estimation of the model).
#'
#' @param object A \code{btergm} or \code{mtergm} object, resulting from a call
#'   of the \code{\link{btergm}} or \code{\link{mtergm}} function.
#' @param nsim The number of networks to be simulated. Note that for values
#'   greater than one, a \code{network.list} object is returned, which can be
#'   indexed just like a \code{list} object, for example \code{mynetworks[[1]]}
#'   for the first simulated network in the object \code{mynetworks}.
#' @param seed Random number integer seed. See \link[base:Random]{set.seed}.
#' @param formula A model formula from which the new network(s) should be
#'   simulated. By default, the formula is taken from the \code{btergm} object.
#' @param index Index of the network from which the new network(s) should be
#'   simulated. The index refers to the list of response networks on the
#'   left-hand side of the model formula. Note that more recent networks are
#'   located at the end of the list. By default, the last (= most recent)
#'   network is used.
#' @param coef A vector of parameter estimates. By default, the coefficients are
#'   extracted from the given \code{btergm} object.
#' @param verbose Print additional details while running the simulations?
#' @param ... Further arguments are handed over to the
#'   \code{\link[ergm]{simulate_formula}} function in the \pkg{ergm} package.
#'
#' @references
#' Desmarais, Bruce A. and Skyler J. Cranmer (2012): Statistical Mechanics of
#' Networks: Estimation and Uncertainty. \emph{Physica A} 391: 1865--1876.
#' \doi{10.1016/j.physa.2011.10.018}.
#'
#' Leifeld, Philip, Skyler J. Cranmer and Bruce A. Desmarais (2018): Temporal
#' Exponential Random Graph Models with btergm: Estimation and Bootstrap
#' Confidence Intervals. \emph{Journal of Statistical Software} 83(6): 1--36.
#' \doi{10.18637/jss.v083.i06}.
#'
#' @examples
#' \dontrun{
#' # fit a TERGM to some toy data
#' library("network")
#' set.seed(5)
#' networks <- list()
#' for(i in 1:10){            # create 10 random networks with 10 actors
#'   mat <- matrix(rbinom(100, 1, .25), nrow = 10, ncol = 10)
#'   diag(mat) <- 0           # loops are excluded
#'   nw <- network(mat)       # create network object
#'   networks[[i]] <- nw      # add network to the list
#' }
#' covariates <- list()
#' for (i in 1:10) {          # create 10 matrices as covariate
#'   mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
#'   covariates[[i]] <- mat   # add matrix to the list
#' }
#' fit <- btergm(networks ~ edges + istar(2) +
#'                 edgecov(covariates), R = 100)
#'
#' # simulate 12 new networks from the last (= 10th) time step
#' sim1 <- simulate(fit, nsim = 12)
#'
#' # simulate 1 network from the first time step
#' sim2 <- simulate(fit, index = 1)
#'
#' # simulate network from t = 5 with larger covariate coefficient
#' coefs <- coef(fit)
#' coefs["edgecov.covariates[[i]]"] <- 0.5
#' sim3 <- simulate(fit, index = 5, coef = coefs)
#' }
#'
#' @importFrom ergm simulate_formula
#' @importFrom network is.network
#' @importFrom stats simulate
#' @export
simulate.btergm <- function(object, nsim = 1, seed = NULL, index = NULL,
    formula = getformula(object), coef = object@coef, verbose = TRUE, ...) {

  # retrieve data and integrate results locally
  for (i in 1:length(object@data)) {
    assign(names(object@data)[i], object@data[[i]])
  }
  form <- as.formula(object@formula2, env = environment())

  # check and correct index argument
  if (is.null(index)) {
    index <- object@time.steps
    if (verbose == TRUE) {
      message("\nNo index provided. Simulating from the last time step.")
    }
  } else if (!"numeric" %in% class(index)) {
    stop(paste("The 'index' argument must contain a numeric time point from",
        "which to simulate new networks."))
  } else if (index > object@time.steps) {
    index <- object@time.steps
    message(paste("Index larger than the number of time steps. Simulating",
        "from the last time step."))
  }
  i <- index

  # print formula from which networks are simulated
  if (verbose == TRUE) {
    f.i <- gsub("\\[\\[i\\]\\]", paste0("[[", index, "]]"),
        paste(deparse(form), collapse = ""))
    f.i <- gsub("\\s+", " ", f.i)
    f.i <- gsub("^networks", deparse(object@formula[[2]]), f.i)
    message(paste("Simulating", nsim, "networks from the following formula:\n",
        f.i, "\n"))
  }

  # simulate
  if (object@offset == TRUE) {
    coef <- c(coef, -Inf)
  }
  s <- ergm::simulate_formula(form, nsim = nsim, seed = seed, coef = coef,
      verbose = verbose, ...)
  if ("btergm" %in% class(object)) {
    return(s)
  } else if ("mtergm" %in% class(object)) {
    if (index == 1) {
      r1 <- 1
      c1 <- 1
    } else {
      r1 <- sum(object@nvertices[1, 1:(index - 1)]) + 1
      c1 <- sum(object@nvertices[2, 1:(index - 1)]) + 1
    }
    r2 <- sum(object@nvertices[1, 1:index])
    c2 <- sum(object@nvertices[2, 1:index])

    if (is.network(s)) {
      s <- list(s)
    }
    s <- lapply(s, function(sim) network(as.matrix(sim)[r1:r2, c1:c2],
        bipartite = object@bipartite, directed = object@directed))
    if (length(s) == 1) {
      s <- s[[1]]
    }
    return(s)
  }
}

#' Simulate Networks from an \code{mtergm} Object
#'
#' Simulate networks from an \code{mtergm} object using MCMC sampler.
#'
#' @rdname simulate.btergm
#' @export
simulate.mtergm <- simulate.btergm  # create a copy for mtergm objects
