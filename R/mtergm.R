#' Redefine S3 as S4 class for proper handling as part of the \code{mtergm}
#' class.
#' 
#' @noRd
setOldClass(c("ergm", "ergm"))

#' An S4 Class to represent a fitted TERGM by MCMC-MLE
#' 
#' An S4 class to represent a fitted TERGM by MCMC-MLE.
#' 
#' \code{mtergm} objects result from MCMC-MLE-based estimation of a TERGM via
#' the \code{\link{mtergm}} function. They contain the coefficients, standard
#' errors, and p-values, among other details.
#' 
#' @slot coef Object of class \code{"numeric"}. The coefficients.
#' @slot se Object of class \code{"numeric"}. The standard errors.
#' @slot pval Object of class \code{"numeric"}. The p-values.
#' @slot nobs Object of class \code{"numeric"}. Number of observations.
#' @slot time.steps Object of class \code{"numeric"}. Number of time steps.
#' @slot formula Object of class \code{"formula"}. The original model formula
#'   (without indices for the time steps).
#' @slot formula2 The revised formula with the object references after applying
#'   the \code{\link{tergmprepare}} function.
#' @slot auto.adjust Object of class \code{"logical"}. Indicates whether
#'   automatic adjustment of dimensions was done before estimation.
#' @slot offset Object of class \code{"logical"}. Indicates whether an offset
#'   matrix with structural zeros was used.
#' @slot directed Object of class \code{"logical"}. Are the dependent networks
#'   directed?
#' @slot bipartite Object of class \code{"logical"}. Are the dependent networks
#'   bipartite?
#' @slot estimate Estimate: either MLE or MPLE.
#' @slot loglik Log likelihood of the MLE.
#' @slot aic Akaike's Information Criterion.
#' @slot bic Bayesian Information Criterion.
#' @slot ergm The original \code{ergm} object as estimated by the
#'   \code{\link[ergm]{ergm}} function in the \pkg{ergm} package.
#' @slot nvertices Number of vertices.
#' @slot data The data after processing by the \code{\link{tergmprepare}}
#'   function.
#'
#' @author Philip Leifeld
#' 
#' @family tergm-classes
#' 
#' @export
setClass(Class = "mtergm", 
    slots = c(
        coef = "numeric", 
        se = "numeric", 
        pval = "numeric", 
        nobs = "numeric", 
        time.steps = "numeric",
        formula = "formula",
        formula2 = "character", 
        auto.adjust = "logical", 
        offset = "logical", 
        directed = "logical", 
        bipartite = "logical", 
        estimate = "character", 
        loglik = "numeric", 
        aic = "numeric", 
        bic = "numeric", 
        ergm = "ergm", 
        nvertices = "matrix", 
        data = "list"
    ), 
    validity = function(object) {
        if (!"numeric" %in% class(object@coef)) {
          stop("'coef' must be a 'numeric' vector.")
        }
        if (!"numeric" %in% class(object@se)) {
          stop("'se' must be a 'numeric' vector.")
        }
        if (!"numeric" %in% class(object@pval)) {
          stop("'pval' must be a 'numeric' vector.")
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
        if (length(object@coef) != length(object@se)) {
          stop("Number of terms differs between 'coef' and 'se'")
        }
        if (length(object@coef) != length(object@pval)) {
          stop("Number of terms differs between 'coef' and 'pval'")
        }
        if (length(object@loglik) > 1) {
          stop("'loglik' must be a numeric value of length 1.")
        }
        if (length(object@aic) > 1) {
          stop("'aic' must be a numeric value of length 1.")
        }
        if (length(object@bic) > 1) {
          stop("'bic' must be a numeric value of length 1.")
        }
        return(TRUE)
    }
)


#' Constructor for \linkS4class{mtergm} objects
#'
#' Constructor for \linkS4class{mtergm} objects.
#'
#' Create an S4 \linkS4class{mtergm} object using this constructor function.
#'
#' @param se Standard errors.
#' @param pval The p-values.
#' @param estimate Estimate: either MCMC MLE or MPLE.
#' @param loglik Log likelihood of the MLE.
#' @param aic Akaike's Information Criterion.
#' @param bic Bayesian Information Criterion.
#' @param ergm The original \code{ergm} object as estimated by the
#'   \code{\link[ergm]{ergm}} function in the \pkg{ergm} package.
#' @inheritParams createBtergm
#'
#' @author Philip Leifeld
#' 
#' @family tergm-classes
createMtergm <- function(coef, se, pval, nobs, time.steps, formula, formula2, 
    auto.adjust, offset, directed, bipartite, estimate, loglik, aic, bic, 
    ergm, nvertices, data) {
  new("mtergm", coef = coef, se = se, pval = pval, nobs = nobs, 
      time.steps = time.steps, formula = formula, formula2 = formula2, 
      auto.adjust = auto.adjust, offset = offset, directed = directed, 
      bipartite = bipartite, estimate = estimate, loglik = loglik, aic = aic, 
      bic = bic, ergm = ergm, nvertices = nvertices, data = data)
}

#' @describeIn mtergm-class Show the coefficients of an \code{mtergm} object.
#' 
#' @param object An \code{mtergm} object.
#' 
#' @export
setMethod(f = "show", signature = "mtergm", definition = function(object) {
    message("MLE Coefficients:")
    print(object@coef)
  }
)

#' @describeIn mtergm-class Return the coefficients of an \code{mtergm} object.
#' 
#' @param object An \code{mtergm} object.
#' @param invlogit Apply inverse logit transformation to the estimates and/or
#'   confidence intervals? That is, \eqn{\frac{1}{1 + \exp(-x)}}, where \eqn{x}
#'   is the respective value.
#' 
#' @export
setMethod(f = "coef", signature = "mtergm", definition = function(object, 
      invlogit = FALSE, ...) {
    if (invlogit == FALSE) {
      return(object@coef)
    } else {
      return(1 / (1 + exp(-object@coef)))
    }
  }
)

#' @describeIn mtergm-class Return the coefficients of an \code{mtergm} object.
#' 
#' @param object An \code{mtergm} object.
#' 
#' @export
setMethod(f = "nobs", signature = "mtergm", definition = function(object) {
    n <- object@nobs
    t <- object@time.steps
    return(c("Number of time steps" = t, "Number of observations" = n))
  }
)

#' @describeIn mtergm-class Return the number of time steps saved in an
#'   \code{mtergm} object.
#' 
#' @param object An \code{mtergm} object.
#' 
#' @export
timesteps.mtergm <- function(object) {
  return(object@time.steps)
}

#' @describeIn mtergm-class Return the coefficients of an \code{mtergm} object.
#' 
#' @param object An \code{mtergm} object.
#' @param ... Currently not in use.
#' 
#' @export
setMethod(f = "summary", signature = "mtergm", definition = function(object, 
    ...) {
    message(paste(rep("=", 26), collapse = ""))
    message("Summary of model fit")
    message(paste(rep("=", 26), collapse = ""))
    message(paste("\nFormula:  ", gsub("\\s+", " ", 
        paste(deparse(object@formula), collapse = "")), "\n"))
    message(paste("Time steps:", object@time.steps, "\n"))
    
    message("Monte Carlo MLE Results:")
    cmat <- cbind(object@coef, object@se, object@coef / object@se, object@pval)
    colnames(cmat) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    printCoefmat(cmat, cs.ind = 1:2, tst.ind = 3)
  }
)

#' Estimate a TERGM by MCMC-MLE
#'
#' Estimate a TERGM by Markov Chain Monte Carlo Maximum Likelihood Estimation
#'
#' The \code{mtergm} function computes TERGMs by MCMC MLE (or MPLE with
#' uncorrected standard errors) via blockdiagonal matrices and structural zeros.
#' It acts as a wrapper for the \pkg{ergm} package. The \code{btergm} function
#' is faster than the \code{mtergm} function but is only asymptotically unbiased
#' the longer the time series. The \code{mtergm} function yields unbiased
#' estimates and standard errors but may suffer from degeneracy if the model is
#' not specified in good keeping with the true data-generating process.
#'
#' @param constraints Constraints of the ERGM. See \code{\link[ergm]{ergm}} for
#'   details.
#' @param ... Further arguments to be handed over to the
#'   \code{\link[ergm]{ergm}} function.
#' @inheritParams btergm
#'
#' @author Philip Leifeld, Skyler J. Cranmer, Bruce A. Desmarais
#'
#' @seealso \code{\link{btergm}} \code{\link{tbergm}}
#'
#' @references Leifeld, Philip, Skyler J. Cranmer and Bruce A. Desmarais (2017):
#'   Temporal Exponential Random Graph Models with btergm: Estimation and
#'   Bootstrap Confidence Intervals. \emph{Journal of Statistical Software}
#'   83(6): 1-36. \doi{10.18637/jss.v083.i06}.
#'
#' @examples
#' library("network")
#' set.seed(5)
#' 
#' networks <- list()
#' for (i in 1:10) {          # create 10 random networks with 10 actors
#'   mat <- matrix(rbinom(100, 1, .25), nrow = 10, ncol = 10)
#'   diag(mat) <- 0           # loops are excluded
#'   nw <- network(mat)       # create network object
#'   networks[[i]] <- nw      # add network to the list
#' }
#' 
#' covariates <- list()
#' for (i in 1:10) {          # create 10 matrices as covariate
#'   mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
#'   covariates[[i]] <- mat   # add matrix to the list
#' }
#' 
#' \dontrun{
#' fit2 <- mtergm(networks ~ edges + istar(2) + edgecov(covariates))
#' summary(fit2)
#' }
#' 
#' # For examples with real data, see help("knecht") or help("alliances").
#' 
#' @importFrom ergm ergm
#' @export
mtergm <- function(formula, constraints = ~ ., returndata = FALSE, 
    verbose = TRUE, ...) {
  
  # call tergmprepare and integrate results as a child environment in the chain
  l <- tergmprepare(formula = formula, offset = FALSE, blockdiag = TRUE, 
      verbose = verbose)
  for (i in 1:length(l$covnames)) {
    assign(l$covnames[i], l[[l$covnames[i]]])
  }
  assign("offsmat", l$offsmat)
  form <- as.formula(l$form, env = environment())
  
  # compile data for creating an mtergm object later; return if necessary
  data <- list()
  for (i in 1:length(l$covnames)) {
    data[[l$covnames[i]]] <- l[[l$covnames[i]]]
  }
  data$offsmat <- l$offsmat
  if (returndata == TRUE) {
    message("Returning a list with data.")
    return(data)
  }
  
  if (verbose == TRUE) {
    message("Estimating...")
    e <- ergm(form, offset.coef = -Inf, constraints = constraints, 
        eval.loglik = TRUE, ...)
  } else {
    e <- suppressMessages(ergm(form, offset.coef = -Inf, 
        constraints = constraints, eval.loglik = TRUE, ...))
  }
  
  # get coefficients and other details
  cf <- coef(e)
  mat <- as.matrix(l$networks)
  if (l$bipartite == TRUE) {
    dyads <- sum(1 - l$offsmat)
  } else {
    dyads <- sum(1 - l$offsmat) - nrow(mat)
  }
  rdf <- dyads - length(cf)
  asyse <- vcov(e, sources = "all")
  se <- sqrt(diag(asyse))
  tval <- coef(e) / se
  pval <- 2 * pt(q = abs(tval), df = rdf, lower.tail = FALSE)
  
  # create mtergm object
  object <- createMtergm(
      coef = cf[-length(cf)],  # do not include NA value for offset matrix
      se = se[-length(se)], 
      pval = pval[-length(pval)], 
      nobs = dyads, 
      time.steps = l$time.steps,
      formula = formula, 
      formula2 = l$form, 
      auto.adjust = l$auto.adjust, 
      offset = TRUE, 
      directed = l$directed, 
      bipartite = l$bipartite, 
      estimate = e$estimate,  # MLE or MPLE
      loglik = e$mle.lik[1], 
      aic = AIC(e), 
      bic = BIC(e), 
      ergm = e, 
      nvertices = l$nvertices, 
      data = data
  )
  
  if (verbose == TRUE) {
    message("Done.")
  }
  
  return(object)
}