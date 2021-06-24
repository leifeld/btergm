# redefine S3 as S4 classes for proper handling as part of the 'tbergm' class
setOldClass(c("bergm", "bergm"))

#' An S4 class to represent a fitted TERGM using Bayesian estimation
#' 
#' An S4 class to represent a fitted TERGM using Bayesian estimation.
#' 
#' \code{tbergm} objects result from Bayesian estimation of a TERGM using the
#' \code{\link{tbergm}} function. They contain the original \code{bergm} object
#' and some additional information.
#'
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
#' @slot estimate Estimate: \code{"bergm"} for Bayesian estimation.
#' @slot bergm The original \code{bergm} object as estimated by the
#'   \code{\link[Bergm]{bergm}} function in the \pkg{Bergm} package.
#' @slot nvertices Number of vertices.
#' @slot data The data after processing by the \code{\link{tergmprepare}}
#'   function.
#'
#' @author Philip Leifeld
#' 
#' @family tergm-classes
#' 
#' @export
setClass(Class = "tbergm",
         slots = c(
           time.steps = "numeric",
           formula = "formula",
           formula2 = "character",
           auto.adjust = "logical",
           offset = "logical",
           directed = "logical",
           bipartite = "logical",
           estimate = "character",
           bergm = "bergm",
           nvertices = "matrix",
           data = "list"
         ),
         validity = function(object) {
           if (!is.numeric(object@time.steps)) {
             stop("'time.steps' must be a numeric value of length 1.")
           }
           if (!"formula" %in% class(object@formula)) {
             stop("'formula' is not a 'formula' object.")
           }
           return(TRUE)
         }
)

#' Constructor for \linkS4class{tbergm} objects
#'
#' Constructor for \linkS4class{tbergm} objects.
#'
#' Create an S4 \linkS4class{tbergm} object using this constructor function.
#'
#' @param estimate Estimate: \code{"bergm"} for Bayesian estimation.
#' @param bergm The original \code{bergm} object as estimated by the
#'   \code{\link[Bergm]{bergm}} function in the \pkg{Bergm} package.
#' @inheritParams createBtergm
#'
#' @author Philip Leifeld
#' 
#' @family tergm-classes
createTbergm <- function(time.steps, formula, formula2, auto.adjust, offset,
                         directed, bipartite, estimate, bergm, nvertices,
                         data) {
  new("tbergm", time.steps = time.steps, formula = formula, formula2 = formula2,
      auto.adjust = auto.adjust, offset = offset, directed = directed,
      bipartite = bipartite, estimate = estimate, bergm = bergm,
      nvertices = nvertices, data = data)
}

#' @describeIn tbergm-class Show the coefficients of a \code{tbergm} object.
#' 
#' @param object A \code{tbergm} object.
#' 
#' @export
setMethod(f = "show", signature = "tbergm", definition = function(object) {
  summary(object@bergm)
}
)

#' @describeIn tbergm-class Return the number of observations saved in a
#'   \code{tbergm} object.
#' 
#' @param object A \code{tbergm} object.
#' 
#' @export
setMethod(f = "nobs", signature = "tbergm", definition = function(object) {
  if (object@bipartite == TRUE) {
    n <- sum(1 - object@data$offsmat)
  } else {
    n <- sum(1 - object@data$offsmat) - nrow(object@data$offsmat)
  }
  t <- object@time.steps
  return(c("Number of time steps" = t, "Number of observations" = n))
}
)

#' @describeIn tbergm-class Return the number of time steps saved in a
#'   \code{tbergm} object.
#' 
#' @param object A \code{tbergm} object.
#' 
#' @export
timesteps.tbergm <- function(object) {
  return(object@time.steps)
}

#' @describeIn tbergm-class Summary of a fitted \code{tbergm} object.
#' 
#' @param object A \code{tbergm} object.
#' @param ... Further arguments for the \code{summary} function in the
#'   \pkg{Bergm} package.
#' 
#' @export
setMethod(f = "summary", signature = "tbergm", definition = function(object,
                                                                     ...) {
  summary(object@bergm, ...)
}
)

#' Estimate a TERGM using Bayesian estimation
#'
#' Estimate a TERGM using Bayesian estimation.
#'
#' The \code{tbergm} function computes TERGMs by Bayesian estimation via
#' blockdiagonal matrices and structural zeros. It acts as a wrapper for the
#' \code{\link[Bergm]{bergm}} function in the \pkg{Bergm} package.
#'
#' @param ... Further arguments to be handed over to the
#'   \code{\link[Bergm]{bergm}} function in the \pkg{Bergm} package.
#' @inheritParams btergm
#'
#' @author Philip Leifeld
#'
#' @seealso \code{\link{btergm}} \code{\link{mtergm}}
#'
#' @references Caimo, Alberto and Nial Friel (2012): Bergm: Bayesian Exponential
#'   Random Graphs in R. \emph{Journal of Statistical Software} 61(2): 1-25.
#'   \doi{10.18637/jss.v061.i02}.
#'
#' @export
tbergm <- function(formula, returndata = FALSE, verbose = TRUE, ...) {
  
  if (!requireNamespace("Bergm", quietly = TRUE)) {
    stop("tbergm requires the 'Bergm' package to be installed.\n",
         "To do this, enter 'install.packages(\"Bergm\")'.")
  }
  
  # call tergmprepare and integrate results as a child environment in the chain
  l <- tergmprepare(formula = formula, offset = FALSE, blockdiag = TRUE,
                    verbose = verbose)
  for (i in 1:length(l$covnames)) {
    assign(l$covnames[i], l[[l$covnames[i]]])
  }
  assign("offsmat", l$offsmat)
  form <- as.formula(l$form, env = environment())

  # compile data for creating a tbergm object later; return if necessary
  dta <- list()
  for (i in 1:length(l$covnames)) {
    dta[[l$covnames[i]]] <- l[[l$covnames[i]]]
  }
  dta$offsmat <- l$offsmat
  if (returndata == TRUE) {
    message("Returning a list with data.")
    return(dta)
  }

  if (verbose == TRUE) {
    message("Estimating...")
  }
  b <- Bergm::bergm(form, offset.coef = -1000, ...)

  # get coefficients and other details
  mat <- as.matrix(l$networks)
  if (l$bipartite == TRUE) {
    dyads <- sum(1 - l$offsmat)
  } else {
    dyads <- sum(1 - l$offsmat) - nrow(mat)
  }
  
  # create tbergm object
  object <- createTbergm(
    time.steps = l$time.steps,
    formula = formula,
    formula2 = l$form,
    auto.adjust = l$auto.adjust,
    offset = TRUE,
    directed = l$directed,
    bipartite = l$bipartite,
    estimate = "bergm",
    bergm = b,
    nvertices = l$nvertices,
    data = dta
  )

  if (verbose == TRUE) {
    message("Done. The Bergm object is stored in the @bergm slot of the fitted object. The data are in @data.")
  }

  return(object)
}