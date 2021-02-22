#' Extract the formula from a model
#' 
#' Extract the model formula from a fitted object.
#' 
#' The \code{getformula} function will extract the formula from a fitted model.
#' 
#' @param x A fitted model.
#'
#' @docType methods
#' @aliases getformula-methods
#' @export
setGeneric("getformula", function(x) standardGeneric("getformula"),
           package = "btergm")

#' @describeIn getformula Extract the formula from an \code{ergm} object.
#' @export
setMethod("getformula", signature = className("ergm", "ergm"), 
          definition = function(x) x$formula)

#' @describeIn getformula Extract the formula from a \code{btergm} object.
#' @export
setMethod("getformula", signature = className("btergm", "btergm"), 
          definition = function(x) x@formula)

#' @describeIn getformula Extract the formula from an \code{mtergm} object.
#' @export
setMethod("getformula", signature = className("mtergm", "btergm"), 
          definition = function(x) x@formula)

#' @describeIn getformula Extract the formula from a \code{tbergm} object.
#' @export
setMethod("getformula", signature = className("tbergm", "btergm"), 
          definition = function(x) x@formula)