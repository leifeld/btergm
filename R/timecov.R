#' Temporal dependencies for TERGMs
#'
#' Network statistics that span multiple time points.
#' 
#' In addition to the ERGM user terms that can be estimated within a single
#' network (see \link[ergm]{ergm-terms}), the \pkg{btergm} package provides
#' additional model terms that can be used within a formula. These additional
#' statistics span multiple time periods and are therefore called "temporal
#' dependencies." Examples include memory terms (i.e., positive autoregression,
#' dyadic stability, edge innovation, or edge loss), delayed reciprocity or
#' mutuality, and time covariates (i.e., functions of time or interactions with
#' time):
#' \describe{
#'   \item{\code{delrecip(mutuality = FALSE, lag = 1)}}{The \code{delrecip} term
#'     checks for delayed reciprocity. For example, if node \code{j} is tied to
#'     node \code{i} at \code{t = 1}, does this lead to a reciprocation of that
#'     tie back from \code{i} to \code{j} at \code{t = 2}? If
#'     \code{mutuality = TRUE} is set, this extends not only to ties, but also
#'     non-ties. That is, if \code{i} is not tied to \code{j} at \code{t = 1},
#'     will this lead to \code{j} not being tied to \code{i} at \code{t = 2}, in
#'     addition to positively reciprocal patterns over time? The \code{lag}
#'     argument controls the size of the temporal lag: with \code{lag = 1},
#'     reciprocity over one consecutive time period is checked. Note that as
#'     \code{lag} increases, the number of time steps on the dependent variable
#'     decreases.}
#'   \item{\code{memory(type = "stability", lag = 1)}}{Memory terms control for
#'     the impact of a previous network on the current network. Four different
#'     types of memory terms are available: positive autoregression
#'     (\code{type = "autoregression"}) checks whether previous ties are carried
#'     over to the current network; dyadic stability (\code{type = "stability"})
#'     checks whether both edges and non-edges are stable between the previous
#'     and the current network; edge loss (\code{type = "loss"}) checks whether
#'     ties in the previous network have been dissolved and no longer exist in
#'     the current network; and edge innovation (\code{type = "innovation"})
#'     checks whether previously unconnected nodes have the tendency to become
#'     tied in the current network. The \code{lag} argument accepts integer
#'     values and controls whether the comparison is made with the previous
#'     network (\code{lag = 1}), the pre-previous network (\code{lag = 2}) etc.
#'     Note that as \code{lag} increases, the number of time steps on the
#'     dependent variable decreases.}
#'   \item{\code{timecov(x = NULL, minimum = 1, maximum = NULL,
#'     transform = function(t) t)}}{The \code{timecov} model term checks for
#'     linear or non-linear time trends with regard to edge formation.
#'     Optionally, this can be combined with a covariate to create an
#'     interaction effect between a dyadic covariate and time in order to test
#'     whether the importance of a covariate increases or decreases over time.
#'     In the default case, edges modeled as being linearly increasingly
#'     important over time. By tweaking the \code{transform} function,
#'     arbitrary functional forms of time can be tested. For example,
#'     \code{transform = sqrt} (for a geometrically decreasing time effect),
#'     \code{transform = function(x) x^2} (for a geometrically increasing time
#'     effect), \code{transform = function(t) t} (for a linear time trend) or
#'     polynomial functional forms (e.g., \code{0 + (1 * t) + (1 * t^2)}) can
#'     be used. For time steps below the \code{minimum} value and above the
#'     \code{maximum} value, the time covariate is set to 0. These arguments
#'     can be used to create step-wise, discrete effects, for example to use a
#'     value of 0 up to an external event and 1 from that event onwards in
#'     order to control for influences of external events.}
#' }
#'
#' @references Leifeld, Philip, Skyler J. Cranmer and Bruce A. Desmarais (2017):
#'   Temporal Exponential Random Graph Models with btergm: Estimation and
#'   Bootstrap Confidence Intervals. \emph{Journal of Statistical Software}
#'   83(6): 1-36. \doi{10.18637/jss.v083.i06}.
#'
#' @name tergm-terms
#' @aliases tergm-terms memory timecov delrecip
NULL

#' Transform a covariate using a function of time
#' 
#' Transform a covariate using a function of time.
#' 
#' The \code{timecov} model term checks for linear or non-linear time trends
#' with regard to edge formation. Optionally, this can be combined with a
#' covariate to create an interaction effect between a dyadic covariate and time
#' in order to test whether the importance of a covariate increases or decreases
#' over time. The function can either be used in a formula with
#' \code{\link{btergm}}, \code{\link{mtergm}}, or \code{\link{tbergm}}, or it
#' can be executed directly for manual inclusion of the results as a covariate.
#'
#' @param covariate The list of networks or matrices for which to create a time
#'   covariate. This can be the list of networks on the left-hand side of the
#'   formula, in which case a time trend is created as a covariate list of
#'   matrices, or it can be a list of networks or matrices that is used as a
#'   dyadic covariate on the right-hand side of the formula, in which case an
#'   interaction effect between the time trend and the covariate is created. If
#'   used as a model term inside a formula, \code{covariate = NULL} is
#'   permitted, in which case the networks on the left-hand side will be used to
#'   form a time trend.
#' @param minimum,maximum For time steps below the \code{minimum} value and
#'   above the \code{maximum} value, the time covariate is set to 0. These
#'   arguments can be used to create step-wise, discrete effects, for example to
#'   use a value of 0 up to an external event and 1 from that event onwards in
#'   order to control for influences of external events.
#' @param transform In the default case, edges are modeled as being linearly
#'   increasingly important over time (i.e., a linear time trend). By tweaking
#'   the \code{transform} function, arbitrary functional forms of time can be
#'   tested. For example, \code{transform = sqrt} (for a geometrically
#'   decreasing time effect), \code{transform = function(x) x^2} (for a
#'   geometrically increasing time effect), \code{transform = function(t) t}
#'   (for a linear time trend) or polynomial functional forms (e.g.,
#'   \code{transform = function(t) 0 + (1 * t) + (1 * t^2)}) can be used.
#' @param onlytime If \code{TRUE}, return a time trend only. If \code{FALSE},
#'   return an interaction between the time trend and the covariate. Note that
#'   the model term may need to be called twice or more inside a formula: one
#'   time to create the time trend main effect and one time for each
#'   interaction term; you also need to include the main effects for the
#'   covariates separately using \code{edgecov} or similar terms.
#'
#' @describeIn tergm-terms Time trends and temporal covariate interactions
#' @aliases timecov
#' 
#' @importFrom network is.network
#' @export
timecov <- function(covariate, minimum = 1, maximum = length(covariate), 
    transform = function(t) 1 + (0 * t) + (0 * t^2), onlytime = FALSE) {
  if (!"list" %in% class(covariate)) {
    stop("'covariate' must be a list of matrices or network objects.")
  }
  for (i in 1:length(covariate)) {
    if (is.network(covariate[[i]])) {
      covariate[[i]] <- as.matrix(covariate[[i]])
    } else if (!is.matrix(covariate[[i]])) {
      stop("'covariate' must be a list of matrices or network objects.")
    }
  }
  if (is.null(minimum) || is.null(maximum) || !is.numeric(minimum) || 
      !is.numeric(maximum) || length(minimum) > 1 || length(maximum) > 1) {
    stop("'minimum' and 'maximum' must be single numeric values.")
  }
  if (is.null(transform)) {
    transform <- function(t) 1 + (0 * t) + (0 * t^2)
  } else if (!is.function(transform)) {
    stop("'transform' must be a function.")
  }
  l <- 1:length(covariate)
  values <- transform(l)  # apply transformation of time
  if (is.null(values) || any(is.null(values)) || any(!is.finite(values)) || 
      any(is.na(values)) || any(!is.numeric(values))) {
    stop("The 'transform' function produces non-numeric values.")
  }
  values <- values * (l >= minimum) * (l <= maximum)  # interval dummy
  timecov <- list()
  for (i in 1:length(l)) {  # create matrix at each time step
    if (onlytime == FALSE) {
      timecov[[i]] <- covariate[[i]] * matrix(values[i], nrow = 
          nrow(covariate[[i]]), ncol = ncol(covariate[[i]]))
    } else {
      timecov[[i]] <- matrix(values[i], nrow = nrow(covariate[[i]]), 
          ncol = ncol(covariate[[i]]))
    }
  }
  for (i in 1:length(timecov)) {
    rownames(timecov[[i]]) <- rownames(covariate[[i]])
    colnames(timecov[[i]]) <- colnames(covariate[[i]])
  }
  return(timecov)
}
