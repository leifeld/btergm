#' Check for degeneracy in fitted TERGMs
#'
#' Check for degeneracy in fitted TERGMs.
#' 
#' The methods for the generic \code{degeneracy} function implement a degeneracy
#' check for \code{btergm} and \code{mtergm} objects. For \code{btergm}, this
#' works by comparing the global statistics of simulated networks to those of
#' the observed networks at each observed time step. If the global statistics
#' differ significantly, this is indicated by small p-values. If there are many
#' significant results, this indicates degeneracy. For \code{mtergm}, the
#' \code{mcmc.diagnostics} function from the \pkg{ergm} package is used.
#' 
#' @param object A \code{btergm} or \code{mtergm} object, as estimated using the
#'   \code{btergm} or \code{mtergm} function.
#' @param nsim The number of networks to be simulated at each time step. This
#'   number should be sufficiently large for a meaningful comparison. If
#'   possible, much more than 1,000 simulations.
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
#' @param verbose Print details?
#' @param x A \code{degeneracy} object created by the \code{checkdegeneracy}
#'   function.
#' @param center If \code{TRUE}, print/plot the simulated minus the target
#'   statistics, with an expected value of 0 in a non-degenerate model. If
#'   \code{FALSE}, print/plot the distribution of simulated statistics and show
#'   the target statistic separately.
#' @param t Time indices to include, e.g., \code{t = 2:4} for time steps 2 to 4.
#' @param terms Indices of the model terms to include, e.g., \code{terms = 1:3}
#'   includes the first three statistics.
#' @param vbar Show vertical bar for target statistic in histogram.
#' @param main Main title of the plot.
#' @param xlab Label on the x-axis. Defaults to the name of the statistic.
#' @param target.col Color of the vertical bar for the target statistic.
#'   Defaults to red.
#' @param target.lwd Line width of the vertical bar for the target statistic.
#'   Defaults to 3.
#' @param ... Arbitrary further arguments for subroutines.
#'
#' @return A list with target statistics and simulations.
#' 
#' @references
#' Hanneke, Steve, Wenjie Fu and Eric P. Xing (2010): Discrete Temporal Models
#' of Social Networks. \emph{Electronic Journal of Statistics} 4: 585--605.
#' \doi{10.1214/09-EJS548}.
#' 
#' Leifeld, Philip, Skyler J. Cranmer and Bruce A. Desmarais (2018): Temporal
#' Exponential Random Graph Models with btergm: Estimation and Bootstrap
#' Confidence Intervals. \emph{Journal of Statistical Software} 83(6): 1-36.
#' \doi{10.18637/jss.v083.i06}.
#'
#' @docType methods
#' @aliases checkdegeneracy-methods
#' @importFrom ergm simulate_formula control.simulate.formula mcmc.diagnostics
#' @export
setGeneric("checkdegeneracy", function(object, ...)
  standardGeneric("checkdegeneracy"), package = "btergm")

#' @noRd
checkdegeneracy.mtergm <- function(object, ...) {
  ergm::mcmc.diagnostics(object@ergm, ...)
}

#' @rdname checkdegeneracy
setMethod("checkdegeneracy", signature = className("mtergm", "btergm"),
          definition = checkdegeneracy.mtergm)

#' @noRd
checkdegeneracy.btergm <- function(object, nsim = 1000, MCMC.interval = 1000,
    MCMC.burnin = 10000, verbose = FALSE) {
  if (nsim < 2) {
    stop("The 'nsim' argument must be greater than 1.")
  }

  # call tergmprepare and integrate results in local environment
  l <- tergmprepare(formula = getformula(object), offset = object@offset,
      verbose = verbose)
  for (i in 1:length(l$covnames)) {
    assign(l$covnames[i], l[[l$covnames[i]]])
  }
  assign("offsmat", l$offsmat)
  form <- as.formula(l$form)
  offset <- object@offset
  target <- l$networks

  # extract coefficients from object
  if ("btergm" %in% class(object) && offset == TRUE) {
    coefs <- c(coef(object), -Inf)  # -Inf for offset matrix
  } else {
    coefs <- coef(object)
  }

  # adjust formula at each step, and simulate networks
  sim <- list()
  target.stats <- list()
  degen <- list()
  for (index in 1:l$time.steps) {
    i <- index
    if (verbose == TRUE) {
      f.i <- gsub("\\[\\[i\\]\\]", paste0("[[", index, "]]"),
          paste(deparse(form), collapse = ""))
      f.i <- gsub("\\s+", " ", f.i)
      f.i <- gsub("^networks", l$lhs.original, f.i)
      message(paste("Simulating", nsim,
          "networks from the following formula:\n", f.i, "\n"))
    }
    target.stats[[index]] <- summary(statnet.common::filter_rhs.formula(form),
        response = NULL)
    degen[[index]] <- ergm::simulate_formula(form,
                                             nsim = nsim,
                                             coef = coefs,
                                             statsonly = TRUE,
                                             control = control.simulate.formula(MCMC.interval = MCMC.interval,
                                                                                MCMC.burnin = MCMC.burnin))
    if (offset == TRUE || "mtergm" %in% class(object)) {
      degen[[i]] <- degen[[i]][, -ncol(degen[[i]])]  # remove offset statistic
    }
  }

  if (verbose == TRUE) {
    message("Checking degeneracy...")
  }
  object <- list()
  object$target.stats <- target.stats
  object$sim <- degen
  class(object) <- "degeneracy"
  if (verbose == TRUE) {
    message("Done.")
  }
  return(object)
}

#' @rdname checkdegeneracy
setMethod("checkdegeneracy", signature = className("btergm", "btergm"),
    definition = checkdegeneracy.btergm)


#' @rdname checkdegeneracy
print.degeneracy <- function(x, center = FALSE, t = 1:length(x$sim),
    terms = 1:length(x$target.stats[[1]]), ...) {
  for (i in t) {
    message(paste0("\nDegeneracy check for network ", i, ":"))
    if (center == TRUE) {
      sm <- coda::as.mcmc.list(coda::as.mcmc(x$sim[[i]]))
      sm <- statnet.common::sweep.mcmc.list(sm, x$target.stats[[i]], "-")[[1]] # diff
      q <- t(apply(as.matrix(sm), 2, function(x) {
        quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
      }))
    } else {
      q <- t(apply(x$sim[[i]], 2, function(x) {
        quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
      }))
      q <- cbind("obs" = x$target.stats[[i]], q)
    }
    rn <- rownames(q)[terms]
    cn <- colnames(q)
    q <- q[terms, ]
    if (!is.matrix(q)) {
      q <- matrix(q, nrow = 1)
      rownames(q) <- rn
      colnames(q) <- cn
    }
    printCoefmat(q)
  }
}

#' @rdname checkdegeneracy
plot.degeneracy <- function(x, center = TRUE, t = 1:length(x$sim),
    terms = 1:length(x$target.stats[[1]]), vbar = TRUE, main = NULL,
    xlab = NULL, target.col = "red", target.lwd = 3, ...) {
  for (i in t) {
    for (j in terms) {
      if (is.null(main)) {
         m <- paste0(colnames(x$sim[[i]])[j], " at t = ", i)
      } else {
         m <- main
      }
      if (is.null(xlab)) {
         xl <- colnames(x$sim[[i]])[j]
      } else {
         xl <- xlab
      }
      if (center == FALSE) {
        hist(x$sim[[i]][, j], main = m, xlab = xl, ...)
        if (vbar == TRUE) {
          abline(v = x$target.stats[[i]][j], col = target.col, lwd = target.lwd)
        }
      } else {
        centered <- x$sim[[i]][, j] - x$target.stats[[i]][j]
        hist(centered, main = m, xlab = xl, ...)
        if (vbar == TRUE) {
          abline(v = 0, col = target.col, lwd = target.lwd)
        }
      }
    }
  }
}