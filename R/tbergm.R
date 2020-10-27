# redefine S3 as S4 classes for proper handling as part of the 'tbergm' class
setOldClass(c("bergm", "bergm"))

# an S4 class for tbergm objects
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


# constructor for tbergm objects
createTbergm <- function(time.steps, formula, formula2, auto.adjust, offset,
                         directed, bipartite, estimate, bergm, nvertices,
                         data) {
  new("tbergm", time.steps = time.steps, formula = formula, formula2 = formula2,
      auto.adjust = auto.adjust, offset = offset, directed = directed,
      bipartite = bipartite, estimate = estimate, bergm = bergm,
      nvertices = nvertices, data = data)
}


# define show method for pretty output of tbergm objects
setMethod(f = "show", signature = "tbergm", definition = function(object) {
  summary(object)
}
)


# define nobs method for extracting number of observations from tbergm objects
setMethod(f = "nobs", signature = "tbergm", definition = function(object) {
  n <- object@nobs
  t <- object@time.steps
  return(c("Number of time steps" = t, "Number of observations" = n))
}
)


# function which can extract the number of time steps
timesteps.tbergm <- function(object) {
  return(object@time.steps)
}


# define summary method for pretty output of mtergm objects
setMethod(f = "summary", signature = "tbergm", definition = function(object,
                                                                     ...) {
  summary(object, ...)
}
)


# Bayesian estimation (wrapper for the bergm function in the Berg package)
tbergm <- function(formula, returndata = FALSE, verbose = TRUE, ...) {

  # call tergmprepare and integrate results as a child environment in the chain
  l <- tergmprepare(formula = formula, offset = FALSE, blockdiag = TRUE,
                    verbose = verbose)
  for (i in 1:length(l$covnames)) {
    assign(l$covnames[i], l[[l$covnames[i]]])
  }
  assign("offsmat", l$offsmat)
  form <- as.formula(l$form, env = environment())

  # compile data for creating a tbergm object later; return if necessary
  data <- list()
  for (i in 1:length(l$covnames)) {
    data[[l$covnames[i]]] <- l[[l$covnames[i]]]
  }
  data$offsmat <- l$offsmat
  if (returndata == TRUE) {
    message("Returning a list with data.")
    return(data)
  }

  if (!requireNamespace("Bergm", quietly = TRUE)) {
    stop("tbergm requires the 'Bergm' package to be installed.\n",
         "To do this, enter 'install.packages(\"Bergm\")'.")
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
    data = data
  )

  if (verbose == TRUE) {
    message("Done.")
  }

  return(object)
}
