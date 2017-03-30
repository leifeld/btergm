
# display version number and date when the package is loaded
.onAttach <- function(libname, pkgname) {
  desc  <- packageDescription(pkgname, libname)
  packageStartupMessage(
    'Package:  btergm\n', 
    'Version:  ', desc$Version, '\n', 
    'Date:     ', desc$Date, '\n', 
    'Authors:  Philip Leifeld (University of Glasgow)\n',
    '          Skyler J. Cranmer (The Ohio State University)\n',
    '          Bruce A. Desmarais (Pennsylvania State University)\n'
  )
}


# redefine S3 as S4 classes for proper handling as part of the 'btergm' class
setOldClass(c("boot", "boot"))


# an S4 class for btergm objects
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
          stop("The sample size does not correspond to the 'R' parameter.")
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


# constructor for btergm objects
createBtergm <- function(coef, boot, R, nobs, time.steps, formula, 
    formula2, response, effects, weights, auto.adjust, offset, 
    directed, bipartite, nvertices, data) {
  new("btergm", coef = coef, boot = boot, R = R, nobs = nobs, 
      time.steps = time.steps, formula = formula, formula2 = formula2, 
      response = response, effects = effects, weights = weights, 
      auto.adjust = auto.adjust, offset = offset, directed = directed, 
      bipartite = bipartite, nvertices = nvertices, data = data)
}


# define show method for pretty output of btergm objects
setMethod(f = "show", signature = "btergm", definition = function(object) {
    message("MLE Coefficients:")
    print(object@coef)
  }
)


# define coef method for extracting coefficients from btergm objects
setMethod(f = "coef", signature = "btergm", definition = function(object, 
      invlogit = FALSE, ...) {
    if (invlogit == FALSE) {
      return(object@coef)
    } else {
      return(1 / (1 + exp(-object@coef)))
    }
  }
)


# define nobs method for extracting number of observations from btergm objects
setMethod(f = "nobs", signature = "btergm", definition = function(object) {
    n <- object@nobs
    t <- object@time.steps
    rep <- object@R
    return(c("Number of time steps" = t, "Number of observations" = n, 
        "Bootstrap replications" = rep))
  }
)


# function which can extract a coefficient matrix with SEs and p values
btergm.se <- function(object, print = FALSE) {
  co <- object@coef
  #sdev <- apply(object@boot$t, 2, sd) # old; now use deviation from estimate:
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


# confint method for btergm objects
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
    ci <- cbind(cf, t(ci))
    if (class(ci) == "numeric") {
      ci.nam <- names(ci)
      ci <- matrix(ci, nrow = 1)
      colnames(ci) <- ci.nam
      rownames(ci) <- names(cf)
    }
    ci <- ci[parm, ]
    if (class(ci) != "matrix") {
      ci <- matrix(ci, ncol = 3)
      rownames(ci) <- parm
    }
    label1 <- paste0(100 * (1 - level) / 2, "%")
    label2 <- paste0(100 * (1 - (1 - level) / 2), "%")
    colnames(ci) <- c("Estimate", label1, label2)
    return(ci)
  }
)


# function which can extract the number of time steps
timesteps.btergm <- function(object) {
  return(object@time.steps)
}


# define summary method for pretty output of btergm objects
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
    printCoefmat(cmat, cs.ind = 1, tst.ind = 2:3)
  }
)


# TERGM by bootstrapped pseudolikelihood
btergm <- function(formula, R = 500, offset = FALSE, 
    returndata = FALSE, parallel = c("no", "multicore", 
    "snow"), ncpus = 1, cl = NULL, verbose = TRUE, ...) {
  
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
      model <- ergm::ergm.getmodel(form, nw, initialfit = TRUE)
      Clist <- ergm::ergm.Cprepare(nw, model)
      Clist.miss <- ergm::ergm.design(nw, model, verbose = FALSE)
      pl <- ergm::ergm.pl(Clist, Clist.miss, model, theta.offset = 
          c(rep(FALSE, length(l$rhs.terms) - 1), TRUE), verbose = FALSE, 
          control = ergm::control.ergm(init = c(rep(NA, 
          length(l$rhs.terms) - 1), 1)))
      Y <- c(Y, pl$zy[pl$foffset == 0])
      X <- rbind(X, cbind(data.frame(pl$xmat[pl$foffset == 0, ], 
          check.names = FALSE), i))
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
      mpli <- ergm::ergmMPLE(form)
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
  
  if (returndata == TRUE) {
    return(cbind(Y, x))
  }
  
  # create sparse matrix and compute start values for GLM
  xsparse <- Matrix(as.matrix(x), sparse = TRUE)
  if (ncol(xsparse) == 1) {
    stop("At least two model terms must be provided to estimate a TERGM.")
  }
  est <- speedglm.wfit(y = Y, X = xsparse, weights = W, offset = O, 
      family = binomial(link = logit), sparse = TRUE)
  startval <- coef(est)
  nobs <- est$n
  # define function for bootstrapping and estimation
  estimate <- function(unique.time.steps, bsi, Yi = Y, xsparsei = xsparse, 
      Wi = W, Oi = O, timei = X$time, startvali = startval) {
    indic <- unlist(lapply(bsi, function(x) which(timei == x)))
    tryCatch(
      expr = {
        return(coef(speedglm.wfit(y = Yi[indic], X = xsparsei[indic, ], 
            weights = Wi[indic], offset = Oi[indic], 
            family = binomial(link = logit), sparse = TRUE, start = startvali)))
      }, 
      error = function(e) {
        # when fitted probabilities of 0 or 1 occur or when the algorithm does 
        # not converge, use glm because it only throws a warning, not an error
        return(coef(glm.fit(y = Yi[indic], x = as.matrix(x)[indic, ], 
            weights = Wi[indic], offset = Oi[indic], 
            family = binomial(link = logit))))
      }, 
      warning = function(w) {
        warning(w)
      }, 
      finally = {}
    )
  }
  
  # run the estimation (single-core or parallel)
  coefs <- boot(unique.time.steps, estimate, R = R, Yi = Y, xsparsei = xsparse, 
      Wi = W, Oi = O, timei = X$time, startvali = startval, 
      parallel = parallel, ncpus = ncpus, cl = cl, ...)
  rm(X)
  #if (nrow(coefs$t) == 1) { # in case there is only one model term
  #  coefs <- t(coefs)
  #}
  if (ncol(coefs$t) == 1 && length(term.names) > 1 
      && coefs$t[1, 1] == "glm.fit: algorithm did not converge") {
    stop(paste("Algorithm did not converge. There might be a collinearity ", 
        "between predictors and/or dependent networks at one or more time", 
        "steps."))
  }
  
  # create and return btergm object
  colnames(coefs$t) <- term.names[1:(length(term.names) - 1)]
  names(startval) <- colnames(coefs$t)
  data <- list()
  for (i in 1:length(l$covnames)) {
    data[[l$covnames[i]]] <- l[[l$covnames[i]]]
  }
  data$offsmat <- l$offsmat
  
  btergm.object <- createBtergm(startval, coefs, R, nobs, l$time.steps, 
      formula, l$form, Y, x, W, l$auto.adjust, offset, l$directed, l$bipartite, 
      nvertices = l$nvertices, data)
  if (verbose == TRUE) {
    message("Done.")
  }
  return(btergm.object)
}


# simulation of new networks based on a btergm fit
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
  } else if (!is.numeric(index)) {
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
  s <- simulate.formula(form, nsim = nsim, seed = seed, coef = coef, 
      verbose = verbose, ...)
  if ("btergm" %in% class(object)) {
    return(s)
  } else if ("mtergm" %in% class(object)) {
    r1 <- sum(object@nvertices[1, 1:(index - 1)]) + 1
    c1 <- sum(object@nvertices[2, 1:(index - 1)]) + 1
    r2 <- sum(object@nvertices[1, 1:index])
    c2 <- sum(object@nvertices[2, 1:index])
    s <- lapply(s, function(sim) network(as.matrix(sim)[r1:r2, c1:c2], 
        bipartite = object@bipartite, directed = object@directed))
    return(s)
  }
}
simulate.mtergm <- simulate.btergm  # create a copy for mtergm objects
