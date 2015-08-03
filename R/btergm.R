
# display version number and date when the package is loaded
.onAttach <- function(libname, pkgname) {
  desc  <- packageDescription(pkgname, libname)
  packageStartupMessage(
    'Package:  btergm\n', 
    'Version:  ', desc$Version, '\n', 
    'Date:     ', desc$Date, '\n', 
    'Authors:  Philip Leifeld (Eawag (ETH) and University of Bern)\n',
    '          Skyler J. Cranmer (The Ohio State University)\n',
    '          Bruce A. Desmarais (Penn State University)\n'
  )
}


# an S4 class for btergm objects
setClass(Class = "btergm", 
    representation = representation(
        coef = "numeric",
        bootsamp = "matrix",
        R = "numeric",
        nobs = "numeric", 
        time.steps = "numeric",
        formula = "formula",
        response = "integer",
        effects = "data.frame", 
        weights = "numeric", 
        auto.adjust = "logical", 
        offset = "logical", 
        directed = "logical", 
        bipartite = "logical"
    ), 
    validity = function(object) {
        if (!"numeric" %in% class(object@coef)) {
          stop("'coef' must be a 'numeric' object.")
        }
        if (!"matrix" %in% class(object@bootsamp)) {
          stop("'bootsamp' must be a 'matrix' object.")
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
        if (nrow(object@bootsamp) != object@R) {
          stop("The sample size does not correspond to the 'R' parameter.")
        }
        if (length(object@coef) != ncol(object@bootsamp)) {
          stop("Number of terms differs between 'bootsamp' and 'coef'")
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
createBtergm <- function(coef, bootsamp, R, nobs, time.steps, 
    formula, response, effects, weights, auto.adjust, offset, directed, 
    bipartite) {
  new("btergm", coef = coef, bootsamp = bootsamp,
      R = R, nobs = nobs, time.steps = time.steps, formula = formula, 
      response = response, effects = effects, weights = weights, 
      auto.adjust = auto.adjust, offset = offset, directed = directed, 
      bipartite = bipartite)
}


# define show method for pretty output of btergm objects
setMethod(f = "show", signature = "btergm", definition = function(object) {
    message("MLE Coefficients:")
    print(object@coef)
  }
)


# define coef method for extracting coefficients from btergm objects
setMethod(f = "coef", signature = "btergm", definition = function(object, ...) {
    return(object@coef)
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
  #sdev <- apply(object@bootsamp, 2, sd) # old; now use deviation from estimate:
  sdev <- numeric()
  for (i in 1:ncol(object@bootsamp)) {
    currentcol <- numeric()
    for (j in 1:nrow(object@bootsamp)) {
      currentcol[j] <- (object@bootsamp[j, i] - co[i])^2
    }
    sdev[i] <- sqrt(sum(currentcol) / length(currentcol))
  }
  zval <- (0 - apply(object@bootsamp, 2, mean)) / sdev
  pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
  cmat <- cbind(co, sdev, zval, pval)
  colnames(cmat) <- c("Estimate", "Std.Err", "Z value", "Pr(>z)")
  warning(paste("Standard errors and p values may be misleading because the",
      "distribution of bootstrapped thetas may not be normal. Please rely on",
      "the confidence intervals instead or make sure the thetas are normally",
      "distributed (e.g., using qqnorm(object@bootsamp[, 1]) etc."))
  if (print == TRUE) {
    printCoefmat(cmat)
  } else {
    return(cmat)
  }
}


# confint method for btergm objects
setMethod(f = "confint", signature = "btergm", definition = function(object, 
    parm, level = 0.95, ...) {
    cf <- coef(object)
    pnames <- names(cf)
    if (missing(parm)) {
      parm <- pnames
    } else if (is.numeric(parm)) {
      parm <- pnames[parm]
    }
    samples <- object@bootsamp[complete.cases(object@bootsamp), ]
    if (class(samples) == "numeric") {  # only one model term
      samples <- as.matrix(samples, ncol = 1)
    }
    n.orig <- nrow(object@bootsamp)
    n.ret <- nrow(samples)
    perc <- 100 * (n.orig - n.ret) / n.orig
    if (nrow(samples) != nrow(object@bootsamp)) {
      warning(paste0("Too little variation in the model. ", n.orig - n.ret, 
          " replications (", perc, "%) are dropped from CI estimation."))
    }
    ci <- t(apply(samples, 2, function(object) quantile(object, 
        c(((1 - level) / 2), 1 - ((1 - level) / 2)))))
    ci <- cbind(cf, ci)[parm, ]
    if (class(ci) == "numeric") {
      ci.nam <- names(ci)
      ci <- matrix(ci, nrow = 1)
      colnames(ci) <- ci.nam
      rownames(ci) <- names(cf)
    }
    colnames(ci)[1] <- "Estimate"
    return(ci)
  }
)


# function which can extract the number of time steps
btergm.timesteps <- function(object) {
  return(object@time.steps)
}


# define summary method for pretty output of btergm objects
setMethod(f = "summary", signature = "btergm", definition = function(object, 
    level = 0.95, ...) {
    message(paste(rep("=", 26), collapse=""))
    message("Summary of model fit")
    message(paste(rep("=", 26), collapse=""))
    message(paste("\nFormula:  ", gsub("\\s+", " ", 
        paste(deparse(object@formula), collapse = "")), "\n"))
    message(paste("Time steps:", object@time.steps, "\n"))
    message(paste("Bootstrapping sample size:", object@R, "\n"))
    
    message(paste0("Estimates and ", 100 * level, "% confidence intervals:"))
    cmat <- confint(object, level = level, ...)
    printCoefmat(cmat, cs.ind = 1, tst.ind = 2:3)
  }
)


# helper function for adjusting matrix dimensions and creating an offset matrix
tergmprepare <- function(formula, offset = TRUE, blockdiag = FALSE, 
    verbose = TRUE) {
  
  # extract response networks and make sure they are saved in a list
  env <- new.env()
  env$lhs.original <- deparse(formula[[2]])  # for reporting purposes later on
  env$networks <- eval(parse(text = deparse(formula[[2]])))
  if (class(env$networks) == "list" || class(env$networks) == "network.list") {
    # do nothing
  } else {
    env$networks <- list(env$networks)
  }
  
  # convert list elements to matrices if unknown data type
  for (i in 1:length(env$networks)) {
    if (!class(env$networks[[i]]) %in% c("network", "matrix", "list")) {
      tryCatch(
        {
          env$networks[[i]] <- as.matrix(env$networks[[i]])
        }, 
        error = function(cond) {
          stop(paste("Object", i, "could not be converted to a matrix."))
        }
      )
    }
  }
  
  # extract additional information
  env$num.vertices <- max(sapply(env$networks, function(x) 
      network::get.network.attribute(network::network(x), "n")))  # number nodes
  if (is.network(env$networks[[1]])) {
    env$directed <- network::is.directed(env$networks[[1]])  # directed?
    env$bipartite <- network::is.bipartite(env$networks[[1]])  # bipartite?
  } else {
    if (xergm.common::is.mat.directed(as.matrix(env$networks[[1]]))) {
      env$directed <- TRUE
    } else {
      env$directed <- FALSE
      #stop(paste("The dependent networks seem to be undirected. In this", 
      #    "case, please store them as a list of network objects."))
    }
    if (xergm.common::is.mat.onemode(as.matrix(env$networks[[1]]))) {
      env$bipartite <- FALSE
    } else {
      env$bipartite <- TRUE
    }
  }
  
  # adjust and disassemble formula
  env$form <- update.formula(formula, networks[[i]] ~ .)
  env$time.steps <- length(env$networks)
  tilde <- deparse(env$form[[1]])
  lhs <- deparse(env$form[[2]])
  rhs <- paste(deparse(env$form[[3]]), collapse = "")  # for long formulae
  rhs <- gsub("\\s+", " ", rhs)
  
  # parse rhs of formula and add indices to edgecov and dyadcov terms
  env$rhs.terms <- strsplit(rhs, "\\s*(\\+|\\*)\\s*")[[1]]
  rhs.indices <- gregexpr("\\+|\\*", rhs)[[1]]
  if (length(rhs.indices) == 1 && rhs.indices < 0) {
    rhs.operators <- character()
  } else {
    rhs.operators <- substring(rhs, rhs.indices, rhs.indices)
  }
  
  # preprocess dyadcov and edgecov terms
  covnames <- character()
  for (k in 1:length(env$rhs.terms)) {
    if (grepl("((edge)|(dyad))cov", env$rhs.terms[k])) {
      
      # if edgecov or dyadcov, split up into components
      if (grepl(",\\s*?((attr)|\\\")", env$rhs.terms[k])) { # with attrib arg.
        s <- "((?:offset\\()?((edge)|(dyad))cov\\()([^\\)]+)((,\\s*a*.*?)\\)(?:\\))?)"
      } else { # without attribute argument
        s <- "((?:offset\\()?((edge)|(dyad))cov\\()([^\\)]+)((,*\\s*a*.*?)\\)(?:\\))?)"
      }
      x1 <- sub(s, "\\1", env$rhs.terms[k], perl = TRUE)  # before the covariate
      x2 <- sub(s, "\\5", env$rhs.terms[k], perl = TRUE)  # name of the cov.
      x3 <- sub(s, "\\6", env$rhs.terms[k], perl = TRUE)  # after the covariate
      x.current <- eval(parse(text = x2))
      type <- class(x.current)
      env$covnames <- c(env$covnames, x2)
      env[[x2]] <- x.current
      
      # add brackets if necessary, convert to list, and reassemble rhs term
      if (grepl("[^\\]]\\]$", x2)) {
        # time-varying covariate with given indices (e.g., formula[1:5])
        env$rhs.terms[k] <- paste(x1, x2, x3, sep = "")
        if (type %in% c("matrix", "network", "dgCMatrix", "dgTMatrix", 
            "dsCMatrix", "dsTMatrix", "dgeMatrix")) {
          x.current <-list(x.current)
          env[[x2]] <- x.current
        }
        if (length(x.current) != env$time.steps) {
          stop(paste(x2, "has", length(x.current), "elements, but there are", 
              env$time.steps, "networks to be modeled."))
        }
        x2 <- paste(x2, "[[i]]", sep = "")
      } else if (type %in% c("matrix", "network", "dgCMatrix", "dgTMatrix", 
          "dsCMatrix", "dsTMatrix", "dgeMatrix")) {
        # time-independent covariate
        if (!type %in% c("matrix", "network")) {
          x.current <- as.matrix(x.current)
        }
        env[[x2]] <- list()
        for (i in 1:env$time.steps) {
          env[[x2]][[i]] <- x.current
        }
        x2 <- paste(x2, "[[i]]", sep = "")
        env$rhs.terms[k] <- paste(x1, x2, x3, sep = "")
      } else if (type == "list" || type == "network.list") {
        # time-varying covariate
        if (length(x.current) != env$time.steps) {
          stop(paste(x2, "has", length(get(x2)), "elements, but there are", 
              env$time.steps, "networks to be modeled."))
        }
        x2 <- paste(x2, "[[i]]", sep = "")
        env$rhs.terms[k] <- paste(x1, x2, x3, sep = "")
      } else {  # something else --> try to convert to matrix list
        tryCatch(
          {
            env[[x2]] <- list(rep(as.matrix(x.current)), env$time.steps)
          }, 
          error = function(cond) {
            stop(paste0("Object '", x2, 
                "' could not be converted to a matrix."))
          }
        )
      }
    }
  }
  
  # determine and report initial dimensions of networks and covariates
  env$covnames <- c("networks", env$covnames)
  if (verbose == TRUE) {
    if (length(env$covnames) > 1) {
      dimensions <- lapply(lapply(env$covnames, function(x) env[[x]]), 
          function(y) sapply(y, function(z) dim(as.matrix(z))))
      rownames(dimensions[[1]]) <- paste(env$lhs.original, c("(row)", "(col)"))
      for (i in 2:length(dimensions)) {
        rownames(dimensions[[i]]) <- c(paste(env$covnames[i], "(row)"), 
            paste(env$covnames[i], "(col)"))
      }
      dimensions <- do.call(rbind, dimensions)
      colnames(dimensions) <- paste0("t=", 1:length(env$networks))
      message("\nInitial dimensions of the network and covariates:")
      print(dimensions)
    } else {
      message("\nNo covariates provided.")
    }
  }
  
  # determine whether covariate dimensions need to be automatically adjusted
  env$auto.adjust <- FALSE
  if (length(env$covnames) > 1) {
    # check number of rows and columns
    nr <- lapply(lapply(env$covnames, function(x) env[[x]]), 
        function(y) sapply(y, function(z) nrow(as.matrix(z))))
    nr <- do.call(rbind, nr)
    nc <- lapply(lapply(env$covnames, function(x) env[[x]]), 
        function(y) sapply(y, function(z) ncol(as.matrix(z))))
    nc <- do.call(rbind, nc)
    for (i in 1:ncol(nr)) {
      if (length(unique(nr[, i])) > 1) {
        env$auto.adjust <- TRUE
      }
    }
    for (i in 1:ncol(nc)) {
      if (length(unique(nc[, i])) > 1) {
        env$auto.adjust <- TRUE
      }
    }
    if (verbose == TRUE && env$auto.adjust == TRUE) {
       message(paste("\nDimensions differ across networks within time steps."))
    }
    # check if labels are present
    if (env$auto.adjust == TRUE) {
      for (i in 1:length(env$covnames)) {
        for (t in 1:env$time.steps) {
          if (is.null(rownames(as.matrix(env[[env$covnames[i]]][[t]]))) || 
              is.null(colnames(as.matrix(env[[env$covnames[i]]][[t]])))) {
            stop(paste0("The dimensions of the covariates differ, but ", 
                "covariate '", env$covnames[i], 
                " does not have node labels at t = ", t, 
                ". Automatic adjustment of dimensions is therefore not ", 
                "possible."))
          }
        }
      }
    }
    # check if there are different labels despite identical dimensions
    if (env$auto.adjust == FALSE) {
      for (t in 1:env$time.steps) {
        rlabels.i <- list()
        clabels.i <- list()
        for (i in 1:length(env$covnames)) {
          rlabels.i[[i]] <- rownames(as.matrix(env[[env$covnames[i]]][[t]]))
          clabels.i[[i]] <- colnames(as.matrix(env[[env$covnames[i]]][[t]]))
        }
        rlabels.i <- do.call(rbind, rlabels.i)
        clabels.i <- do.call(rbind, clabels.i)
        flag <- FALSE
        if (!is.null(rlabels.i)) {
          for (j in 1:ncol(rlabels.i)) {
            if (length(unique(rlabels.i[, j])) > 1) {
              env$auto.adjust <- TRUE
              flag <- TRUE
              break
            }
          }
        }
        if (!is.null(clabels.i)) {
          for (j in 1:ncol(clabels.i)) {
            if (length(unique(clabels.i[, j])) > 1) {
              env$auto.adjust <- TRUE
              flag <- TRUE
              break
            }
          }
        }
      }
      if (verbose == TRUE && flag == TRUE) {
        message(paste("\nSame dimensions but different labels across", 
            "networks within time steps."))
      }
    }
  }
  if (verbose == TRUE && env$auto.adjust == TRUE) {
    message("Trying to auto-adjust the dimensions of the networks. ", 
        "If this fails, provide conformable matrices or network objects.")
  } else if (verbose == TRUE) {
    message("\nAll networks are conformable.")
  }
  
  # do mutual adjustment of networks and covariates at each time step
  structzero.df <- data.frame(label = character(), time = integer(), 
      object = character(), where = character())
  if (length(env$covnames) > 0 && env$auto.adjust == TRUE) {
    for (i in 1:env$time.steps) {
      for (j in 1:length(env$covnames)) {
        for (k in 1:length(env$covnames)) {
          if (j != k) {
            nw.j <- env[[env$covnames[j]]][[i]]
            rn.j <- rownames(as.matrix(nw.j))
            cn.j <- colnames(as.matrix(nw.j))
            nr.j <- nrow(as.matrix(nw.j))
            nc.j <- ncol(as.matrix(nw.j))
            nw.k <- env[[env$covnames[k]]][[i]]
            rn.k <- rownames(as.matrix(nw.k))
            cn.k <- colnames(as.matrix(nw.k))
            nr.k <- nrow(as.matrix(nw.k))
            nc.k <- ncol(as.matrix(nw.k))
            if (is.null(rn.j) || is.null(rn.k)) {
              stop(paste0("Missing row or column labels in object '", 
                  env$covnames[j], "'. Provide row and column ", 
                  "labels for all networks and covariates."))
            } else if (is.null(cn.j) || is.null(cn.k)) {
              stop(paste0("Missing row or column labels in object '", 
                  env$covnames[k], "'. Provide row and column ", 
                  "labels for all networks and covariates."))
            } else {
              if (is.null(rn.j) && !is.null(rn.k) && nr.j == nr.k) {
                if (class(nw.j) == "network") {
                  network::set.vertex.attribute(nw.j, "vertex.names", rn.k)
                } else {
                  rownames(nw.j) <- rn.k
                }
              } else if (is.null(rn.k) && !is.null(rn.j) && nr.j == nr.k) {
                if (class(nw.k) == "network") {
                  network::set.vertex.attribute(nw.k, "vertex.names", rn.j)
                } else {
                  rownames(nw.k) <- rn.j
                }
              } else if ((is.null(rn.k) || is.null(rn.j)) && nr.j != nr.k) {
                stop(paste0("Object '", env$covnames[j], 
                    "' is incompatible with object '", env$covnames[k], 
                    "' at t = ", i, "."))
              }
              # adjust j to k
              nw.j.labels <- adjust(nw.j, nw.k, remove = FALSE, 
                  value = 1, returnlabels = TRUE)
              nw.j <- adjust(nw.j, nw.k, remove = FALSE, value = 1)
              env[[env$covnames[j]]][[i]] <- nw.j
              ro <- nw.j.labels$added.row
              co <- nw.j.labels$added.col
              if (length(ro) > 0) {
                ro <- data.frame(label = ro, time = rep(i, length(ro)), 
                    object = rep(env$covnames[j], length(ro)), 
                    where = rep("row", length(ro)))
                structzero.df <- rbind(structzero.df, ro)
              }
              if (length(co) > 0) {
                co <- data.frame(label = co, time = rep(i, length(co)), 
                    object = rep(env$covnames[j], length(co)), 
                    where = rep("col", length(co)))
                structzero.df <- rbind(structzero.df, co)
              }
              # adjust k back to j
              nw.k.labels <- adjust(nw.j, nw.k, remove = FALSE, 
                  value = 1, returnlabels = TRUE)
              nw.k <- adjust(nw.k, nw.j, remove = FALSE, value = 1)
              env[[env$covnames[k]]][[i]] <- nw.k
              ro <- nw.k.labels$added.row
              co <- nw.k.labels$added.col
              if (length(ro) > 0) {
                ro <- data.frame(label = ro, time = rep(i, length(ro)), 
                    object = rep(env$covnames[j], length(ro)), 
                    where = rep("row", length(ro)))
                structzero.df <- rbind(structzero.df, ro)
              }
              if (length(co) > 0) {
                co <- data.frame(label = co, time = rep(i, length(co)), 
                    object = rep(env$covnames[j], length(co)), 
                    where = rep("col", length(co)))
                structzero.df <- rbind(structzero.df, co)
              }
            }
          }
        }
      }
    }
  }
  
  # check whether all dimensions are cross-sectionally conformable now
  nr.net <- sapply(env$networks, function(x) nrow(as.matrix(x)))
  for (i in 1:length(env$covnames)) {
    nr <- sapply(env[[env$covnames[i]]], function(x) {
      nrow(as.matrix(x))
    })
    for (j in 1:env$time.steps) {
      if (nr[j] != nr.net[j]) {
        stop(paste0("Covariate object '", env$covnames[i], 
            "' does not have the same number of rows as the dependent ", 
            "network at time step ", j, "."))
      }
    }
  }
  nc.net <- sapply(env$networks, function(x) ncol(as.matrix(x)))
  for (i in 1:length(env$covnames)) {
    nc <- sapply(env[[env$covnames[i]]], function(x) {
      ncol(as.matrix(x))
    })
    for (j in 1:env$time.steps) {
      if (nc[j] != nc.net[j]) {
        stop(paste0("Covariate object '", env$covnames[i], 
            "' does not have the same number of columns as the dependent ", 
            "network at time step ", j, "."))
      }
    }
  }
  
  # reporting
  if (verbose == TRUE) {
    if (env$auto.adjust == TRUE) {
      sz.row <- unique(structzero.df[structzero.df$where == "row", -3])
      szrownum <- numeric(length(env$networks))
      for (i in 1:length(env$networks)) {
        szrownum[i] <- nrow(sz.row[sz.row$time == i, ])
      }
      sz.col <- unique(structzero.df[structzero.df$where == "col", -3])
      szcolnum <- numeric(length(env$networks))
      for (i in 1:length(env$networks)) {
        szcolnum[i] <- nrow(sz.col[sz.col$time == i, ])
      }
      totrow <- sapply(env$networks, function(x) nrow(as.matrix(x)))
      totcol <- sapply(env$networks, function(x) ncol(as.matrix(x)))
      if (offset == TRUE) {
        dimensions <- rbind(totrow, totcol, szrownum, szcolnum, 
            totrow - szrownum, totcol - szcolnum)
        rownames(dimensions) <- c("total number of rows", 
            "total number of columns", "row-wise structural zeros", 
            "column-wise structural zeros", "remaining rows", 
            "remaining columns")
      } else {
        dimensions <- rbind(szrownum, szcolnum, totrow - szrownum, 
            totcol - szcolnum)
        rownames(dimensions) <- c("maximum deleted nodes (row)", 
            "maximum deleted nodes (col)", "remaining rows", 
            "remaining columns")
      }
      colnames(dimensions) <- paste0("t=", 1:length(env$networks))
      message("\nNumber of nodes per time step after adjustment:")
      print(dimensions)
      if (nrow(structzero.df) > 0) {
        if (offset == TRUE) {
          message("\nNodes affected completely by structural zeros:")
        } else {
          message("\nNodes removed completely during adjustment:")
        }
        print(unique(structzero.df))
      } else {
        message("\nAll nodes are retained.")
      }
    }
  }
  
  # create list of offset matrices (required both for offset and node removal)
  env$offsmat <- list()
  for (i in 1:env$time.steps) {
    mat <- matrix(0, nrow = nrow(as.matrix(env$networks[[i]])), 
        ncol = ncol(as.matrix(env$networks[[i]])))
    rownames(mat) <- rownames(as.matrix(env$networks[[i]]))
    colnames(mat) <- colnames(as.matrix(env$networks[[i]]))
    env$offsmat[[i]] <- mat
  }
  if (nrow(structzero.df) > 0) {
    for (i in 1:nrow(structzero.df)) {
      if (structzero.df$where[i] == "row") {
        index <- which(rownames(env$offsmat[[structzero.df$time[i]]]) == 
            structzero.df$label[i])
        env$offsmat[[structzero.df$time[i]]][index, ] <- 1
      } else {
        index <- which(colnames(env$offsmat[[structzero.df$time[i]]]) == 
            structzero.df$label[i])
        env$offsmat[[structzero.df$time[i]]][, index] <- 1
      }
    }
  }
  
  # offset preparation or node removal for MPLE
  if (offset == TRUE) {
    # add offset to formula and reassemble formula
    env$rhs.terms[length(env$rhs.terms) + 1] <- "offset(edgecov(offsmat[[i]]))"
    rhs.operators[length(rhs.operators) + 1] <- "+"
  } else {
    # delete nodes with structural zeros
    if (env$auto.adjust == TRUE) {
      env$offsmat <- suppressMessages(handleMissings(env$offsmat, na = 1, 
          method = "remove"))
      for (j in 1:length(env$covnames)) {
        env[[env$covnames[j]]] <- adjust(env[[env$covnames[j]]], env$offsmat)
      }
    }
  }
  
  # assemble formula
  rhs <- env$rhs.terms[1]
  if (length(rhs.operators) > 0) {
    for (i in 1:length(rhs.operators)) {
      rhs <- paste(rhs, rhs.operators[i], env$rhs.terms[i + 1])
    }
  }
  f <- paste(lhs, tilde, rhs)
  env$form <- as.formula(f, env = env)
  
  # for mtergm estimation using MCMC: create block-diagonal matrices
  if (blockdiag == TRUE) {
    if (env$bipartite == TRUE) {
      stop(paste("MCMC estimation is currently only supported for one-mode", 
          "networks. Use the btergm function instead."))
    }
    # use original formula without time indices
    env$form <- update.formula(formula, networks ~ .)
    env$form <- paste(deparse(env$form), collapse = "")
    env$form <- paste(env$form, "+ offset(edgecov(offsmat))")
    env$form <- as.formula(env$form, env = env)
    # make covariates block-diagonal
    if (length(env$covnames) > 1) {
      for (i in 2:length(env$covnames)) {
        env[[env$covnames[i]]] <- as.matrix(Matrix::bdiag(as.matrix(
            env[[env$covnames[i]]])))
      }
    }
    # create block-diagonal offset matrix and merge with existing offsmat term
    env$offsmat <- as.matrix(Matrix::bdiag(env$offsmat))  # make block-diagonal
    bdoffset <- lapply(env$networks, as.matrix)
    for (i in 1:length(bdoffset)) {
      bdoffset[[i]][, ] <- 1
    }
    bdoffset <- as.matrix((Matrix::bdiag(bdoffset) - 1) * -1)  # off-diagonal
    env$offsmat <- env$offsmat + bdoffset
    rm(bdoffset)
    env$offsmat[env$offsmat > 0] <- 1
    # make dependent network block-diagonal
    if (class(env$networks[[1]]) == "network") {  # network
      attrnames <- network::list.vertex.attributes(env$networks[[1]])
      attributes <- list()
      for (i in 1:length(env$networks)) {
        attrib <- list()
        for (j in 1:length(attrnames)) {
          attrib[[j]] <- network::get.vertex.attribute(env$networks[[i]], 
              attrnames[j])
        }
        attributes[[i]] <- attrib
        env$networks[[i]] <- as.matrix(env$networks[[i]])
      }
      env$networks <- network::network(as.matrix(Matrix::bdiag(env$networks)), 
          directed = env$directed, bipartite = env$bipartite)
      for (i in 1:length(attrnames)) {  # collate attributes and merge back in
        attrib <- unlist(lapply(attributes, function(x) x[[i]]))
        network::set.vertex.attribute(env$networks, attrnames[i], attrib)
      }
    } else {  # matrix
      env$networks <- network::network(as.matrix(Matrix::bdiag(env$networks)), 
          directed = env$directed, bipartite = env$bipartite)
    }
    if (verbose == TRUE) {
      cat("\n")  # to get a blank line before the MCMC MLE output starts
    }
  }
  
  return(env)  # return the environment with all the data
}


# TERGM by bootstrapped pseudolikelihood
btergm <- function(formula, R = 500, offset = FALSE, parallel = c("no", 
    "multicore", "snow"), ncpus = 1, cl = NULL, verbose = TRUE, ...) {
  
  # call tergmprepare and integrate results as a child environment in the chain
  env <- tergmprepare(formula = formula, offset = offset, verbose = verbose)
  parent.env(env) <- environment()
  
  # check number of time steps
  if (env$time.steps == 1) {
    warning(paste("The confidence intervals and standard errors are",
        "meaningless because only one time step was provided."))
  }
  
  # verbose reporting
  if (verbose == TRUE) {
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
  }
  
  # create MPLE data structures
  if (offset == TRUE) {  # via structural zeros and an offset term
    # create the data for MPLE with structural zeros
    Y <- NULL  # dependent variable
    X <- NULL  # independent variables data frame
    W <- NULL  # weights for each observation
    O <- NULL  # offset term
    for (i in 1:length(env$networks)) {
      nw <- ergm::ergm.getnetwork(env$form)
      model <- ergm::ergm.getmodel(env$form, nw, initialfit = TRUE)
      Clist <- ergm::ergm.Cprepare(nw, model)
      Clist.miss <- ergm::ergm.design(nw, model, verbose = FALSE)
      pl <- ergm::ergm.pl(Clist, Clist.miss, model, theta.offset = 
          c(rep(FALSE, length(env$rhs.terms) - 1), TRUE), verbose = FALSE, 
          control = ergm::control.ergm(init = c(rep(NA, 
          length(env$rhs.terms) - 1), 1)))
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
    for (i in 1:length(env$networks)) {
      mpli <- ergm::ergmMPLE(env$form)
      Y <- c(Y, mpli$response)
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
      parallel = parallel, ncpus = ncpus, cl = cl, ...)$t
  rm(X)
  if (nrow(coefs) == 1) { # in case there is only one model term
    coefs <- t(coefs)
  }
  
  # create and return btergm object
  colnames(coefs) <- term.names[1:(length(term.names) - 1)]
  names(startval) <- colnames(coefs)
  
  btergm.object <- createBtergm(startval, coefs, R, nobs, env$time.steps, 
      formula, Y, x, W, env$auto.adjust, offset, env$directed, env$bipartite)
  if (verbose == TRUE) {
    message("Done.")
  }
  return(btergm.object)
}


# MCMC MLE estimation function (basically a wrapper for the ergm function)
mtergm <- function(formula, offset = FALSE, constraints = ~ ., 
    estimate = c("MLE", "MPLE"), verbose = TRUE, ...) {
  
  # call tergmprepare and integrate results as a child environment in the chain
  env <- tergmprepare(formula = formula, offset = offset, blockdiag = TRUE, 
      verbose = verbose)
  parent.env(env) <- environment()
  
  # estimate an ERGM
  e <- ergm(env$form, offset.coef = -Inf, constraints = constraints, 
      estimate = estimate, ...)
  
  if (verbose == TRUE) {
    message(paste("\nNote: The infinite 'edgecov.offsmat' model term contains", 
        "structural zeros and can be ignored."))
  }
  
#  # remove information about the offset term from the ergm object (here: MPLE)
#  e$coef <- e$coef[-length(e$coef)]
#  e$MCMCtheta <- e$MCMCtheta[-length(e$MCMCtheta)]
#  e$gradient <- e$gradient[-length(e$gradient)]
#  e$covar <- e$covar[, -ncol(e$covar)]
#  e$mc.se <- e$mc.se[-length(e$mc.se)]
#  e$offset <- e$offset[-length(e$offset)]
#  e$drop <- e$drop[-length(e$drop)]
#  e$estimable <- e$estimable[-length(e$estimable)]
#  e$offset <- e$offset[-length(e$offset)]
#  e$formula <- f
#  e$control$init <- e$control$init[-length(e$control$init)]
  
  return(e)
}


# simulation of new networks based on a btergm fit
simulate.btergm <- function(object, nsim = 1, seed = NULL, index = NULL, 
    formula = getformula(object), coef = object@coef, verbose = TRUE, ...) {
  
  # call tergmprepare and integrate results as a child environment in the chain
  env <- tergmprepare(formula = formula, offset = object@offset, 
      verbose = FALSE)
  parent.env(env) <- environment()
  
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
        paste(deparse(env$form), collapse = ""))
    f.i <- gsub("\\s+", " ", f.i)
    f.i <- gsub("^networks", env$lhs.original, f.i)
    message(paste("Simulating", nsim, "networks from the following formula:\n", 
        f.i, "\n"))
  }
  
  # simulate
  ergm::simulate.formula(env$form, nsim = nsim, seed = seed, coef = coef, 
      verbose = verbose, ...)
}

