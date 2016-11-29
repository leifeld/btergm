
# checkdegeneracy function for mtergm objects
checkdegeneracy.mtergm <- function(object, ...) {
  mcmc.diagnostics(object@ergm, ...)
}

setMethod("checkdegeneracy", signature = className("mtergm", "btergm"), 
    definition = checkdegeneracy.mtergm)


# checkdegeneracy function for btergm objects
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
  if (class(object)[1] == "btergm" && offset == TRUE) {
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
    target.stats[[index]] <- summary(ergm::remove.offset.formula(form), 
        response = NULL)
    degen[[index]] <- simulate.formula(form, nsim = nsim, 
        coef = coefs, statsonly = TRUE, 
        control = control.simulate.formula(MCMC.interval = 
        MCMC.interval, MCMC.burnin = MCMC.burnin))
    if (offset == TRUE || "mtergm" %in% class(object)) {
      degen[[i]] <- degen[[i]][, -ncol(degen[[i]])]  # remove offset statistic
    }
  }
  
  if (verbose == TRUE) {
    message("Checking degeneracy...")
  }
  class(mat) <- "degeneracy"
  object <- list()
  object$target.stats <- target.stats
  object$sim <- degen
  class(object) <- "degeneracy"
  if (verbose == TRUE) {
    message("Done.")
  }
  return(object)
}

setMethod("checkdegeneracy", signature = className("btergm", "btergm"), 
    definition = checkdegeneracy.btergm)


# print method for 'degeneracy' objects
print.degeneracy <- function(x, center = FALSE, t = 1:length(x$sim), 
    terms = 1:length(x$target.stats[[1]]), ...) {
  for (i in t) {
    message(paste0("\nDegeneracy check for network ", i, ":"))
    if (center == TRUE) {
      sm <- coda::as.mcmc.list(coda::as.mcmc(x$sim[[i]]))
      sm <- ergm::sweep.mcmc.list(sm, x$target.stats[[i]], "-")[[1]] # diff
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
    if (class(q) != "matrix") {
      q <- matrix(q, nrow = 1)
      rownames(q) <- rn
      colnames(q) <- cn
    }
    printCoefmat(q)
  }
}


# plot method for 'degeneracy' objects
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