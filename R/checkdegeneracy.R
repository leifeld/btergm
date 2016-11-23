
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
  mat <- list()
  for (i in 1:l$time.steps) {
    sm <- coda::as.mcmc.list(coda::as.mcmc(degen[[i]]))
    sm <- ergm::sweep.mcmc.list(sm, target.stats[[i]], "-")  # difference
    ds <- ergm::colMeans.mcmc.list(sm)  # mean difference for each statistic
    sds <- apply(degen[[i]], 2, sd)  # standard deviations of sample
    ns <- coda::effectiveSize(sm)  # instead of dividing by sample size...
    se <- sds / sqrt(ns)   # divide by effective sample size
    z <- ds / se  # divide means by SEs
    p.z <- pnorm(abs(z), lower.tail = FALSE) * 2  # p-values
    mat[[i]] <- cbind("obs" = target.stats[[i]], "sim" = colMeans(degen[[i]]), 
        "est" = ds, "se" = se, "zval" = z, "pval" = p.z)
  }
  class(mat) <- "degeneracy"
  if (verbose == TRUE) {
    message("Done.")
  }
  return(mat)
}

setMethod("checkdegeneracy", signature = className("btergm", "btergm"), 
    definition = checkdegeneracy.btergm)


# print method for 'degeneracy' objects
print.degeneracy <- function(x, ...) {
  for (i in 1:length(x)) {
    message(paste0("\nDegeneracy check for network ", i, ":"))
    printCoefmat(x[[i]], digits = 3, P.values = TRUE, has.Pvalue = TRUE, 
        cs.ind = 3:4, ts.ind = 5)
  }
  message("\nSmall p-values indicate degenerate results.")
}
