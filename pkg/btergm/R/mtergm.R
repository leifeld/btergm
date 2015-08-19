# MCMC MLE estimation function (basically a wrapper for the ergm function)
mtergm <- function(formula, offset = FALSE, constraints = ~ ., 
    estimate = c("MLE", "MPLE"), verbose = TRUE, ...) {
  
  # call tergmprepare and integrate results as a child environment in the chain
  env <- tergmprepare(formula = formula, offset = offset, blockdiag = TRUE, 
      verbose = verbose)
  parent.env(env) <- environment()
  
  if (verbose == TRUE) {
    message("Estimating...")
  }
  
  # estimate an ERGM
  e <- ergm(env$mtergmestform, offset.coef = -Inf, constraints = constraints, 
      estimate = estimate[1], ...)
  
  e$formula <- formula  # insert original formula
  class(e) <- c(class(e), "mtergm")  # add mtergm class label
  
  if (verbose == TRUE) {
    message("Done.")
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
