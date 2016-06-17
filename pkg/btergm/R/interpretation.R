
# Internal helper function, which does the actual interpretation computations. 
# Needs an environment as provided by tergmprepare, including a list of 
# networks (env$networks, with possibly just one network in this list) and a 
# formula (env$form, with temporal indices, as prepared by tergmprepare). Also 
# need to supply coefficients vector, type string, i, j, and t (where t can be 
# 1 if there is only one network).
dointerpret <- function(env, coefficients, type, i, j, t) {
  parent.env(env) <- environment()
  
  # prepare i and j
  if (!is.list(i)) {
    i <- rep(list(i), length(env$networks))
    num.actors <- numeric()
    for (k in t) {  # note that t can be a vector of time steps
      num.actors[k] <- nrow(as.matrix(env$networks[[k]]))
    }
    if (length(table(num.actors)) > 1) {
      warning(paste("'i' does not vary across time steps, but the number of",
          "actors does. 'i' can be provided as a list or as a name."))
    }
  }
  if (!is.list(j)) {
    j <- rep(list(j), length(env$networks))
    num.actors <- numeric()
    for (k in t) {
      num.actors[k] <- nrow(as.matrix(env$networks[[k]]))
    }
    if (length(table(num.actors)) > 1) {
      warning(paste("'j' does not vary across time steps, but the number of",
          "actors does. 'j' can be provided as a list or as a name."))
    }
  }
  for (l in 1:length(j)) {
    if (length(j[[l]]) > 1 && (type == "tie" || type == "dyad")) {
      stop(paste("For computing dyadic or tie probabilities, only a single 'j'",
          "node can be specified per time step."))
    }
  }
  node_i <- i
  node_j <- j
  
  if (type == "tie") {
    results <- numeric()
    for (i in t) {
      env$networks[[i]][node_i[[i]], node_j[[i]]] <- 0
      stat0 <- summary(env$form, response = NULL)
      env$networks[[i]][node_i[[i]], node_j[[i]]] <- 1
      stat1 <- summary(env$form, response = NULL)
      chgstat <- stat1 - stat0
      if (length(chgstat) != length(coefficients)) {
        stop(paste("Number of coefficients and statistics differ.",
            "Did you fit a curved model? Curved models with non-fixed",
            "parameters are currently not supported."))
      }
      for (k in 1:length(coefficients)) {
        if (coefficients[k] == -Inf && chgstat[k] == 0) {
          coefficients[k] <- 0  # no struct zero dyad; replace -Inf to avoid NaN
        }
      }
      lp <- t(chgstat) %*% cbind(coefficients)
      result <- c(1 / (1 + exp(-lp)))
      names(result) <- "i->j = 1"
      results[i] <- result
    }
    results <- results[!is.na(results)]
    names(results) <- paste("t =", t)
  } else if (type == "dyad") {
    if (-Inf %in% coefficients) {
      coefficients[coefficients == -Inf] <- -9e8
      warning(paste("There are -Inf coefficients (possibly due to", 
          "offset terms). To yield interpretable results, -Inf is", 
          "approximated by a large negative number (-9e8). Note that", 
          "this may be imprecise."))
    }
    results <- list()
    for (i in t) {
      # print error if undirected
      if ((class(env$networks[[i]]) == "network" && 
          !is.directed(env$networks[[i]])) || 
          (class(env$networks[[i]]) == "matrix" && 
          is.mat.directed(env$networks[[i]]) == FALSE)) {
        stop(paste0("Network at t=", i, " is undirected. Dyadic ", 
            "probabilities do not make sense in undirected networks. Try ", 
            "type = \"tie\" instead!"))
      }
      eta_mat <- matrix(NA, 2, 2)
      for (xi in 0:1) {
        for (xj in 0:1) {
          env$networks[[i]][node_i[[i]], node_j[[i]]] <- xi
          env$networks[[i]][node_j[[i]], node_i[[i]]] <- xj
          stat <- summary(env$form, response = NULL)
          if (length(stat) != length(coefficients)) {
            stop(paste("Number of coefficients and statistics differ.",
                "Did you fit a curved model? Curved models with non-fixed",
                "parameters are currently not supported."))
          }
          eta_mat[xi + 1, xj + 1] <- t(coefficients) %*% cbind(stat)
        }
      }
      prob_mat <- matrix(NA, 2, 2)
      for (xi in 0:1) {
        for (xj in 0:1) {
          etas <- c(eta_mat) - eta_mat[xi + 1, xj + 1]
          prob_mat[xi + 1, xj + 1] <- 1 / (sum(exp(etas)))
        }
      }
      rownames(prob_mat) <- c("i->j = 0", "i->j = 1")
      colnames(prob_mat) <- c("j->i = 0", "j->i = 1")
      results[[i]] <- prob_mat
    }
    results <- results[!sapply(results, is.null)]
    names(results) <- paste("t =", t)
  } else if (type == "node") {
    if (-Inf %in% coefficients) {
      coefficients[coefficients == -Inf] <- -9e8
      warning(paste("There are -Inf coefficients (possibly due to", 
          "offset terms). To yield interpretable results, -Inf is", 
          "approximated by a large negative number (-9e8). Note that", 
          "this may be imprecise."))
    }
    results <- list()
    for (i in t) {
      m <- length(node_i[[i]])
      n <- length(node_j[[i]])
      if (m == 1 && n > 1) {
        labels <- c("Sender", "Receiver")
      } else if (m > 1 && n == 1) {
        labels <- c("Receiver", "Sender")
        j.old <- node_j[[i]]
        node_j[[i]] <- node_i[[i]]
        node_i[[i]] <- j.old
        m <- length(node_i[[i]])
        n <- length(node_j[[i]])
      } else {
        stop(paste("Either 'i' or 'j' must contain more than one node per",
            "time step."))
      }
      vecs <- rbind(rep(0, n), rep(1, n))
      base <- rep(0, n)
      for (l in 1:(n - 1)) {
        places <- t(combn(1:n, l))
        for (r in 1:nrow(places)) {
          veci <- base
          veci[places[r, ]] <- 1
          vecs <- rbind(vecs, veci)
        }
      }
      eta <- numeric(nrow(vecs))
      for (l in 1:nrow(vecs)) {
        ik <- node_i[[i]]
        jk <- node_j[[i]]
        env$networks[[i]][ik, jk] <- vecs[l, ]
        stat <- summary(env$form, response = NULL)
        if (length(stat) != length(coefficients)) {
          stop(paste("Number of coefficients and statistics differ.",
              "Did you fit a curved model? Curved models with non-fixed",
              "parameters are currently not supported."))
        }
        eta[l] <- t(coefficients) %*% cbind(stat)
      }
      prob <- numeric(nrow(vecs))
      for (l in 1:nrow(vecs)) {
        prob[l] <- 1 / sum(exp(eta - eta[l]))
      }
      colnames(vecs) <- paste(labels[2], node_j[[i]])
      rownames(vecs) <- rep(paste(labels[1], node_i[[i]]), nrow(vecs))
      result <- cbind(prob, vecs)
      colnames(result)[1] <- "probability"
      results[[i]] <- result
    }
    results <- results[!sapply(results, is.null)]
    names(results) <- paste("t =", t)
  } else {
    stop("'type' argument undefined or not recognized.")
  }
  return(results)
}


# interpretation method for ergm objects
interpret.ergm <- function(object, formula = getformula(object), 
    coefficients = coef(object), target = NULL, type = "tie", i, j) {
  
  env <- tergmprepare(formula = formula, offset = FALSE, blockdiag = FALSE, 
      verbose = FALSE)
  
  # extract response network and adjust formula
  if (!is.null(target)) {
    env$networks <- list(target)
  }
  
  dointerpret(env, coefficients = coefficients, type = type, i = i, j = j, 
      t = 1)[[1]]
}


# interpretation method for btergm objects
interpret.btergm <- function(object, formula = getformula(object), 
    coefficients = coef(object), target = NULL, type = "tie", i, j, 
    t = 1:object@time.steps) {
  
  env <- tergmprepare(formula = formula, offset = FALSE, blockdiag = FALSE, 
      verbose = FALSE)
  parent.env(env) <- environment()
  
  # extract response networks and adjust formula
  if (!is.null(target)) {
    env$networks <- target
  }
  
  dointerpret(env, coefficients = coefficients, type = type, i = i, j = j, 
      t = t)
}


# register generic methods with ergm, btergm, and mtergm objects
setMethod("interpret", signature = className("ergm", "ergm"), 
    definition = interpret.ergm)

setMethod("interpret", signature = className("btergm", "btergm"), 
    definition = interpret.btergm)

setMethod("interpret", signature = className("mtergm", "btergm"), 
    definition = interpret.btergm)


# a function that creates all tie probabilities along with some other variables
edgeprob <- function(object, parallel = c("no", "multicore", "snow"), 
    ncpus = 1, cl = NULL, verbose = TRUE) {
  
  # determine if ERGM or TERGM
  if (class(object) == "ergm") {
    tergm <- FALSE
  } else if (class(object) %in% c("btergm", "mtergm")) {
    tergm <- TRUE
  } else {
    stop(paste("The tieprob function is only applicable to ergm, btergm, and",
        "mtergm objects."))
  }
  
  # prepare data structures in a local environment
  env <- tergmprepare(formula = getformula(object), offset = FALSE, 
      blockdiag = FALSE, verbose = FALSE)
  parent.env(env) <- environment()
  covnames <- env$covnames[-1]  # leave out the LHS network(s)
  coefs <- coef(object)
  
  # create matrix with MPLE predictors and i, j, and t indices
  if (verbose == TRUE) {
    message("Creating data frame with predictors...")
  }
  Y <- NULL
  dyads <- NULL
  for (i in 1:length(env$networks)) {
    mat <- as.matrix(env$networks[[i]])
    imat <- matrix(rep(1:nrow(mat), ncol(mat)), nrow = nrow(mat))
    if ((class(env$networks[[i]]) == "network" && 
        is.bipartite(env$networks[[i]])) || 
        (class(env$networks[[i]]) == "matrix" && 
        is.mat.onemode(env$networks[[i]]) == FALSE)) {
      mn <- nrow(mat) + 1
      mx <- nrow(mat) + ncol(mat)
      jmat <- matrix(rep(mn:mx, nrow(mat)), nrow = nrow(mat), byrow = TRUE)
    } else {
      jmat <- matrix(rep(1:ncol(mat), nrow(mat)), nrow = nrow(mat), 
          byrow = TRUE)
    }
    form.middle <- as.character(env$form)[[1]]
    form.left <- as.character(env$form)[[2]]
    form.right <- as.character(env$form)[[3]]
    form.right <- paste0(form.right, " + edgecov(imat) + edgecov(jmat)")
    f <- as.formula(paste(form.left, form.middle, form.right), env = env)
    mpli <- ergm::ergmMPLE(f)
    Y <- c(Y, mpli$response)
    dyads <- rbind(dyads, cbind(mpli$predictor, i))
  }
  term.names <- colnames(dyads)[-(length(colnames(dyads)):(length(colnames(
      dyads)) - 2))]
  term.names <- c(term.names, "i", "j", "t")
  dyads <- data.frame(dyads)
  colnames(dyads) <- term.names
  dyads <- cbind(Y, dyads)
  colnames(dyads)[1] <- "tie"
  class(dyads[, length(colnames(dyads))]) <- "integer"
  class(dyads[, length(colnames(dyads)) - 1]) <- "integer"
  class(dyads[, length(colnames(dyads)) - 2]) <- "integer"
  dyads <- dyads[sample(nrow(dyads)), ]
  
  # define interpret wrapper function for rows of the MPLE predictor matrix
  if (tergm == TRUE) {
    interpret.parallel <- function(index, subs, env, coefs) {
      dointerpret(env, coefs, type = "tie", i = subs[index, 1], 
          j = subs[index, 2], t = subs[index, 3])
    }
  } else {
    interpret.parallel <- function(index, subs, env, coefs) {
      dointerpret(env, coefs, type = "tie", i = subs[index, 1], 
          j = subs[index, 2], t = 1)
    }
  }
  
  # subset dyads if necessary
  indices <- 1:nrow(dyads)
  subs <- cbind(dyads$i, dyads$j, dyads$t)
  subs <- subs[indices, ]
  dyads <- dyads[indices, ]
  
  # compute probabilities (possibly in parallel)
  if (is.null(ncpus) || ncpus == 0) {
    ncpus <- 1
  }
  if (!parallel[1] %in% c("no", "multicore", "snow")) {
    parallel <- "no"
    warning("'parallel' argument not recognized. Using 'no' instead.")
  }
  if (parallel[1] == "snow") {
    if (is.null(cl)) {
      created <- TRUE
      cl <- makeCluster(ncpus)
    } else {
      created <- FALSE
    }
    if (verbose == TRUE) {
      message(paste("Using snow parallelization with", ncpus, "cores."))
    }
    clusterEvalQ(cl, library("btergm"))
    a <- proc.time()
    prob <- parSapply(cl, 1:nrow(subs), interpret.parallel, subs, env, coefs)
    b <- proc.time()
    if (created == TRUE) {
      stopCluster(cl)
    }
    if (verbose == TRUE) {
      message(paste("Time elapsed:"))
      print(b - a)
    }
  } else if (parallel[1] == "multicore") {
    if (verbose == TRUE) {
      message(paste("Using multicore parallelization with", ncpus, "cores."))
    }
    a <- proc.time()
    prob <- mclapply(1:nrow(subs), interpret.parallel, subs, env, coefs, 
        mc.cores = ncpus)
    prob <- unlist(prob)
    b <- proc.time()
    if (verbose == TRUE) {
      message(paste("Time elapsed:"))
      print(b - a)
    }
  } else {
    if (verbose == TRUE) {
      message(paste0("Parallelization is switched off. Computing", 
          "probabilities sequentially."))
    }
    a <- proc.time()
    prob <- sapply(1:nrow(subs), interpret.parallel, subs, env, coefs)
    b <- proc.time()
    if (verbose == TRUE) {
      message(paste("Time elapsed:"))
      print(b - a)
    }
  }
  
  # combine dyads/predictors and probabilities, sort, and return
  dyads <- cbind(dyads, prob)
  colnames(dyads)[ncol(dyads)] <- "probability"
  dyads <- dyads[order(dyads$t, dyads$i, dyads$j), ]
  rownames(dyads) <- NULL
  return(dyads)
}