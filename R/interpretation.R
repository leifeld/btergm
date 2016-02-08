# interpretation function for ergm objects
interpret.ergm <- function(object, formula = object$formula, 
    coefficients = coef(object), target = eval(parse(text = 
    deparse(formula[[2]]))), type = "tie", i, j, ...) {
  
  if (length(j) > 1 && (type == "tie" || type == "dyad")) {
    stop(paste("For computing dyadic or tie probabilities, only a single 'j'",
        "node can be specified."))
  }
  
  # extract response network and adjust formula
  nw <- target
  dir <- is.directed(nw)
  if (dir == FALSE && type == "dyad") {
    type <- "tie"
    warning(paste("Dyadic probabilities not available for undirected",
        "networks. Reporting tie probability instead."))
  }
  
  # disassemble formula and preprocess rhs
  tilde <- deparse(formula[[1]])
  lhs <- "nw"
  rhs <- paste(deparse(formula[[3]]), collapse = "")
  rhs <- gsub("\\s+", " ", rhs)
  if (grepl("(gw.{0,2}degree)|(gwdsp)|(gwesp)|(gwnsp)", rhs)) {
    stop("The interpretation functions do not work with curved model terms.")
  }
  
  # reassemble formula
  f <- paste(lhs, tilde, rhs)
  form <- as.formula(f)
  
  if (type == "tie") {
    nw[i, j] <- 0  # set zero and compute stats
    mod <- ergm.getmodel(form, nw, drop = FALSE)
    stat0 <- ergm.getglobalstats(nw, mod)
    
    nw[i, j] <- 1  # set one and compute stats
    mod <- ergm.getmodel(form, nw, drop = FALSE)
    stat1 <- ergm.getglobalstats(nw, mod)
    
    # compute change statistics and ratio
    chgstat <- stat1 - stat0
    if (length(chgstat) != length(coefficients)) {
      stop(paste("Number of coefficients and statistics differ.",
          "Did you fit a curved model? Curved models with non-fixed",
          "parameters are currently not supported."))
    }
    lp <- t(chgstat) %*% cbind(coefficients)
    results <- c(1 / (1 + exp(-lp)))
    names(results) <- "i->j = 1"
  } else if (type == "dyad") {
    eta_mat <- matrix(NA, 2, 2)
    for (xi in 0:1) {
      for (xj in 0:1) {
        nw[i, j] <- xi
        nw[j, i] <- xj
        mod <- ergm.getmodel(form, nw, drop = FALSE)
        stat <- ergm.getglobalstats(ergm.getnetwork(form), mod)
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
    results <- prob_mat
  } else if (type == "node") {
    m <- length(i)
    n <- length(j)
    if (m == 1 && n > 1) {
      labels <- c("Sender", "Receiver")
    } else if (m > 1 && n == 1) {
      labels <- c("Receiver", "Sender")
      j.old <- j
      j <- i
      i <- j.old
      m <- length(i)
      n <- length(j)
    } else {
      stop("Either 'i' or 'j' must contain more than one node.")
    }
    vecs <- rbind(rep(0, n), rep(1, n))
    base <- rep(0, n)
    for (l in 1:(n - 1)) {
      places <- t(combn(1:n, l))
      for (k in 1:nrow(places)) {
        veci <- base
        veci[places[k, ]] <- 1
        vecs <- rbind(vecs, veci)
      }
    }
    eta <- numeric(nrow(vecs))
    for (l in 1:nrow(vecs)) {
      nw[i, j] <- vecs[l, ]
      mod <- ergm.getmodel(form, nw, drop = FALSE)
      stat <- ergm.getglobalstats(ergm.getnetwork(form), mod)
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
    colnames(vecs) <- paste(labels[2], j)
    rownames(vecs) <- rep(paste(labels[1], i), nrow(vecs))
    results <- cbind(prob, vecs)
    colnames(results)[1] <- "probability"
  } else {
    stop("'type' argument undefined or not recognized.")
  }
  return(results)
}


# interpretation method for btergm objects
interpret.btergm <- function(object, formula = getformula(object), 
    coefficients = coef(object), target = NULL, type = "tie", i, j, 
    t = 1:object@time.steps, ...) {
  
  env <- tergmprepare(formula = formula, offset = FALSE, blockdiag = FALSE, 
      verbose = FALSE)
  
  # extract response networks and adjust formula
  if (is.null(target)) {
    networks <- env$networks
  } else {
    networks <- target
  }
  if (class(networks) == "network" || class(networks) == "matrix") {
    networks <- list(networks)
  }
  dir <- is.directed(networks[[1]])
  if (dir == FALSE && type == "dyad") {
    type <- "tie"
    warning(paste("Dyadic probabilities not available for undirected",
        "networks. Reporting tie probabilities instead."))
  }
  if (is.null(t)) {
    t <- 1:length(networks)
  }
  
  # prepare i and j
  if (!is.list(i)) {
    i <- rep(list(i), length(networks))
    num.actors <- numeric()
    for (k in t) {
      num.actors[k] <- nrow(as.matrix(networks[[k]]))
    }
    if (length(table(num.actors)) > 1) {
      warning(paste("'i' does not vary across time steps, but the number of",
          "actors does. 'i' can be provided as a list."))
    }
  }
  if (!is.list(j)) {
    j <- rep(list(j), length(networks))
    num.actors <- numeric()
    for (k in t) {
      num.actors[k] <- nrow(as.matrix(networks[[k]]))
    }
    if (length(table(num.actors)) > 1) {
      warning(paste("'j' does not vary across time steps, but the number of",
          "actors does. 'j' can be provided as a list."))
    }
  }
  for (l in 1:length(j)) {
    if (length(j[[l]]) > 1 && (type == "tie" || type == "dyad")) {
      stop(paste("For computing dyadic or tie probabilities, only a single 'j'",
          "node can be specified per time step."))
    }
  }
  
  # do the computations
  if (type == "tie") {
    results <- numeric()
    for (k in t) {
      networks[[k]][i[[k]], j[[k]]] <- 0
      stat0 <- summary(remove.offset.formula(env$form), response = NULL)
      networks[[k]][i[[k]], j[[k]]] <- 1
      stat1 <- summary(remove.offset.formula(env$form), response = NULL)
      chgstat <- stat1 - stat0
      if (length(chgstat) != length(coefficients)) {
        stop(paste("Number of coefficients and statistics differ.",
            "Did you fit a curved model? Curved models with non-fixed",
            "parameters are currently not supported."))
      }
      lp <- t(chgstat) %*% cbind(coefficients)
      result <- c(1 / (1 + exp(-lp)))
      names(result) <- "i->j = 1"
      results[k] <- result
    }
    results <- results[!is.na(results)]
    names(results) <- paste("t =", t)
  } else if (type == "dyad") {
    results <- list()
    for (k in t) {
      eta_mat <- matrix(NA, 2, 2)
      for (xi in 0:1) {
        for (xj in 0:1) {
          networks[[k]][i[[k]], j[[k]]] <- xi
          networks[[k]][j[[k]], i[[k]]] <- xj
          stat <- summary(remove.offset.formula(env$form), response = NULL)
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
      results[[k]] <- prob_mat
    }
    results <- results[!sapply(results, is.null)]
    names(results) <- paste("t =", t)
  } else if (type == "node") {
    results <- list()
    for (k in t) {
      m <- length(i[[k]])
      n <- length(j[[k]])
      if (m == 1 && n > 1) {
        labels <- c("Sender", "Receiver")
      } else if (m > 1 && n == 1) {
        labels <- c("Receiver", "Sender")
        j.old <- j[[k]]
        j[[k]] <- i[[k]]
        i[[k]] <- j.old
        m <- length(i[[k]])
        n <- length(j[[k]])
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
        ik <- i[[k]]
        jk <- j[[k]]
        networks[[k]][ik, jk] <- vecs[l, ]
        stat <- summary(remove.offset.formula(env$form), response = NULL)
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
      colnames(vecs) <- paste(labels[2], j[[k]])
      rownames(vecs) <- rep(paste(labels[1], i[[k]]), nrow(vecs))
      result <- cbind(prob, vecs)
      colnames(result)[1] <- "probability"
      results[[k]] <- result
    }
    results <- results[!sapply(results, is.null)]
    names(results) <- paste("t =", t)
  } else {
    stop("'type' argument undefined or not recognized.")
  }
  return(results)
}


# register generic methods with ergm and btergm objects
setMethod("interpret", signature = className("ergm", "ergm"), 
    definition = interpret.ergm)

setMethod("interpret", signature = className("btergm", "btergm"), 
    definition = interpret.btergm)

setMethod("interpret", signature = className("mtergm", "btergm"), 
    definition = interpret.btergm)
