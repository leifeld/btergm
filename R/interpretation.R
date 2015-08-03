# interpretation function for ergm objects
interpret.ergm <- function(object, formula = object$formula, 
    coefficients = coef(object), network = eval(parse(text = 
    deparse(formula[[2]]))), type = "tie", i, j, ...) {
  
  if (length(j) > 1 && (type == "tie" || type == "dyad")) {
    stop(paste("For computing dyadic or tie probabilities, only a single 'j'",
        "node can be specified."))
  }
  
  # extract response network and adjust formula
  #nw <- eval(parse(text = deparse(formula[[2]])))  # response network
  nw <- network
  #nw <- ergm.getnetwork(formula)
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
interpret.btergm <- function(object, formula = object@formula, 
    coefficients = coef(object), network = eval(parse(text = 
    deparse(formula[[2]]))), type = "tie", i, j, t = 1:length(network), 
    ...) {
  
  # extract response networks and adjust formula
  networks <- network
  if (class(networks) == "network" || class(networks) == "matrix") {
    networks <- list(networks)
  }
  dir <- is.directed(networks[[1]])
  if (dir == FALSE && type == "dyad") {
    type <- "tie"
    warning(paste("Dyadic probabilities not available for undirected",
        "networks. Reporting tie probabilities instead."))
  }
  
  # disassemble formula and preprocess rhs
  tilde <- deparse(formula[[1]])
  lhs <- "networks[[k]]"
  rhs <- paste(deparse(formula[[3]]), collapse = "")
  rhs <- gsub("\\s+", " ", rhs)
  rhs <- preprocessrhs(rhs, length(networks), iterator = "k")
  
  # reassemble formula
  f <- paste(lhs, tilde, rhs)
  form <- as.formula(f)
  
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
      stat0 <- summary(remove.offset.formula(form), response = NULL)
      networks[[k]][i[[k]], j[[k]]] <- 1
      stat1 <- summary(remove.offset.formula(form), response = NULL)
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
      #names(results)[k] <- paste0("t = ", k)
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
          stat <- summary(remove.offset.formula(form), response = NULL)
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
        stat <- summary(remove.offset.formula(form), response = NULL)
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


# function which preprocesses the right-hand side (rhs) of the formula
# (may be obsolete in the near future when the tergmprepare function is used)
preprocessrhs <- function(rhs, time.steps, iterator = "i", dep = NULL) {
  
  covnames <- character()
  
  # split up rhs terms
  rhs.terms <- strsplit(rhs, "\\s*(\\+|\\*)\\s*")[[1]]
  rhs.indices <- gregexpr("\\+|\\*", rhs)[[1]]
  if (length(rhs.indices) == 1 && rhs.indices < 0) {
    rhs.operators <- character()
  } else {
    rhs.operators <- substring(rhs, rhs.indices, rhs.indices)
  }
  
  # preprocess dyadcov and edgecov terms
  for (k in 1:length(rhs.terms)) {
    if (grepl("((edge)|(dyad))cov", rhs.terms[k])) {
      if (grepl(",\\s*?((attr)|\\\")", rhs.terms[k])) { # with attrib argument
        x1 <- sub("((?:offset\\()?((edge)|(dyad))cov\\()([^\\)]+)((,\\s*a*.*?)\\)(?:\\))?)", "\\1", 
            rhs.terms[k], perl = TRUE)
        x2 <- sub("((?:offset\\()?((edge)|(dyad))cov\\()([^\\)]+)((,\\s*a*.*?)\\)(?:\\))?)", "\\5", 
            rhs.terms[k], perl = TRUE)
        x3 <- sub("((?:offset\\()?((edge)|(dyad))cov\\()([^\\)]+)((,\\s*a*.*?)\\)(?:\\))?)", "\\6", 
            rhs.terms[k], perl = TRUE)
      } else { # without attribute argument
        x1 <- sub("((?:offset\\()?((edge)|(dyad))cov\\()([^\\)]+)((,*\\s*a*.*?)\\)(?:\\))?)", "\\1", 
            rhs.terms[k], perl = TRUE)
        x2 <- sub("((?:offset\\()?((edge)|(dyad))cov\\()([^\\)]+)((,*\\s*a*.*?)\\)(?:\\))?)", "\\5", 
            rhs.terms[k], perl = TRUE)
        x3 <- sub("((?:offset\\()?((edge)|(dyad))cov\\()([^\\)]+)((,*\\s*a*.*?)\\)(?:\\))?)", "\\6", 
            rhs.terms[k], perl = TRUE)
      }
      type <- class(eval(parse(text = x2)))
      covnames <- c(covnames, x2)
      
      if (grepl("[^\\]]\\]$", x2)) {
        # time-varying covariate with given indices (e.g., formula[1:5])
        rhs.terms[k] <- paste(x1, x2, x3, sep = "")
        if (length(eval(parse(text = x2))) != time.steps) {
          stop(paste(x2, "has", length(eval(parse(text = x2))), 
              "elements, but there are", time.steps, "networks to be modeled."))
        }
        
      } else if (type == "matrix" || type == "network") {
        # time-independent covariate
        rhs.terms[k] <- paste(x1, x2, x3, sep = "")
      } else if (type == "list" || type == "network.list") {
        # time-varying covariate
        if (length(eval(parse(text = x2))) != time.steps) {
          stop(paste(x2, "has", length(get(x2)), "elements, but there are", 
              time.steps, "networks to be modeled."))
        }
        x2 <- paste(x2, "[[", iterator, "]]", sep = "")
        rhs.terms[k] <- paste(x1, x2, x3, sep = "")
      } else {
        stop(paste(x2, "is not a matrix, network, or list."))
      }
      
      # check if dimensions at each time step are OK
      if (!is.null(dep)) {
        for (i in 1:length(dep)) {
          cv <- eval(parse(text = x2))
          msg <- paste0("btergm error: The dimensions of covariate '", x2, 
              "' do not match the dimensions of the dependent network ", 
              "at time step ", i, ".")
          if ("list" %in% class(cv)) {
            if (any(dim(as.matrix(dep[[i]])) != dim(as.matrix(cv[[i]])))) {
              stop(msg)
            }
          } else {
            if (any(dim(as.matrix(dep[[i]])) != dim(as.matrix(cv)))) {
              stop(msg)
            }
          }
        }
      }
    }
  }
  
  assign("covnames", covnames, envir = parent.frame(n = 1))
  
  # reassemble rhs
  rhs <- rhs.terms[1]
  if (length(rhs.operators) > 0) {
    for (i in 1:length(rhs.operators)) {
      rhs <- paste(rhs, rhs.operators[i], rhs.terms[i + 1])
    }
  }
  return(rhs)
}

