
# Internal helper function, which does the actual interpretation computations. 
# Needs a list as provided by tergmprepare, including a list of 
# networks (l$networks, with possibly just one network in this list) and a 
# formula (l$form, with temporal indices, as prepared by tergmprepare). Also 
# need to supply coefficients vector, type string, i, j, and t (where t can be 
# 1 if there is only one network).
dointerpret <- function(l, coefficients, type, i, j, t) {
  for (cv in 1:length(l$covnames)) {
    assign(l$covnames[cv], l[[l$covnames[cv]]])
  }
  assign("offsmat", l$offsmat)
  form <- as.formula(l$form)
  
  # prepare i and j
  if (!is.list(i)) {
    i <- rep(list(i), length(l$networks))
    num.actors <- numeric()
    for (k in t) {  # note that t can be a vector of time steps
      num.actors[k] <- nrow(as.matrix(l$networks[[k]]))
    }
    if (length(table(num.actors)) > 1) {
      warning(paste("'i' does not vary across time steps, but the number of",
          "actors does. 'i' can be provided as a list or as a name."))
    }
  }
  if (!is.list(j)) {
    j <- rep(list(j), length(l$networks))
    num.actors <- numeric()
    for (k in t) {
      num.actors[k] <- nrow(as.matrix(l$networks[[k]]))
    }
    if (length(table(num.actors)) > 1) {
      warning(paste("'j' does not vary across time steps, but the number of",
          "actors does. 'j' can be provided as a list or as a name."))
    }
  }
  for (ll in 1:length(j)) {
    if (length(j[[ll]]) > 1 && (type == "tie" || type == "dyad")) {
      stop(paste("For computing dyadic or tie probabilities, only a single 'j'",
          "node can be specified per time step."))
    }
  }
  node_i <- i
  node_j <- j
  
  if (type == "tie") {
    results <- numeric()
    for (i in t) {
      networks[[i]][node_i[[i]], node_j[[i]]] <- 0
      stat0 <- summary(form, response = NULL)
      networks[[i]][node_i[[i]], node_j[[i]]] <- 1
      stat1 <- summary(form, response = NULL)
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
      inf <- FALSE  # next couple of lines: prevent errors with structural ones
      for (k in 1:length(coefficients)) {
        if (coefficients[k] == Inf && chgstat[k] != 0) {
          inf <- TRUE
          break
        }
      }
      if (inf == TRUE) {
        result <- 1
      } else {  # do the actual probability computations
        coefficients[coefficients == Inf] <- 0
        lp <- t(chgstat) %*% cbind(coefficients)
        result <- c(1 / (1 + exp(-lp)))
      }
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
    if (Inf %in% coefficients) {
      coefficients[coefficients == Inf] <- 9e8
      warning(paste("There are +Inf coefficients (possibly due to", 
          "offset terms). To yield interpretable results, +Inf is", 
          "approximated by a large positive number (9e8). Note that", 
          "this may be imprecise."))
    }
    results <- list()
    for (i in t) {
      # print error if undirected
      if ((is.network(l$networks[[i]]) && !is.directed(l$networks[[i]])) || 
          (is.matrix(l$networks[[i]]) && is.mat.directed(l$networks[[i]]) == FALSE)) {
        stop(paste0("Network at t=", i, " is undirected. Dyadic ", 
            "probabilities do not make sense in undirected networks. Try ", 
            "type = \"tie\" instead!"))
      }
      eta_mat <- matrix(NA, 2, 2)
      for (xi in 0:1) {
        for (xj in 0:1) {
          networks[[i]][node_i[[i]], node_j[[i]]] <- xi
          networks[[i]][node_j[[i]], node_i[[i]]] <- xj
          stat <- summary(form, response = NULL)
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
    if (Inf %in% coefficients) {
      coefficients[coefficients == Inf] <- 9e8
      warning(paste("There are +Inf coefficients (possibly due to", 
          "offset terms). To yield interpretable results, +Inf is", 
          "approximated by a large positive number (9e8). Note that", 
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
        stop(paste("Either 'i' or 'j' must contain more than one node (per",
            "time step)."))
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
        networks[[i]][ik, jk] <- vecs[l, ]
        stat <- summary(form, response = NULL)
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
  
  if (!is.null(target) && is.numeric(target)) {
    stop(paste("'target' argument (short: 't') cannot be numeric."))
  }
  
  l <- tergmprepare(formula = formula, offset = FALSE, blockdiag = FALSE, 
      verbose = FALSE)
  
  # extract response network and adjust formula
  if (!is.null(target)) {
    l$networks <- list(target)
  }
  
  # warn about indices with bipartite networks
  if (is.network(l$networks[[1]]) && is.bipartite(l$networks[[1]]) || 
      is.matrix(l$networks[[1]]) && !is.mat.onemode(l$networks[[1]])) {
    if (is.numeric(j) && j < nrow(as.matrix(l$networks[[1]])) + 1) {
      stop(paste0("In this bipartite network, the indices of the 'j' ", 
          "index start with ", nrow(as.matrix(l$networks[[1]])) + 1, 
          ". Smaller indices denote within-block entries and do not make ", 
          "any sense."))
    }
  }
  
  dointerpret(l, coefficients = coefficients, type = type, i = i, j = j, 
      t = 1)[[1]]
}


# interpretation method for btergm objects
interpret.btergm <- function(object, formula = getformula(object), 
    coefficients = coef(object), target = NULL, type = "tie", i, j, 
    t = 1:object@time.steps) {
  
  l <- tergmprepare(formula = formula, offset = FALSE, blockdiag = FALSE, 
      verbose = FALSE)
  
  # extract response networks and adjust formula
  if (!is.null(target)) {
    l$networks <- target
  }
  
  # warn about indices with bipartite networks
  for (h in t) {
    if (is.network(l$networks[[h]]) && is.bipartite(l$networks[[h]]) || 
        is.matrix(l$networks[[h]]) && !is.mat.onemode(l$networks[[h]])) {
      if (is.numeric(j) && j < nrow(as.matrix(l$networks[[h]])) + 1) {
        stop(paste0("In this bipartite network, the indices of the 'j' ", 
            "index start with ", nrow(as.matrix(l$networks[[h]])) + 1, 
            " at t=", h, ". Smaller indices denote within-block entries and ", 
            "do not make any sense."))
      }
    }
  }
  
  dointerpret(l, coefficients = coefficients, type = type, i = i, j = j, 
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
edgeprob <- function (object, verbose = FALSE) {
  if ("ergm" %in% class(object)) {
    tergm <- FALSE
  } else if ("btergm" %in% class(object) || "mtergm" %in% class(object)) {
    tergm <- TRUE
  } else {
    stop(paste("The edgeprob function is only applicable to ergm, btergm, and",
               "mtergm objects."))
  }
  l <- tergmprepare(formula = getformula(object), offset = FALSE,
                    blockdiag = FALSE, verbose = FALSE)
  for (cv in 1:length(l$covnames)) {
    assign(l$covnames[cv], l[[l$covnames[cv]]])
  }
  assign("offsmat", l$offsmat)
  form <- stats::as.formula(l$form)
  covnames <- l$covnames[-1]
  coefs <- coef(object)
  if (verbose == TRUE) {
    message("Creating data frame with predictors...")
  }
  Y <- NULL
  dyads <- NULL
  for (i in 1:length(l$networks)) {
    mat <- as.matrix(l$networks[[i]])
    imat <- matrix(rep(1:nrow(mat), ncol(mat)), nrow = nrow(mat))
    if ((is.network(l$networks[[i]]) && is.bipartite(l$networks[[i]])) ||
        (is.matrix(l$networks[[i]]) && is.mat.onemode(l$networks[[i]]) == FALSE)) {
      mn <- nrow(mat) + 1
      mx <- nrow(mat) + ncol(mat)
      jmat <- matrix(rep(mn:mx, nrow(mat)), nrow = nrow(mat),
                     byrow = TRUE)
    } else {
      jmat <- matrix(rep(1:ncol(mat), nrow(mat)), nrow = nrow(mat),
                     byrow = TRUE)
    }
    f <- stats::as.formula(paste(l$form, " + edgecov(imat) + edgecov(jmat)"))
    mpli <- ergm::ergmMPLE(f)
    Y <- c(Y, mpli$response)
    mpli$predictor <- cbind(mpli$predictor, i)
    # TODO: the previous line might say something like: "Note: Term nodeofactor("edu") skipped because it contributes no statistics." when a covariate is full of zeros.
    if (is.null(dyads)) {
      dyads <- mpli$predictor
    } else { # before rbind can be use, we need to check if some covariate levels were absent either before or now
      newColumns <- which(sapply(colnames(mpli$predictor), function(x) x %in% colnames(dyads)) == FALSE)
      if (length(newColumns) > 0) {
        for (j in 1:length(newColumns)) {
          dyads <- cbind(dyads[, 1:(newColumns[j] - 1)], 0, dyads[, newColumns[j]:ncol(dyads)])
          colnames(dyads)[newColumns[j]] <- names(newColumns)[j]
        }
      }
      notPresentAnymore <- which(sapply(colnames(dyads), function(x) x %in% colnames(mpli$predictor)) == FALSE)
      if (length(notPresentAnymore) > 0) {
        for (j in 1:length(notPresentAnymore)) {
          mpli$predictor <- cbind(mpli$predictor[, 1:(notPresentAnymore[j] - 1)], 0, mpli$predictor[, notPresentAnymore[j]:ncol(mpli$predictor)])
          colnames(mpli$predictor)[notPresentAnymore[j]] <- names(notPresentAnymore)[j]
        }
      }
      dyads <- rbind(dyads, mpli$predictor)
    }
  }
  term.names <- colnames(dyads)[-(length(colnames(dyads)):(length(colnames(dyads)) - 2))]
  term.names <- c(term.names, "i", "j", "t")
  dyads <- data.frame(dyads)
  colnames(dyads) <- term.names
  dyads <- cbind(Y, dyads)
  colnames(dyads)[1] <- "tie"
  class(dyads[, length(colnames(dyads))]) <- "integer"
  class(dyads[, length(colnames(dyads)) - 1]) <- "integer"
  class(dyads[, length(colnames(dyads)) - 2]) <- "integer"
  cf <- coef(object)
  cf.length <- length(cf)
  cf <- cf[!cf %in% c(Inf, -Inf)]
  if (length(cf) != cf.length) {
    warning(paste("There are structural zeros or ones. For these dyads, the",
                  "predicted probabilities are not valid and must be manually replaced",
                  "by 0 or 1, respectively."))
  }
  cbcoef <- cbind(cf)
  chgstat <- dyads[, 2:(ncol(dyads) - 3)]

  # handle decay term in curved ERGMs
  if ("mtergm" %in% class(object) && ergm::is.curved(object@ergm)) {
    curved.term <- vector(length = length(object@ergm$etamap$curved))
    for (i in 1:length(object@ergm$etamap$curved)) {
      curved.term[i] <- object@ergm$etamap$curved[[i]]$from[2]
    }
    cbcoef <- cbcoef[-c(curved.term)]
  } else if ("ergm" %in% class(object) && ergm::is.curved(object)) {
    curved.term <- vector(length = length(object$etamap$curved))
    for (i in 1:length(object$etamap$curved)) {
      curved.term[i] <- object$etamap$curved[[i]]$from[2]
    }
    cbcoef <- cbcoef[-c(curved.term)]
  }

  lp <- apply(chgstat, 1, function(x) t(x) %*% cbcoef)
  result <- c(1/(1 + exp(-lp)))
  i.name <- numeric(nrow(dyads))
  j.name <- numeric(nrow(dyads))
  for (t in 1:length(l$networks)) {
    vnames.t <- colnames(l$networks[[t]][, ])
    dyads.t <- dyads[which(dyads$t == t), ]
    i.name.t <- vnames.t[dyads.t$i]
    j.name.t <- vnames.t[dyads.t$j]
    i.name[which(dyads$t == t)] <- i.name.t
    j.name[which(dyads$t == t)] <- j.name.t
  }
  dyads$i.name <- i.name
  dyads$j.name <- j.name
  dyads <- cbind(dyads, result)
  colnames(dyads)[ncol(dyads)] <- "probability"
  dyads <- dyads[order(dyads$t, dyads$i, dyads$j), ]
  rownames(dyads) <- NULL
  return(dyads)
}


# function for marginal effects plots
marginalplot <- function(model, var1, var2, inter, ci = 0.95, 
    rug = FALSE, point = FALSE, structzeromat = NULL, 
    zeroline = TRUE, color = "black", xlab = NULL, ylab = NULL) {

  # check arguments
  if (!var1 %in% names(coef(model))) {
    stop("'var1' not found.")
  }
  if (!var2 %in% names(coef(model))) {
    stop("'var2' not found.")
  }

  # change statistics (needed for second variable)
  ep <- edgeprob(model)

  # remove '[[i]]' from edge covariate names in edgeprob
  for (i in 1:ncol(ep)) {
    if (grepl("((edge)|(dyad))cov", colnames(ep)[i])) {
      colnames(ep)[i] <- substr(colnames(ep)[i], 1, nchar(colnames(ep)[i]) - 5)
    }
  }
  
  # delete structural zeros from edgeprob output
  if (!is.null(structzeromat)) {
    include <- logical(nrow(ep))
    for (r in 1:length(include)) {
      if (structzeromat[ep$i[r], ep$j[r] - nrow(structzeromat)] == 1) {
        include[r] <- FALSE
      } else {
        include[r] <- TRUE
      }
    }
    ep <- ep[include == TRUE, ]
  }
  
  # unique values of the second variable
  v2 <- sort(unlist(unique(ep[var2])))
  names(v2) <- NULL
  
  # coefficients
  co <- coef(model)
  beta1 <- co[match(var1, names(co))]
  beta3 <- co[match(inter, names(co))]
  
  # marginal effects
  delta1 = beta1 + beta3 * v2
  
  # variance-covariance matrix
  if ("mtergm" %in% class(model)) {
    model <- model@ergm
  }
  vc <- stats::vcov(model)
  
  # variances
  variances1 = vc[var1, var1] + (v2^2) * vc[inter, inter] + 
      (2 * v2 * vc[var1, inter])
  
  # standard errors
  se1 = sqrt(variances1)
  
  # confidence intervals
  z_score = qnorm(1 - ((1 - ci) / 2))
  upper = delta1 + z_score * se1
  lower = delta1 - z_score * se1
  
  # aggregate in data frame
  dta <- data.frame(v2 = v2, delta1 = delta1, upper = upper, lower = lower)
  
  # axis labels
  if (is.null(xlab)) {
    xlab <- var2
  }
  if (is.null(ylab)) {
    ylab <- var1
  }
  
  # plot error bars or line
  if (point == TRUE) {
    if (length(v2) == 2 && v2[1] == 0 && v2[2] == 1 && rug == FALSE) {
      dta$v2 <- as.factor(dta$v2)
    }
    gp <- ggplot(data = dta, aes(x = v2, y = delta1)) + 
      geom_errorbar(data = dta, 
                    aes(x = v2, ymin = lower, ymax = upper), 
                    color = color) + 
      geom_point(data = dta, aes(x = v2, y = delta1)) + 
      ylab(ylab) + 
      xlab(xlab)
  } else {
    gp <- ggplot(data = dta, aes(x = v2, y = delta1)) + 
      geom_line(color = color) + 
      geom_ribbon(aes_string(ymin = "lower", ymax = "upper"), 
                  alpha = 0.15, 
                  fill = color) + 
      ylab(ylab) + 
      xlab(xlab)
  }
  
  # add distribution to the plot
  if (rug == TRUE) {
    gp <- gp + geom_rug(data = ep[var2], 
                        aes_string(x = var2), 
                        inherit.aes = FALSE, 
                        sides = "b", 
                        col = color, 
                        alpha = 0.1)
  }
  
  # add horizontal line for 0
  if (zeroline == TRUE) {
    gp <- gp + geom_hline(yintercept = 0, linetype = "dashed")
  }
  
  return(gp)
}
