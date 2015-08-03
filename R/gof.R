
# function which reduces a statistic x nsim matrix and computes summary stats
# input: two matrices of a certain type of statistics (simulated and observed)
# goal: get rid of empty rows at the end, e.g., where dsp(37) or so is usually
# not observed; return value: object containing the summary statistics
reduce.matrix <- function(sim, obs) {
  
  numsim <- ncol(as.matrix(sim))
  numobs <- ncol(as.matrix(obs))
  
  # if geodist statistic: put aside the last 'Inf' row
  if (is.null(rownames(sim)) || rownames(sim)[nrow(sim)] == "Inf") {
    geo <- TRUE
    inf.sim <- sim[nrow(sim), ]  # put aside this row for now and reuse later
    sim <- sim[-nrow(sim), ]
    if (class(obs) == "matrix") {
      inf.obs <- obs[nrow(obs), ]
      obs <- matrix(obs[-nrow(obs), ], ncol = numobs)
    } else {
      inf.obs <- obs[length(obs)]
      obs <- matrix(obs[-length(obs)], ncol = numobs)
    }
  } else {
    geo <- FALSE
  }
  
  # find first empty row for simulation matrix
  sim.rs <- rowSums(sim)
  sim.remove <- length(sim.rs)  # at which row index can removal start?
  for (i in (length(sim.rs) - 1):1) {
    if (sim.rs[i] == 0 && sim.remove == (i + 1)) {
      sim.remove <- i  # remember which is the first empty row
    }
  }
  #sim.remove <- sim.remove + 1  # keep one row without observations
  if (class(obs) != "matrix") {  # one network is compared
    obs.remove <- length(obs)
    for (i in (length(obs) - 1):1) {
      if (obs[i] == 0 && obs.remove == (i + 1)) {
        obs.remove <- i
      }
    }
  } else {  # several networks are compared
    obs.rs <- rowSums(obs)
    obs.remove <- length(obs.rs)
    for (i in (length(obs.rs) - 1):1) {
      if (obs.rs[i] == 0 && obs.remove == (i + 1)) {
        obs.remove <- i
      }
    }
  }
  rem <- max(c(obs.remove, sim.remove), na.rm = TRUE)  # which one is longer?
  
  # remove unnecessary rows
  if (class(obs) != "matrix") {  # get rid of empty observations or rows of obs
    obs <- matrix(obs[-(rem:length(obs))], ncol = numobs)
  } else {
    obs <- matrix(obs[-(rem:nrow(obs)), ], ncol = numobs)
  }
  sim <- matrix(sim[-(rem:nrow(sim)), ], ncol = numsim)  # same for sim stats
  
  if (nrow(obs) < rem) {
    for (i in (nrow(obs) + 1):rem) {
      obs <- rbind(obs, rep(0, ncol(obs)))
    }
  }
  if (nrow(sim) < rem) {
    for (i in (nrow(sim) + 1):rem) {
      sim <- rbind(sim, rep(0, ncol(sim)))
    }
  }
  
  # for geodist, add Inf row again
  if (geo == TRUE) {
    sim <- rbind(sim, inf.sim)
    obs <- rbind(obs, inf.obs)
    rownames(sim) <- c(1:(nrow(sim) - 1), "Inf")
  }
  
  # create final object which will contain the raw simulations + the comparison
  reducedobject <- list()
  reducedobject$sim <- as.data.frame(sim)
  
  rownames(sim) <- NULL
  rownames(obs) <- NULL
  
  # compute means, p values, etc. and put them in a data frame
  x <- matrix()
  if (ncol(obs) == 1 || ncol(sim) == 1) {  # compute all the summary statistics
    x.obs <- obs
    x.mean <- apply(sim, 1, mean)
    x.min <- apply(sim, 1, min)
    x.max <- apply(sim, 1, max)
    x.median <- apply(sim, 1, median)
    zval <- (x.mean - x.obs) / sd(x.mean)
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
    x.pval <- pval
    x <- data.frame(x.obs, x.mean, x.median, x.min, x.max, x.pval)
    colnames(x) <- c("obs", "sim: mean", "median", "min", "max", "Pr(>z)")
  } else {  # for several target networks, compute distribution
    x.obs.mean <- apply(obs, 1, mean)
    x.obs.min <- apply(obs, 1, min)
    x.obs.max <- apply(obs, 1, max)
    x.obs.median <- apply(obs, 1, median)
    x.mean <- apply(sim, 1, mean)
    x.min <- apply(sim, 1, min)
    x.max <- apply(sim, 1, max)
    x.median <- apply(sim, 1, median)
    x.pval <- numeric()
    for (i in 1:nrow(sim)) {
      tryCatch(
        expr = {
          x.pval[i] <- t.test(obs[i, ], sim[i, ])$p.value  # compare group means
        }, 
        error = function(e) {
          x.pval[i] <- 1  # if both are 0, a t.test cannot be computed...
        }, 
        finally = {}
      )
      
      if (is.nan(x.pval[i])) {  # geodist contains "Inf"
        x.pval[i] <- 1
      }
    }
    x <- data.frame(x.obs.mean, x.obs.median, x.obs.min, x.obs.max, x.mean, 
        x.median, x.min, x.max, x.pval)
    colnames(x) <- c("obs: mean", "median", "min", "max", 
        "sim: mean", "median", "min", "max", "Pr(>z)")
  }
  if (geo == TRUE) {
    rownames(x) <- c(1:(nrow(x) - 1), "Inf")
  } else {
    rownames(x) <- 0:(nrow(x) - 1)
  }
  reducedobject$comparison <- as.data.frame(x)

  return(reducedobject)
}


# classic statnet-like GOF function
statnetgof <- function(gofobject, simulations, target, statistics = c("dsp", 
    "esp", "geodist", "degree", "idegree", "odegree", "kstar", "istar", 
    "ostar"), ...) {
  
  message("\nComputing classic (statnet-style) goodness of fit:")
  directed <- simulations[[1]]$gal$directed
  bipartite <- network::is.bipartite(network::network(target[[1]]))
  
  # create lists for simulated and observed network stats
  stats <- list()
  raw <- list()
  
  num.stats <- length(statistics)
  mx <- max(gofobject$num.vertices)  # maximum number of nodes
  for (i in 1:num.stats) {
    stat.i <- statistics[i]
    sim <- NULL
    obs <- NULL
    if (stat.i == "dsp") {  # dsp
      message(paste0("[", i, "/", num.stats, 
          "] Comparison for dyad-wise shared partners..."))
      sim <- sapply(simulations, function(x) summary(x ~ dsp(0:mx)))
      obs <- sapply(target, function(x) summary(x ~ dsp(0:mx)))
    } else if (stat.i == "esp") {  # esp
      message(paste0("[", i, "/", num.stats, 
          "] Comparison for edge-wise shared partners..."))
      if (bipartite == FALSE) {
        sim <- sapply(simulations, function(x) summary(x ~ esp(0:mx)))
        obs <- sapply(target, function(x) summary(x ~ esp(0:mx)))
      } else {
        message(" -> Bipartite networks. Ignoring the 'esp' statistic.")
      }
    } else if (stat.i == "geodist") {  # geodist
      message(paste0("[", i, "/", num.stats, 
          "] Comparison for geodesic distances..."))
      fillup <- function(x, another.length) {  # fill up x if shorter
        difference <- length(x) - another.length
        inf.value <- x[length(x)]
        if (difference < 0) {  # x is shorter
          x <- x[1:(length(x) - 1)]
          x <- c(x, rep(0, abs(difference)), inf.value)
        } else if (difference > 0) {
          x <- x[1:(length(x) - difference)]
          x <- c(x, inf.value)
        }
        return(x)
      }
      sim <- sapply(simulations, function(x) fillup(ergm.geodistdist(x), mx))
      obs <- sapply(target, function(x) fillup(ergm.geodistdist(x), mx))
    } else if (stat.i == "degree") {  # degree
      message(paste0("[", i, "/", num.stats, "] Comparison for degree..."))
      if (directed == TRUE) {
        message(" -> Directed networks. Ignoring the 'degree' statistic.")
      } else {
        sim <- sapply(simulations, function(x) summary(x ~ degree(0:mx)))
        obs <- sapply(target, function(x) summary(x ~ degree(0:mx)))
      }
    } else if (stat.i == "idegree") {  # idegree
      message(paste0("[", i, "/", num.stats, "] Comparison for indegree..."))
      if (directed == FALSE) {
        message(" -> Undirected networks. Ignoring the 'idegree' statistic.")
      } else {
        sim <- sapply(simulations, function(x) summary(x ~ idegree(0:mx)))
        obs <- sapply(target, function(x) summary(x ~ idegree(0:mx)))
      }
    } else if (stat.i == "odegree") {  # odegree
      message(paste0("[", i, "/", num.stats, "] Comparison for outdegree..."))
      if (directed == FALSE) {
        message(" -> Undirected networks. Ignoring the 'odegree' statistic.")
      } else {
        sim <- sapply(simulations, function(x) summary(x ~ odegree(0:mx)))
        obs <- sapply(target, function(x) summary(x ~ odegree(0:mx)))
      }
    } else if (stat.i == "kstar") {  # kstar
      message(paste0("[", i, "/", num.stats, "] Comparison for k-stars..."))
      if (directed == TRUE) {
        message(" -> Directed networks. Ignoring the 'kstar' statistic.")
      } else {
        sim <- sapply(simulations, function(x) summary(x ~ kstar(1:mx)))
        obs <- sapply(target, function(x) summary(x ~ kstar(1:mx)))
      }
    } else if (stat.i == "istar") {  # istar
      message(paste0("[", i, "/", num.stats, "] Comparison for in-stars..."))
      if (directed == FALSE) {
        message(" -> Undirected networks. Ignoring the 'istar' statistic.")
      } else {
        sim <- sapply(simulations, function(x) summary(x ~ istar(0:mx)))
        obs <- sapply(target, function(x) summary(x ~ istar(0:mx)))
      }
    } else if (stat.i == "ostar") {  # ostar
      message(paste0("[", i, "/", num.stats, "] Comparison for out-stars..."))
      if (directed == FALSE) {
        message(" -> Undirected networks. Ignoring the 'ostar' statistic.")
      } else {
        sim <- sapply(simulations, function(x) summary(x ~ ostar(0:mx)))
        obs <- sapply(target, function(x) summary(x ~ ostar(0:mx)))
      }
    } else {
      message(paste0("[", i, "/", num.stats, "] ", stat.i, 
          " not recognized. Skipping this statistic."))
    }
    if (!is.null(sim) && !is.null(obs)) {
      reduced <- reduce.matrix(sim, obs)
      stats <- c(stats, list(reduced$comparison))
      raw <- c(raw, list(Matrix(as.matrix(reduced$sim))))
      names(stats)[length(stats)] <- stat.i
      names(raw)[length(raw)] <- stat.i
    }
  }
  
  gofobject$raw <- raw
  gofobject$statistics <- stats
  return(gofobject)
}


# AUC-PR function -- modified version of: https://github.com/ipa-tys/ROCR/pull/2
aucpr <- function(pred, precision, recall) {
  falsepos <- pred@fp
  truepos <- pred@tp
  n.positive <- pred@n.pos
  aucvalues <- numeric()
  for (j in 1:length(precision)) {
    fp <- falsepos[[j]]
    tp <- truepos[[j]]
    n.pos <- n.positive[[j]]
    prec <- precision[[j]]
    rec <- recall[[j]]
    
    # if two points are too distant from each other, we need to
    # correctly interpolate between them. This is done according to
    # Davis & Goadrich,
    #"The Relationship Between Precision-Recall and ROC Curves", ICML'06
    for (i in seq_along(rec[-length(rec)])) {
      if (tp[i + 1] - tp[i] > 2) {
        skew = (fp[i + 1] - fp[i]) / (tp[i + 1] - tp[i])
        x = seq(1, tp[i + 1] - tp[i], by = 1)
        rec <- append(rec, (x + tp[i]) / n.pos, after = i)
        prec <- append(prec, (x + tp[i]) / (tp[i] + fp[i] + x + skew * x), 
            after = i)
      }
    }
    
    auc <- 0
    for (i in 2:length(rec)) {
        auc <- auc + 0.5 * (rec[i] - rec[i-1]) * (prec[i] + prec[i-1])
    }
    aucvalues <- c(aucvalues, auc)
  }
  return(aucvalues)
}


# compute ROC and PR curves and the area under the curve (AUC)
rocprgof <- function(gofobject, simulations, target, missingmatrix = NULL, 
    na.method = "remove", rgraph = FALSE, pr.impute = "poly4") {
  
  if (rgraph == FALSE) {
    message(paste("\nComputing classification performance (ROC/PR/AUC)", 
        "goodness of fit..."))
  }
  
  # remove missing data (absent obs in the target nw) from the simulations
  if (!is.null(missingmatrix)) {
    for (i in 1:length(simulations)) {
      sim <- as.matrix(simulations[[i]])
      if (nrow(sim) == nrow(missingmatrix) && ncol(sim) == 
          ncol(missingmatrix)) {
        sim[missingmatrix == TRUE] <- NA
        sim <- suppressMessages(handleMissings(sim, na = NA, 
            method = na.method))
        simulations[[i]] <- network::network(sim)
      } else {
        stop(paste("Missing data matrix and simulations do not have the", 
            "same number of nodes."))
      }
    }
  }
  
  # check if ROC, PR, and AUC can be computed
  compute.rocr <- TRUE
  nrow.sim <- sapply(simulations, function(x) nrow(as.matrix(x)))
  if (length(table(nrow.sim)) > 1) {
    compute.rocr <- FALSE
    if (rgraph == FALSE) {
      message(paste("ROC, PR, and AUC cannot be computed because the number", 
          "of rows differs across simulations."))
    }
  } else {
    nrow.sim <- nrow.sim[1]
  }
  ncol.sim <- sapply(simulations, function(x) ncol(as.matrix(x)))
  if (length(table(ncol.sim)) > 1) {
    compute.rocr <- FALSE
    if (rgraph == FALSE) {
      message(paste("ROC, PR, and AUC cannot be computed because the number", 
          "of columns differs across simulations."))
    }
  } else {
    ncol.sim <- ncol.sim[1]
  }
  nrow.tar <- sapply(target, function(x) nrow(as.matrix(x)))
  if (length(table(nrow.tar)) > 1) {
    compute.rocr <- FALSE
    if (rgraph == FALSE) {
      message(paste("ROC, PR, and AUC cannot be computed because the number", 
          "of rows differs across target networks."))
    }
  } else {
    nrow.tar <- nrow.tar[1]
  }
  ncol.tar <- sapply(target, function(x) ncol(as.matrix(x)))
  if (length(table(ncol.tar)) > 1) {
    compute.rocr <- FALSE
    if (rgraph == FALSE) {
      message(paste("ROC, PR, and AUC cannot be computed because the number", 
          "of columns differs across target networks."))
    }
  } else {
    ncol.tar <- ncol.tar[1]
  }
  if (compute.rocr == TRUE) {
    if ((nrow.sim != nrow.tar) || (ncol.sim != ncol.tar)) {
      compute.rocr <- FALSE
      if (rgraph == FALSE) {
        message(paste("ROC, PR, and AUC cannot be computed because the", 
            "number of nodes differs between simulations and target", 
            "networks."))
      }
    }
  }
  
  # compute ROC, PR, and AUC
  if (compute.rocr == TRUE) {
    target.pr <- list()
    target.y <- list()
    for (j in 1:length(target)) {
      net <- target[[j]]
      sums <- 0
      for (i in 1:length(simulations)) {
        sums <- sums + as.matrix(simulations[[i]])
      }
      sums <- sums / length(simulations)
      if (network::is.directed(target[[j]]) == TRUE || 
          network::is.bipartite(target[[j]]) == TRUE) {
        pr <- c(sums)
        y <- c(as.matrix(net))
      } else {
        pr <- sums[lower.tri(sums)]
        y <- as.matrix(net)[lower.tri(as.matrix(net))]
      }
      pr <- pr[!is.na(y)]
      y <- y[!is.na(y)]
      target.pr[[j]] <- pr
      target.y[[j]] <- y
    }
    pred <- prediction(target.pr, target.y)
    roc <- performance(pred, "tpr", "fpr")  # ROC curve
    pr <- performance(pred, "ppv", "tpr")  # precision-recall curve
    
    # impute the first PR value (usually NaN)
    for (j in 1:length(pr@y.values)) {
      fp <- pred@fp[[j]]
      tp <- pred@tp[[j]]
      if (fp[1] == 0 & tp[1] == 0) {
        if (pr.impute != "no" && rgraph == FALSE) {
          message(paste0("t = ", j, ": imputation of the first PR value not ", 
              "necessary."))
        }
        pr@y.values[[j]][1] <- 1
      } else if (pr.impute == "no" && rgraph == FALSE) {
        message(paste0("t = ", j, ": warning -- the first PR value was not ", 
          "imputed; this may lead to underestimated AUC-PR values."))
        # do nothing
      } else if (is.nan(pr@y.values[[j]][1])) {
        if (pr.impute == "second") {
          if (rgraph == FALSE) {
            message(paste0("t = ", j, ": imputing the first PR value by the ", 
                "next (=adjacent) value."))
          }
          pr@y.values[[j]][1] <- pr@y.values[[j]][2]
        } else if (pr.impute == "one") {
          if (rgraph == FALSE) {
            message(paste0("t = ", j, ": imputing the first PR value by the ", 
                "maximum value of 1."))
          }
          pr@y.values[[j]][1] <- 1
        } else if (grepl("^poly[1-9]", pr.impute)) {
          num <- as.numeric(substr(pr.impute, 5, 5))
          if (rgraph == FALSE) {
            message(paste0("t = ", j, ": imputing the first PR value using a ", 
                "polynomial of order ", num, ". Check the results by plotting ",
                "the GOF object using the \"pr.poly = ", num, "\" argument."))
          }
          p <- data.frame(poly(pr@x.values[[j]], num, raw = TRUE))
          fit <- lm(pr@y.values[[j]] ~ ., data = p)
          pr@y.values[[j]][1] <- predict(fit, newdata = p[1, ])
        } else if (rgraph == FALSE) {
          message(paste0("t = ", j, ": PR imputation method not recognized. ", 
              "Not using any imputation."))
        }
        if (pr@y.values[[j]][1] < 0) {
          pr@y.values[[j]][0] <- 0
        }
        if (pr@y.values[[j]][1] > 1) {
          pr@y.values[[j]][1] <- 1
        }
      }
    }
    
    auc.roc <- unlist(performance(pred, measure = "auc")@y.values)  # ROC-AUC
    auc.pr <- aucpr(pred, precision = pr@y.values, recall = pr@x.values) #PR-AUC
    
    if (rgraph == TRUE) {
      gofobject$rgraph.roc <- roc
      gofobject$rgraph.pr <- pr
      gofobject$rgraph.auc.roc <- auc.roc
      gofobject$rgraph.auc.pr <- auc.pr
      message("Done.")
    } else {
      gofobject$roc <- roc
      gofobject$pr <- pr
      gofobject$auc.roc <- auc.roc
      gofobject$auc.pr <- auc.pr
    }
  } else {
    gofobject$pr <- NULL
    gofobject$roc <- NULL
    gofobject$auc.roc <- NULL
    gofobject$auc.pr <- NULL
  }
  return(gofobject)
}

setMethod("getformula", signature = className("btergm", "btergm"), 
    definition = function(x) x@formula)

setMethod("getformula", signature = className("ergm", "ergm"), 
    definition = function(x) x$formula)


# create random graphs with the corresponding tie probability of each time step
randomgraph <- function(networks, nsim = 100) {
  di <- network::is.directed(networks[[1]])
  bip <- network::is.bipartite(networks[[1]])
  rownum <- sapply(networks, function(x) nrow(as.matrix(x)))
  colnum <- sapply(networks, function(x) ncol(as.matrix(x)))
  dens <- sapply(networks, sna::gden)
  rg <- list()
  for (i in 1:length(networks)) {
    rlist <- list()
    for (j in 1:nsim) {
      r <- rbinom(n = rownum[i] * colnum[i], size = 1, prob = dens[i])
      r <- matrix(r, nrow = rownum[i], ncol = colnum[i])
      if (bip == FALSE) {
        diag(r) <- 0
        if (di == FALSE) {
          r <- symmetrize(r, rule = "upper")
        }
      }
      rlist[[j]] <- r
    }
    rg <- c(rg, rlist)
  }
  return(rg)
}


# GOF function for in-sample and out-of-sample GOF assessment
gof.btergm <- function(object, target = NULL, formula = getformula(object), 
    nsim = 100, MCMC.interval = 1000, MCMC.burnin = 10000, 
    parallel = c("no", "MPI", "SOCK"), ncpus = 1, cl = NULL, 
    classicgof = TRUE, rocprgof = TRUE, checkdegeneracy = TRUE, 
    statistics = c("dsp", "esp", "geodist", "degree", "idegree", "odegree", 
    "kstar", "istar", "ostar"), pr.impute = "poly4", verbose = TRUE, ...) {
  
  if (nsim < 2) {
    stop("The 'nsim' argument must be greater than 1.")
  }
  
  simulations <- list()
  
  if (checkdegeneracy == TRUE) {
    target.stats <- list()
    degensim <- list()
  }
  
  # call tergmprepare and integrate results as a child environment in the chain
  if (class(object) == "btergm") {
    env <- tergmprepare(formula = formula, offset = object@offset, 
        verbose = verbose)
    parent.env(env) <- environment()
    offset <- object@offset
  } else {
    env <- tergmprepare(formula = formula, offset = FALSE, verbose = FALSE)
    parent.env(env) <- environment()
    offset <- FALSE
  }
  
  # check and rearrange target network(s)
  if (is.null(target)) {
    message(paste("\nNo 'target' network(s) provided. Using networks on the",
        "left-hand side of the model formula as observed networks.\n"))
    target <- env$networks
  } else if (class(target) == "network" || class(target) == "matrix") {
    target <- list(target)
    message("\nOne observed ('target') network was provided.\n")
  } else if (class(target) == "list") {
    message(paste("\n", length(target), "observed ('target') networks were",
        "provided.\n"))
  } else {
    stop("'target' must be a network, matrix, or list of matrices or networks.")
  }
  
  # extract coefficients from object
  if (class(object) == "btergm" && offset == TRUE) {
    coefs <- c(coef(object), -Inf)  # -Inf for offset matrix
  } else {
    coefs <- coef(object)
  }
  
  # adjust formula at each step, and simulate networks
  sim <- list()
  tstats <- list()
  degen <- list()
  for (index in 1:env$time.steps) {
    # simulations for statnet-style and rocpr GOF
    if (classicgof == TRUE || rocprgof == TRUE) {
      if (verbose == TRUE) {
        if (class(object) == "ergm") {
          f.i <- paste(deparse(formula), collapse = "")
          f.i <- gsub("\\s+", " ", f.i)
        } else {
          f.i <- gsub("\\[\\[i\\]\\]", paste0("[[", index, "]]"), 
              paste(deparse(env$form), collapse = ""))
          f.i <- gsub("\\s+", " ", f.i)
          f.i <- gsub("^networks", env$lhs.original, f.i)
        }
        message(paste("Simulating", nsim, 
            "networks from the following formula:\n", f.i, "\n"))
      }
      if (parallel[1] %in% c("no", "SOCK", "MPI")) {
        if (parallel[1] == "no") {
          ncpus <- 0
          p <- NULL
        }
        i <- index
        sim[[index]] <- simulate.formula(env$form, nsim = nsim, 
            coef = coefs, constraints = ~ ., 
            control = control.simulate.formula(MCMC.interval = 
            MCMC.interval, MCMC.burnin = MCMC.burnin, parallel = ncpus, 
            parallel.type = p))
      } else {
        stop(paste0("For this type of object, \"SOCK\" and \"MPI\" parallel ", 
            "processing are allowed. \"", parallel, 
            "\" is not a valid value for the \"parallel\" argument."))
      }
    }
    
    # compute target stats if degeneracy check is switched on
    if (checkdegeneracy == TRUE) {
      tstats[[index]] <- summary(ergm::remove.offset.formula(env$form), 
          response = NULL)
      degen[[index]] <- simulate.formula(env$form, nsim = nsim, 
          coef = coefs, statsonly = TRUE, 
          control = control.simulate.formula(MCMC.interval = 
          MCMC.interval, MCMC.burnin = MCMC.burnin))
    }
  }
  
  # check basis network(s)
  if (env$time.steps == 1) {
    message("One network from which simulations are drawn was provided.")
  } else {
    message(paste(object@time.steps, "networks from which simulations are",
        "drawn were provided."))
  }
  
  # unpack nested lists of simulations and target statistics
  simulations <- list()
  if (classicgof == TRUE || rocprgof == TRUE) {
    for (i in 1:length(sim)) {
      for (j in 1:length(sim[[i]])) {
        simulations[[length(simulations) + 1]] <- sim[[i]][[j]]
      }
    }
    rm(sim)
  }
  if (checkdegeneracy == TRUE) {
    degensim <- matrix(nrow = 0, ncol = ncol(degen[[1]]))
    target.stats <- list()
    for (i in 1:length(degen)) {
      degensim <- rbind(degensim, degen[[i]])
      target.stats[[i]] <- tstats[[i]]
    }
    if (offset == TRUE) {
      degensim <- degensim[, -ncol(degensim)]  # get rid of offset statistic
    }
    rm(tstats)
    rm(degen)
  }
  
  # if NA in target networks, put them in the base network, too, and vice-versa
  if (length(env$networks) == length(target)) {
    for (i in 1:env$time.steps) {
      env$networks[[i]] <- as.matrix(env$networks[[i]])
      env$networks[[i]][is.na(as.matrix(target[[i]]))] <- NA
      env$networks[[i]] <- network::network(env$networks[[i]], 
          directed = env$directed, bipartite = env$bipartite)
      target[[i]] <- as.matrix(target[[i]])
      target[[i]][is.na(as.matrix(env$networks[[i]]))] <- NA
      target[[i]] <- network::network(target[[i]], directed = env$directed, 
          bipartite = env$bipartite)
    }
  }
  
  # create an object where the final results are stored
  gofobject <- list()
  class(gofobject) <- "btergmgof"
  gofobject$numbasis <- length(env$networks)
  gofobject$numtarget <- length(target)
  num.vertices <- max(sapply(env$networks, function(object) 
      network::get.network.attribute(network::network(object, 
      directed = env$directed, bipartite = env$bipartite), "n")))
  gofobject$num.vertices <- num.vertices
  
  # comparison of simulated and observed network statistics (statnet-style GOF)
  if (classicgof == TRUE) {
    gofobject <- statnetgof(gofobject, simulations, target, 
        statistics = statistics)
  } else {
    gofobject$statistics <- NULL
    gofobject$raw <- NULL
  }
  
  # ROC and PR curves and AUC measures
  if (rocprgof == TRUE) {
    gofobject <- rocprgof(gofobject, simulations, target, rgraph = FALSE, 
        pr.impute = pr.impute)
    randomgraphs <- randomgraph(env$networks, nsim = nsim)
    gofobject <- rocprgof(gofobject, randomgraphs, target, rgraph = TRUE, 
        pr.impute = pr.impute)
  } else {
    gofobject$roc <- NULL
    gofobject$pr <- NULL
    gofobject$auc.roc <- NULL
    gofobject$auc.pr <- NULL
    gofobject$rgraph.roc <- NULL
    gofobject$rgraph.pr <- NULL
    gofobject$rgraph.auc.roc <- NULL
    gofobject$rgraph.auc.pr <- NULL
  }
  
  # degeneracy check
  if (checkdegeneracy == TRUE) {
    message("\nChecking degeneracy...")
    mat <- list()
    for (i in 1:env$time.steps) {
      sm <- coda::as.mcmc.list(coda::as.mcmc(degensim))
      sm <- ergm::sweep.mcmc.list(sm, target.stats[[i]], "-")
      center <- TRUE
      ds <- ergm::colMeans.mcmc.list(sm) - if (!center) target.stats[[i]] else 0
      sds <- apply(degensim, 2, sd)
      ns <- coda::effectiveSize(sm)
      se <- sds * sqrt(ns)
      z <- ds / se
      p.z <- pnorm(abs(z), lower.tail = FALSE) * 2
      mat[[i]] <- cbind("obs" = target.stats[[i]], "sim" = colMeans(degensim), 
          "est" = ds, "se" = se, "zval" = z, "pval" = p.z)
    }
    gofobject$degeneracy <- mat
    message("Done.")
  } else {
    gofobject$degeneracy <- NULL
  }
  
  if ((length(gofobject$auc.roc == 1) && gofobject$auc.roc >= 0.98) || 
      (length(gofobject$auc.roc > 1) && any(gofobject$auc.roc >= 0.98))) {
    warning(paste("Extremely high AUC-ROC value(s)! This may be a sign of", 
        "degenerate MCMC simulations. You should try to tweak the 'nsim',", 
        "'MCMC.interval', and 'MCMC.burnin' arguments."))
  }
  
  if ((length(gofobject$auc.pr == 1) && gofobject$auc.pr >= 0.98) || 
      (length(gofobject$auc.pr > 1) && any(gofobject$auc.pr >= 0.98))) {
    warning(paste("Extremely high AUC-PR value(s)! This may be a sign of", 
        "degenerate MCMC simulations. You should try to tweak the 'nsim',", 
        "'MCMC.interval', and 'MCMC.burnin' arguments."))
  }
  
  return(gofobject)
}


# GOF function for SIENA (creates btergm-compatible GOF objects)
gof.sienaAlgorithm <- function(object, siena.data, siena.effects, 
    predict.period = NULL, nsim = 50, parallel = c("no", "multicore", 
    "snow"), ncpus = 1, cl = NULL, target.na = NA, 
    target.na.method = "remove", target.structzero = 10, 
    classicgof = TRUE, rocprgof = TRUE, statistics = c("dsp", "esp", 
    "geodist", "degree", "idegree", "odegree", "kstar", "istar", 
    "ostar"), pr.impute = "poly4", ...) {
  
  # check RSiena version
  if (!requireNamespace("RSiena", quietly = TRUE)) {
    stop("Please install the RSiena package to use this method.")
  }
  if (packageVersion("RSiena") < as.package_version("1.0.12.169")) {
    stop("RSiena (>= 1.0.12.169) is required.")
  }
  
  # check and prepare arguments for SIENA
  if ((!"sienaModel" %in% class(object) && packageVersion("RSiena") < 
      as.package_version("1.1-227")) || (!"sienaAlgorithm" %in% 
      class(object) && packageVersion("RSiena") >= 
      as.package_version("1.1-227"))) {
    if (packageVersion("RSiena") < as.package_version("1.1-227")) {
      stop(paste("'object' must be an object of class 'sienaModel'.", 
          "Please use the sienaModelCreate() function to create such an", 
          "object."))
    } else {
      stop(paste("'object' must be an object of class 'sienaAlgorithm'.", 
          "Please use the sienaAlgorithmCreate() function to create such an", 
          "object."))
    }
  }
  if (!"siena" %in% class(siena.data)) {
    stop(paste("'siena.data' must be an object of class 'siena'.", 
        "Please use the sienaDataCreate() function to create such an object."))
  }
  if (!"sienaEffects" %in% class(siena.effects)) {
    stop(paste("'siena.effects' must be an object of class 'sienaEffects'.", 
        "Please use the getEffects() and includeEffects() functions to create", 
        "such an object."))
  }
  if (nsim < 2) {
    stop("The 'nsim' argument must be greater than 1.")
  }
  if (is.null(predict.period)) {
    base <- siena.data$observations - 1
  } else {
    base <- predict.period - 1
  }
  message(paste0("The network at time step ", base + 1, 
      " is predicted based on the last simulation at time step ", base, "."))
  message(paste("Simulating", nsim, "networks. This may take a long time."))
  
  # When an old RSiena version is installed, some internal helper functions 
  # are not available yet. In this case, these functions are embedded here. 
  # They were copied from sienaGOF.r in RSiena revision r267:
  if (packageVersion("RSiena") < as.package_version("1.1-231")) {
    
    changeToStructural <- function(X, S) {
      if (any(S >= 10, na.rm = TRUE)) {
        S[is.na(S)] <- 0
        S0 <- Matrix(S == 10)
        S1 <- Matrix(S == 11)
        X <- 1 * ((X - S0 + S1) >= 1)
      }
      X[is.na(X)] <- 0
      drop0(X)
    }
    
    changeToNewStructural <- function(X, SBefore, SAfter) {
      SB <- Matrix(SBefore >= 10)
      SA <- Matrix(SAfter >= 10)
      if (any(SA > SB, na.rm = TRUE)) {
        S0 <- (SA > SB) * Matrix(SAfter == 10)
        S1 <- (SA > SB) * Matrix(SAfter == 11)
        X <- 1 * ((X - S0 + S1) >= 1)
      }
      X[is.na(X)] <- 0
      drop0(X)
    }
  
    sparseMatrixExtraction <- function(i, obsData, sims, period, groupName, 
        varName) {
      dimsOfDepVar<- attr(obsData[[groupName]]$depvars[[varName]], "netdims")
      if (attr(obsData[[groupName]]$depvars[[varName]], "sparse")) {
        missings <- (is.na(obsData[[groupName]]$depvars[[varName]][[period]]) | 
            is.na(obsData[[groupName]]$depvars[[varName]][[period + 1]])) * 1
      } else {
        missings <- Matrix(
            (is.na(obsData[[groupName]]$depvars[[varName]][, , period]) |
            is.na(obsData[[groupName]]$depvars[[varName]][, , period + 1])) * 1)
      }
      if (is.null(i)) {
        if (attr(obsData[[groupName]]$depvars[[varName]], "sparse")) {
          returnValue <- drop0(Matrix(
              obsData[[groupName]]$depvars[[varName]][[period + 1]] %% 10))
          returnValue[is.na(returnValue)] <- 0
          returnValue <- changeToStructural(returnValue, 
              Matrix(obsData[[groupName]]$depvars[[varName]][[period]]))
        } else {
          returnValue <- Matrix(
              obsData[[groupName]]$depvars[[varName]][, , period + 1] %% 10)
          returnValue[is.na(returnValue)] <- 0
          returnValue <- changeToStructural(returnValue, 
              Matrix(obsData[[groupName]]$depvars[[varName]][, , period]))
        }
      } else {
        returnValue <- sparseMatrix(
            sims[[i]][[groupName]][[varName]][[period]][, 1],
            sims[[i]][[groupName]][[varName]][[period]][, 2],
            x = sims[[i]][[groupName]][[varName]][[period]][, 3],
            dims = dimsOfDepVar[1:2]
        )
        if (attr(obsData[[groupName]]$depvars[[varName]], "sparse")) {
          returnValue <- changeToNewStructural(returnValue,
              Matrix(obsData[[groupName]]$depvars[[varName]][[period]]),
              Matrix(obsData[[groupName]]$depvars[[varName]][[period + 1]]))
        } else {
          returnValue <- changeToNewStructural(returnValue,
              Matrix(obsData[[groupName]]$depvars[[varName]][, , period]),
              Matrix(obsData[[groupName]]$depvars[[varName]][, , period + 1]))
        }
      }
      1 * drop0((returnValue - missings) > 0)
    }
    
    networkExtraction <- function (i, obsData, sims, period, groupName, 
        varName) {
      dimsOfDepVar <- attr(obsData[[groupName]]$depvars[[varName]], "netdims")
      isbipartite <- (attr(obsData[[groupName]]$depvars[[varName]], "type")	
          == "bipartite")
      bipartiteOffset <- ifelse (isbipartite, 1 + dimsOfDepVar[1], 1)
      if (isbipartite) {
        emptyNetwork <- network.initialize(dimsOfDepVar[1] + dimsOfDepVar[2], 
            bipartite = dimsOfDepVar[1])
      } else {
        emptyNetwork <- network.initialize(dimsOfDepVar[1], bipartite = NULL)
      }
      matrixNetwork <- sparseMatrixExtraction(i, obsData, sims, period, 
          groupName, varName)
      sparseMatrixNetwork <- as(matrixNetwork, "dgTMatrix")
      if (sum(matrixNetwork) <= 0) {
        returnValue <- emptyNetwork
      } else {
        returnValue <- network.edgelist(
            cbind(sparseMatrixNetwork@i + 1,
            sparseMatrixNetwork@j + bipartiteOffset, 1),
            emptyNetwork
        )
      }
      returnValue
    }
  
  }
  
  # save the target object in a list and remove/handle missing data
  dvname <- attr(siena.data$depvars, "name")[1]
  dv <- eval(parse(text = dvname))
  if (!"sienaDependent" %in% class(dv) && !"sienaNet" %in% class(dv)) {
    stop(paste(dvname, "is not a sienaDependent or sienaNet object."))
  }
  dv <- dv[, , base + 1]
  missings.1 <- suppressMessages(handleMissings(dv, na = target.na, 
      method = target.na.method, logical = TRUE))
  missings.2 <- suppressMessages(handleMissings(dv, na = target.structzero, 
      method = "remove", logical = TRUE))
  missings <- missings.1
  missings[missings.2 == TRUE] <- TRUE
  dv[missings] <- NA
  dv <- suppressMessages(handleMissings(dv, na = NA, method = target.na.method, 
      logical = FALSE))
  target <- list(network::network(dv))
  
  # this function carries out one simulation at a time (for parallelization)
  simSiena <- function(q, mymodel, mydata, myeffects, mybase, mydvname, ...) {
    ans <- RSiena::siena07(mymodel, data = mydata, effects = myeffects, 
        batch = TRUE, verbose = FALSE, silent = TRUE, returnDeps = TRUE, ...)
    simul <- RSiena::networkExtraction(i = length(ans$sims), obsData = ans$f,
        sims = ans$sims, period = mybase, groupName = "Data1", 
        varName = mydvname)
    message(paste0("Completed simulation ", q, "."))
    return(simul)
  }
  
  # run the simulations, possibly in parallel
  if (parallel[1] == "snow") {
    if (is.null(cl)) {
      cl <- makeCluster(ncpus)
    }
    message(paste("Using snow parallelization with", ncpus, "cores."))
    simulations <- parLapply(cl, 1:nsim, simSiena, mymodel = object, 
        mydata = siena.data, myeffects = siena.effects, mybase = base, 
        mydvname = dvname)
    stopCluster(cl)
  } else if (parallel[1] == "multicore") {
    message(paste("Using multicore parallelization with", ncpus, "cores."))
    simulations <- mclapply(1:nsim, simSiena, mymodel = object, 
        mydata = siena.data, myeffects = siena.effects, mybase = base, 
        mydvname = dvname, mc.cores = ncpus)
  } else {
    message("Parallelization is switched off. Simulating sequentially.")
    simulations <- lapply(1:nsim, simSiena, mymodel = object, 
        mydata = siena.data, myeffects = siena.effects, mybase = base, 
        mydvname = dvname)
  }
  
  # correct directed = TRUE --> FALSE
  isbip <- sapply(simulations, network::is.bipartite)
  if (!any(isbip == TRUE)) {
    isdir <- sapply(simulations, network::is.directed)
    issym <- sapply(simulations, function(x) isSymmetric(as.matrix(x)))
    if (any(isdir && issym) && length(table(isdir && issym)) == 1 && 
        !any(sapply(target, network::is.directed))) {
      for (i in 1:length(simulations)) {
        if (issym[i] && isdir[i]) {
          simulations[[i]] <- network::as.network(as.matrix(simulations[[i]]), 
              directed = FALSE)
        }
      }
    }
  }
  
  # correct properties of the target network
  if (network::is.directed(target[[1]]) && 
      !network::is.directed(simulations[[1]])) {
    target[[1]] <- network::network(as.matrix(target[[1]]), directed = FALSE)
  } else if (!network::is.directed(target[[1]]) && 
      network::is.directed(simulations[[1]])) {
    target[[1]] <- network::network(as.matrix(target[[1]]), directed = TRUE)
  }
  if (network::is.bipartite(target[[1]]) && 
      !network::is.bipartite(simulations[[1]])) {
    target[[1]] <- network::network(as.matrix(target[[1]]), bipartite = FALSE)
  } else if (!network::is.bipartite(target[[1]]) && 
      network::is.bipartite(simulations[[1]])) {
    target[[1]] <- network::network(as.matrix(target[[1]]), bipartite = TRUE)
  }
  
  # create an object where the final results are stored
  gofobject <- list()
  class(gofobject) <- "btergmgof"
  gofobject$numbasis <- 1
  gofobject$numtarget <- 1
  gofobject$num.vertices <- simulations[[1]]$gal$n
  
  # compute classic goodness of fit
  if (classicgof == TRUE) {
    gofobject <- statnetgof(gofobject, simulations, target = target, 
        statistics = statistics)
  } else {
    gofobject$statistics <- NULL
    gofobject$raw <- NULL
  }
  
  # ROC and PR curves and AUC measures
  if (rocprgof == TRUE) {
    gofobject <- rocprgof(gofobject, simulations, target, missings, 
        rgraph = FALSE, pr.impute = pr.impute)
    randomgraphs <- randomgraph(target, nsim = nsim)
    gofobject <- rocprgof(gofobject, randomgraphs, target, rgraph = TRUE, 
        pr.impute = pr.impute)
  } else {
    gofobject$roc <- NULL
    gofobject$pr <- NULL
    gofobject$auc.roc <- NULL
    gofobject$auc.pr <- NULL
    gofobject$rgraph.roc <- NULL
    gofobject$rgraph.pr <- NULL
    gofobject$rgraph.auc.roc <- NULL
    gofobject$rgraph.auc.pr <- NULL
  }
  
  gofobject$degeneracy <- NULL
  return(gofobject)
}

# gof method for dyadic-independence models with custom data and coefficients
gof.network <- function(object, covariates, coef, target = NULL, 
    nsim = 100, mcmc = FALSE, MCMC.interval = 1000, 
    MCMC.burnin = 10000, parallel = c("no", "MPI", "SOCK"), 
    ncpus = 1, cl = NULL, classicgof = TRUE, rocprgof = TRUE, 
    statistics = c("dsp", "esp", "geodist", "degree", "idegree", 
    "odegree", "kstar", "istar", "ostar"), pr.impute = "poly4", 
    verbose = TRUE, ...) {
  
  if (nsim < 2) {
    stop("The 'nsim' argument must be greater than 1.")
  }
  
  # check dependent network
  nw <- object
  if (class(nw) == "network") {
    directed <- network::is.directed(nw)
    bipartite <- network::is.bipartite(nw)
  } else if (class(nw) == "matrix") {
    directed <- !isSymmetric(nw)
    if (nrow(nw) == ncol(nw)) {
      bipartite <- FALSE
    } else {
      bipartite <- TRUE
    }
    nw <- network::network(nw, bipartite = bipartite, directed = directed)
  } else {
    stop("'object' must be a network object or a matrix.")
  }
  time.steps <- 1
  num.vertices <- nrow(as.matrix(nw))
  
  # check and rearrange target network(s)
  if (is.null(target)) {
    message(paste("\nNo 'target' network(s) provided. Using networks on the",
        "left-hand side of the model formula as observed networks.\n"))
    target <- nw
  } else if (class(target) == "network" || class(target) == "matrix") {
    # do nothing
    message("\nOne observed ('target') network was provided.\n")
  } else if (class(target) == "list") {
    message(paste("\n", length(target), "observed ('target') networks were",
        "provided. Using the first network\n"))
    target <- target[[1]]
    if (class(target) != "matrix" && class(target) != network) {
      stop("First target network was not a matrix or network object.")
    }
  } else {
    stop("'target' must be a network, matrix, or list of matrices or networks.")
  }
  
  # check predictors and assemble formula
  if (class(covariates) != "list") {
    stop("Covariates must be provided as a list of matrices.")
  }
  numcov <- length(covariates)
  if (numcov + 1 != length(coef)) {
    stop(paste("The 'coef' vector ought to have a coefficient for edges", 
        "plus the same number of coefficients as there are covariates.", 
        "Right now, there are", length(coef), "coefficients and", numcov, 
        "covariates."))
  }
  rhs <- "edges"
  for (i in 1:numcov) {
    if (!class(covariates[[i]]) %in% c("network", "matrix")) {
      stop(paste("Covariate", i, "is not a matrix or network object."))
    }
    if (nrow(as.matrix(covariates[[i]])) != nrow(as.matrix(nw))) {
      stop(paste("Number of row nodes of covariate", i, "is not compatible."))
    }
    if (ncol(as.matrix(covariates[[i]])) != ncol(as.matrix(nw))) {
      stop(paste("Number of column nodes of covariate", i, 
          "is not compatible."))
    }
    rhs <- paste(rhs, "+ edgecov(covariates[[", i, "]])")
  }
  form <- as.formula(paste("nw ~", rhs))
  
  # simulations for statnet-style and rocpr GOF
  if (classicgof == TRUE || rocprgof == TRUE) {
    message(paste("Simulating", nsim, 
        "networks from the following formula:\n", 
        gsub("\\s+", " ", paste(deparse(form), collapse = "")), "\n"))
  if (mcmc == TRUE && parallel[1] == "no") {
    simulations <- simulate.formula(form, nsim = nsim, coef = coef, 
        constraints = ~ ., 
        control = control.simulate.formula(MCMC.interval = MCMC.interval, 
        MCMC.burnin = MCMC.burnin))
  } else if (mcmc == TRUE && parallel[1] == "SOCK" || parallel[1] == "MPI") {
    simulations <- simulate.formula(form, nsim = nsim, coef = coef, 
        constraints = ~ ., 
        control = control.simulate.formula(MCMC.interval = MCMC.interval, 
        MCMC.burnin = MCMC.burnin, parallel = ncpus, 
        parallel.type = parallel))
  } else if (mcmc == FALSE) {
    if (parallel[1] != "no") {
      warning("Parallel computation is only possible with 'mcmc = TRUE'.")
    }
    dat <- sapply(covariates, function(x) c(as.matrix(x)))
    dat <- cbind(rep(1, nrow(dat)), dat)
    prob <- plogis(coef %*% t(dat))
    simval <- t(sapply(prob, function(x) rbinom(nsim, 1, x)))
    simulations <- apply(simval, 2, function(x) network::network(matrix(x, 
        nrow = num.vertices, byrow = FALSE), bipartite = bipartite, 
        directed = directed))
  } else {
    stop(paste0("For this type of object, \"SOCK\" and \"MPI\" parallel ", 
        "processing are allowed. \"", parallel, 
        "\" is not a valid value for the \"parallel\" argument."))
    }
  }
  
  # if NA in target networks, put them in the base network, too, and vice-versa
  nw <- as.matrix(nw)
  nw[is.na(as.matrix(target))] <- NA
  nw <- network::network(nw, directed = directed, bipartite = bipartite)
  target <- as.matrix(target)
  target[is.na(as.matrix(nw))] <- NA
  target <- network::network(target, directed = directed, bipartite = bipartite)
  
  # create an object where the final results are stored
  gofobject <- list()
  class(gofobject) <- "btergmgof"
  gofobject$numbasis <- 1
  gofobject$numtarget <- 1
  gofobject$num.vertices <- num.vertices
  
  # comparison of simulated and observed network statistics (statnet-style GOF)
  target <- list(target)
  if (classicgof == TRUE) {
    gofobject <- statnetgof(gofobject, simulations, target, 
        statistics = statistics)
  } else {
    gofobject$statistics <- NULL
    gofobject$raw <- NULL
  }
  
  # ROC and PR curves and AUC measures
  if (rocprgof == TRUE) {
    gofobject <- rocprgof(gofobject, simulations, target, rgraph = FALSE, 
        pr.impute = pr.impute)
    randomgraphs <- randomgraph(list(nw), nsim = nsim)
    gofobject <- rocprgof(gofobject, randomgraphs, target, rgraph = TRUE, 
        pr.impute = pr.impute)
  } else {
    gofobject$roc <- NULL
    gofobject$pr <- NULL
    gofobject$auc.roc <- NULL
    gofobject$auc.pr <- NULL
    gofobject$rgraph.roc <- NULL
    gofobject$rgraph.pr <- NULL
    gofobject$rgraph.auc.roc <- NULL
    gofobject$rgraph.auc.pr <- NULL
  }
  
  if ((length(gofobject$auc.roc == 1) && gofobject$auc.roc >= 0.98) || 
      (length(gofobject$auc.roc > 1) && any(gofobject$auc.roc >= 0.98))) {
    warning(paste("Extremely high AUC-ROC value(s)! This may be a sign of", 
        "degenerate MCMC simulations. You should try to tweak the 'nsim',", 
        "'MCMC.interval', and 'MCMC.burnin' arguments."))
  }
  
  if ((length(gofobject$auc.pr == 1) && gofobject$auc.pr >= 0.98) || 
      (length(gofobject$auc.pr > 1) && any(gofobject$auc.pr >= 0.98))) {
    warning(paste("Extremely high AUC-PR value(s)! This may be a sign of", 
        "degenerate MCMC simulations. You should try to tweak the 'nsim',", 
        "'MCMC.interval', and 'MCMC.burnin' arguments."))
  }
  
  return(gofobject)
}


# generic methods for goodness-of-fit assessment
setMethod("gof", signature = className("btergm", "btergm"), 
    definition = gof.btergm)

setMethod("gof", signature = className("ergm", "ergm"), 
    definition = gof.btergm)

setMethod("gof", signature = className("sienaAlgorithm", "RSiena"), 
    definition = gof.sienaAlgorithm)

setMethod("gof", signature = className("sienaModel", "RSiena"), 
    definition = gof.sienaAlgorithm)

setMethod("gof", signature = className("network", "network"), 
    definition = gof.network)

setMethod("gof", signature = className("matrix", "base"), 
    definition = gof.network)


# plot method for btergmgof objects
plot.btergmgof <- function(x, boxplot = TRUE, boxplot.mfrow = TRUE, 
    boxplot.dsp = TRUE, boxplot.esp = TRUE, boxplot.geodist = TRUE, 
    boxplot.degree = TRUE, boxplot.idegree = TRUE, boxplot.odegree = TRUE, 
    boxplot.kstar = TRUE, boxplot.istar = TRUE, boxplot.ostar = TRUE, 
    boxplot.dsp.max = NULL, boxplot.esp.max = NULL, boxplot.geodist.max = NULL, 
    boxplot.degree.max = NULL, boxplot.idegree.max = NULL, 
    boxplot.odegree.max = NULL, boxplot.kstar.max = NULL, 
    boxplot.istar.max = NULL, boxplot.ostar.max = NULL, 
    boxplot.transform = function(x) x, boxplot.border = "darkgray", 
    boxplot.mean.col = "black", boxplot.median.col = "black", 
    boxplot.lwd = 0.8, boxplot.outline = FALSE, boxplot.ylab = "Frequency", 
    boxplot.main = NULL, boxplot.ylim = NULL, roc = TRUE, pr = TRUE, 
    rocpr.add = FALSE, rocpr.avg = c("none", "horizontal", "vertical", 
    "threshold"), rocpr.spread = c("boxplot", "stderror", "stddev"), 
    rocpr.lwd = 3, roc.main = NULL, roc.random = FALSE, roc.col = "#bd0017", 
    roc.random.col = "#bd001744", pr.main = NULL, pr.random = FALSE, 
    pr.col = "#5886be", pr.random.col = "#5886be44", pr.poly = 0, ...) {
  
  # boxplots of simulations vs. observed statistics (classic, statnet-style GOF)
  if (boxplot == TRUE && !is.null(x$statistics) && !is.null(x$raw)) {
    nam <- names(x$raw)
    xlength <- length(x$statistics)
    
    actual <- 0  # actual number of boxplots to be plotted
    if ("dsp" %in% nam && boxplot.dsp == TRUE) actual <- actual + 1
    if ("esp" %in% nam && boxplot.esp == TRUE) actual <- actual + 1
    if ("geodist" %in% nam && boxplot.geodist == TRUE) actual <- actual + 1
    if ("degree" %in% nam && boxplot.degree == TRUE) actual <- actual + 1
    if ("idegree" %in% nam && boxplot.idegree == TRUE) actual <- actual + 1
    if ("odegree" %in% nam && boxplot.odegree == TRUE) actual <- actual + 1
    if ("kstar" %in% nam && boxplot.kstar == TRUE) actual <- actual + 1
    if ("istar" %in% nam && boxplot.istar == TRUE) actual <- actual + 1
    if ("ostar" %in% nam && boxplot.ostar == TRUE) actual <- actual + 1
    
    if (boxplot.mfrow == TRUE) {
      if (actual == 1) {
        par(mfrow = c(1, 1))
      } else if (actual == 2) {
        par(mfrow = c(1, 2))
      } else if (actual == 3 || actual == 4) {
        par(mfrow = c(2, 2))
      } else if (actual == 5 || actual == 6) {
        par(mfrow = c(3, 2))
      } else if (actual > 6 && actual <= 9) {
        par(mfrow = c(3, 3))
      } else {
        par(mfrow = c(3, 4))
      }
    }
    
    for (i in 1:xlength) {
      mat <- t(x$raw[[i]])
      mat <- apply(mat, 1:2, boxplot.transform)
      if (ncol(x$statistics[[i]]) == 9) { # several observed networks
        obs.mean <- sapply(x$statistics[[i]][, 1], boxplot.transform)
        obs <- sapply(x$statistics[[i]][, 2], boxplot.transform)
        max.y <- max(mat, obs.mean, obs)
        min.y <- min(mat, obs.mean, obs)
      } else { # only one observed network
        obs <- sapply(x$statistics[[i]][, 1], boxplot.transform)
        max.y <- max(mat, obs)
        min.y <- min(mat, obs)
      }
      if (!is.null(boxplot.ylim)) {
        max.y <- boxplot.ylim
      }
      if (any(is.infinite(c(mat, obs)))) {
        stop(paste("Simulated or observed values contain infinite values.", 
            "Check the 'boxplot.transform' argument."))
      }
      ok <- FALSE
      if (nam[i] == "geodist" && boxplot.geodist == TRUE) {
        ok <- TRUE
        numbers <- 1:ncol(mat)
        numbers <- paste(numbers)
        numbers[length(numbers)] <- "Inf"
        colnames(mat) <- numbers
        if (!is.null(boxplot.geodist.max)) {
          mat <- mat[, c(1:boxplot.geodist.max, ncol(mat))]
          obs <- obs[c(1:boxplot.geodist.max, length(obs))]
          if (exists("obs.mean")) {
            obs.mean <- obs.mean[c(1:boxplot.geodist.max, length(obs.mean))]
          }
        }
        if (is.null(boxplot.main)) {
          boxplot(mat, ylim = c(min.y, max.y), border = boxplot.border, 
              lwd = boxplot.lwd, xlab = "Geodesic distance", 
              main = "Geodesic distance", outline = boxplot.outline, 
              ylab = boxplot.ylab, ...)
        } else {
          boxplot(mat, ylim = c(min.y, max.y), border = boxplot.border, 
              lwd = boxplot.lwd, xlab = "Geodesic distance", 
              main = boxplot.main, outline = boxplot.outline, 
              ylab = boxplot.ylab, ...)
        }
      } else {
        numbers <- 0:(ncol(mat) - 1)
        numbers <- paste(numbers)
        colnames(mat) <- numbers
        ok <- TRUE
        if (nam[i] == "dsp" && boxplot.dsp == TRUE) {
          xlab <- "Number of dyad-wise shared partners"
          main <- "Dyad-wise shared partners"
          if (!is.null(boxplot.dsp.max)) {
            mat <- mat[, c(1:(boxplot.dsp.max + 1))]
            obs <- obs[c(1:(boxplot.dsp.max + 1))]
            if (exists("obs.mean")) {
              obs.mean <- obs.mean[c(1:(boxplot.dsp.max + 1))]
            }
          }
        } else if (nam[i] == "esp" && boxplot.esp == TRUE) {
          xlab <- "Number of edge-wise shared partners"
          main <- "Edge-wise shared partners"
          if (!is.null(boxplot.esp.max)) {
            mat <- mat[, c(1:(boxplot.esp.max + 1))]
            obs <- obs[c(1:(boxplot.esp.max + 1))]
            if (exists("obs.mean")) {
              obs.mean <- obs.mean[c(1:(boxplot.esp.max + 1))]
            }
          }
        } else if (nam[i] == "degree" && boxplot.degree == TRUE) {
          xlab <- "Degree"
          main <- "Degree distribution"
          if (!is.null(boxplot.degree.max)) {
            mat <- mat[, c(1:(boxplot.degree.max + 1))]
            obs <- obs[c(1:(boxplot.degree.max + 1))]
            if (exists("obs.mean")) {
              obs.mean <- obs.mean[c(1:(boxplot.degree.max + 1))]
            }
          }
        } else if (nam[i] == "idegree" && boxplot.idegree == TRUE) {
          xlab <- "In-degree"
          main <- "In-degree distribution"
          if (!is.null(boxplot.idegree.max)) {
            mat <- mat[, c(1:(boxplot.idegree.max + 1))]
            obs <- obs[c(1:(boxplot.idegree.max + 1))]
            if (exists("obs.mean")) {
              obs.mean <- obs.mean[c(1:(boxplot.idegree.max + 1))]
            }
          }
        } else if (nam[i] == "odegree" && boxplot.odegree == TRUE) {
          xlab <- "Out-degree"
          main <- "Out-degree distribution"
          if (!is.null(boxplot.odegree.max)) {
            mat <- mat[, c(1:(boxplot.odegree.max + 1))]
            obs <- obs[c(1:(boxplot.odegree.max + 1))]
            if (exists("obs.mean")) {
              obs.mean <- obs.mean[c(1:(boxplot.odegree.max + 1))]
            }
          }
        } else if (nam[i] == "kstar" && boxplot.kstar == TRUE) {
          xlab <- "k"
          main <- "k-star distribution"
          if (!is.null(boxplot.kstar.max)) {
            mat <- mat[, c(1:(boxplot.kstar.max + 1))]
            obs <- obs[c(1:(boxplot.kstar.max + 1))]
            if (exists("obs.mean")) {
              obs.mean <- obs.mean[c(1:(boxplot.kstar.max + 1))]
            }
          }
        } else if (nam[i] == "istar" && boxplot.istar == TRUE) {
          xlab <- "In-star"
          main <- "In-star distribution"
          if (!is.null(boxplot.istar.max)) {
            mat <- mat[, c(1:(boxplot.istar.max + 1))]
            obs <- obs[c(1:(boxplot.istar.max + 1))]
            if (exists("obs.mean")) {
              obs.mean <- obs.mean[c(1:(boxplot.istar.max + 1))]
            }
          }
        } else if (nam[i] == "ostar" && boxplot.ostar == TRUE) {
          xlab <- "Out-star"
          main <- "Out-star distribution"
          if (!is.null(boxplot.ostar.max)) {
            mat <- mat[, c(1:(boxplot.ostar.max + 1))]
            obs <- obs[c(1:(boxplot.ostar.max + 1))]
            if (exists("obs.mean")) {
              obs.mean <- obs.mean[c(1:(boxplot.ostar.max + 1))]
            }
          }
        } else {
          ok <- FALSE
        }
        if (ok == TRUE) {
          if (is.null(boxplot.main)) {
            boxplot(mat, ylim = c(min.y, max.y), border = boxplot.border, 
                lwd = boxplot.lwd, xlab = xlab, ylab = boxplot.ylab, 
                main = main, outline = boxplot.outline, ...)
          } else {
            boxplot(mat, ylim = c(min.y, max.y), border = boxplot.border, 
                lwd = boxplot.lwd, xlab = xlab, ylab = boxplot.ylab, 
                main = boxplot.main, outline = boxplot.outline, ...)
          }
        }
      }
      if (ok == TRUE) {
        if (ncol(x$statistics[[i]]) == 9) {
          lines(obs.mean, lwd = 2 * boxplot.lwd, type = "l", lty = "dashed", 
              col = boxplot.mean.col)
        }
        lines(obs, lwd = 3 * boxplot.lwd, type = "l", col = boxplot.median.col)
      }
    }
    if (ncol(x$statistics[[1]]) == 9) {
      message(paste("Note: The solid line in the boxplots represents the",
          "median of the observed networks and the dashed line represents", 
          "the mean."))
    }
  } else if (boxplot == TRUE) {
    message("The object does not contain any classic GOF data.")
  }
  
  if (boxplot == TRUE && (roc == TRUE || pr == TRUE)) {
    par(ask = TRUE, mfrow = c(1, 1))
  }
  
  # receiver operating characteristics (ROC), precision recall (PR), and AUC
  if (!is.null(x$roc) && roc == TRUE && pr == FALSE) {  # plot only ROC
    plot(x$roc, avg = rocpr.avg[1], spread.estimate = rocpr.spread[1], 
        add = rocpr.add, col = roc.col, main = roc.main, lwd = rocpr.lwd, 
        boxplot.boxcol = roc.col, ylim = c(0, 1), ...)
    if (roc.random == TRUE) {
      plot(x$rgraph.roc, avg = rocpr.avg[1], spread.estimate = "none", 
          col = roc.random.col, add = TRUE, lwd = rocpr.lwd, ylim = c(0, 1), 
          ...)
    }
  } else if (!is.null(x$pr) && roc == FALSE && pr == TRUE) {  # plot only PR
    plot(x$pr, avg = rocpr.avg[1], spread.estimate = rocpr.spread[1], 
        add = rocpr.add, col = pr.col, main = pr.main, lwd = rocpr.lwd, 
        boxplot.boxcol = pr.col, ylim = c(0, 1), ...)
    if (pr.random == TRUE) {
      plot(x$rgraph.pr, avg = rocpr.avg[1], spread.estimate = "none", 
          col = pr.random.col, add = TRUE, lwd = rocpr.lwd, ylim = c(0, 1), ...)
    }
  } else if (!is.null(x$roc) && !is.null(x$pr) && roc == TRUE && pr == TRUE) {
    plot(x$roc, avg = rocpr.avg[1], spread.estimate = rocpr.spread[1], 
        add = rocpr.add, col = roc.col, main = roc.main, lwd = rocpr.lwd, 
        boxplot.boxcol = roc.col, ylim = c(0, 1), ...)
    if (roc.random == TRUE) {
      plot(x$rgraph.roc, avg = rocpr.avg[1], spread.estimate = "none", 
          col = roc.random.col, add = TRUE, lwd = rocpr.lwd, ylim = c(0, 1), 
          ...)
    }
    par(ask = TRUE)
    plot(x$pr, avg = rocpr.avg[1], spread.estimate = rocpr.spread[1], 
        add = rocpr.add, col = pr.col, main = pr.main, lwd = rocpr.lwd, 
        boxplot.boxcol = pr.col, ylim = c(0, 1), ...)
    if (pr.random == TRUE) {
      plot(x$rgraph.pr, avg = rocpr.avg[1], spread.estimate = "none", 
          col = pr.random.col, add = TRUE, lwd = rocpr.lwd, ylim = c(0, 1), ...)
    }
  }
  if (pr == TRUE && !is.null(x$pr) && pr.poly > 0) {
    for (i in 1:length(x$pr@y.values)) {
      prec <- x$pr@y.values[[i]]
      prec[1] <- NaN
      rec <- x$pr@x.values[[i]]
      p <- data.frame(poly(rec, pr.poly, raw = TRUE))
      fit <- lm(prec ~ ., data = p)
      yhat <- predict(fit, newdata = p)
      for (j in 1:length(yhat)) {
        if (yhat[j] > 1) {
          yhat[j] <- 1
        }
        if (yhat[j] < 0) {
          yhat[j] <- 0
        }
      }
      lines(rec, yhat, type = "l", lty = "dashed", col = "red", 
          lwd = rocpr.lwd)
      points(rec[1], yhat[1], pch = 1, col = "red", cex = 2)
      points(x$pr@x.values[[i]][1], x$pr@y.values[[i]][1], pch = 1, 
          col = pr.col, cex = 2)
    }
  }
}


# print method for btergmgof objects for nicely formatted output
print.btergmgof <- function(x, classicgof = TRUE, rocprgof = TRUE, 
    degeneracy = TRUE, ...) {
  
  # print statnet-like GOF output
  if (classicgof == TRUE && !is.null(x$statistics) && !is.null(x$raw)) {
    nam <- names(x$statistics)
    for (j in 1:length(x$statistics)) {
      mat <- x$statistics[[j]]
      if (nam[j] == "dsp") {
        message("Goodness of fit for dyad-wise shared partners:")
      } else if (nam[j] == "esp") {
        message("Goodness of fit for edge-wise shared partners:")
      } else if (nam[j] == "geodist") {
        message("Goodness of fit for geodesic distances:")
      } else if (nam[j] == "degree") {
        message("Goodness of fit for degree:")
      } else if (nam[j] == "idegree") {
        message("Goodness of fit for in-degree:")
      } else if (nam[j] == "odegree") {
        message("Goodness of fit for out-degree:")
      } else {
        message(paste("Goodness of fit for ", nam[j], ":", sep = ""))
      }
      if (ncol(x$statistics[[1]]) == 9) {
        printCoefmat(mat, has.Pvalue = TRUE, cs.ind=numeric(0), 
            signif.legend = FALSE)
      } else {
        printCoefmat(mat, has.Pvalue = TRUE, cs.ind=numeric(0), 
            signif.legend = FALSE)
      }
      message("")
    }
    message(paste("Note: Small p values indicate a significant difference",
        "between\n      simulations and observed network(s)."))
    if (ncol(x$statistics[[1]]) == 9) {
      message("      The p values are based on a two-sample t-test.")
    }
    message(paste("      Results from ", x$numbasis, 
        " 'basis' network(s) (with ", ncol(x$raw[[1]]) / x$numbasis, 
        " simulations\n      each) and ", x$numtarget, 
        " observed 'target' network(s).", sep = ""))
  }
  
  # print ROC AUC
  if (!is.null(x$auc.roc) && rocprgof == TRUE) {
    message("\nArea under the ROC/PR curve:")
    mat <- cbind(t = 1:length(x$auc.roc), "AUC ROC" = x$auc.roc, 
        "AUC PR" = x$auc.pr)
    row.names(mat) <- rep("", length(x$auc.roc))
    printCoefmat(mat, tst.ind = 2:3)
  }
  
  # print degeneracy check
  if (degeneracy == TRUE && !is.null(x$degeneracy)) {
    for (i in 1:length(x$degeneracy)) {
      message(paste0("\nDegeneracy check for network ", i, ":"))
      printCoefmat(x$degeneracy[[i]], digits = 3, P.values = TRUE, 
          has.Pvalue = TRUE, cs.ind = 3:4, ts.ind = 5)
    }
  }
}

# summary method for btergm.gof objects
summary.btergmgof <- function(object, classicgof = TRUE, rocprgof = TRUE, 
    degeneracy = TRUE, ...) {
  print.btergmgof(object, classicgof = classicgof, rocprgof = rocprgof, 
      degeneracy = degeneracy, ...)
}

