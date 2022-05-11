#' Micro-Level Interpretation of (T)ERGMs
#'
#' Micro-level interpretation of (T)ERGMs.
#'
#' The \code{interpret} function facilitates interpretation of ERGMs and TERGMs
#' at the micro level, as described in Desmarais and Cranmer (2012). There are
#' methods for \code{ergm} objects, \code{btergm} objects, and \code{mtergm}
#' objects. The function can be used to interpret these models at the tie or
#' edge level, dyad level, and block level. For example, what is the probability
#' that two specific nodes \code{i} (the sender) and \code{j} (the receiver) are
#' connected given the rest of the network and given the model? Or what is the
#' probability that any two nodes are tied at \code{t = 2} if they were tied (or
#' disconnected) at \code{t = 1} (i.e., what is the amount of tie stability)?
#' These tie- or edge-level questions can be answered if the \code{type = "tie"}
#' argument is used.
#'
#' Another example: What is the probability that node \code{i} has a tie to node
#' \code{j} but not vice-versa? Or that \code{i} and \code{j} maintain a
#' reciprocal tie? Or that they are disconnected? How much more or less likely
#' are \code{i} and \code{j} reciprocally connected if the \code{mutual} term in
#' the model is fixed at \code{0} (compared to the model that includes the
#' estimated parameter for reciprocity)? See example below. These dyad-level
#' questions can be answered if the \code{type = "dyad"} argument is used.
#'
#' Or what is the probability that a specific node \code{i} is connected to
#' nodes \code{j1} and \code{j2} but not to \code{j5} and \code{j7}? And how
#' likely is any node \code{i} to be connected to exactly four \code{j} nodes?
#' These node-level questions (focusing on the ties of node \code{i} or node
#' \code{j}) can be answered by using the \code{type = "node"} argument.
#'
#' The typical procedure is to manually enumerate all dyads or
#' sender-receiver-time combinations with certain properties and repeat the same
#' thing with some alternative properties for contrasting the two groups. Then
#' apply the \code{interpret} function to the two groups of dyads and compute a
#' measure of central tendency (e.g., mean or median) and possibly some
#' uncertainy measure (i.e., confidence intervals) from the distribution of
#' dyadic probabilities in each group. For example, if there is a gender
#' attribute, one can sample male-male or female-female dyads, compute the
#' distributions of edge probabilities for the two sets of dyads, and create
#' boxplots or barplots with confidence intervals for the two types of dyads in
#' order to contrast edge probabilities for male versus female same-sex dyads.
#'
#' See also the \code{\link{edgeprob}} function for automatic computation of all
#' dyadic edge probabilities.
#'
#' @param object An \code{ergm}, \code{btergm}, or \code{mtergm} object.
#' @param formula The formula to be used for computing probabilities. By
#'   default, the formula embedded in the model object is retrieved and used.
#' @param coefficients The estimates on which probabilities should be based. By
#'   default, the coefficients from the model object are retrieved and used.
#'   Custom coefficients can be handed over, for example, in order to compare
#'   versions of the model where the reciprocity term is fixed at \code{0}
#'   versus versions of the model where the reciprocity term is left as in the
#'   empirical result. This is one of the examples described in Desmarais and
#'   Cranmer (2012).
#' @param target The response network on which probabilities are based.
#'   Depending on whether the function is applied to an \code{ergm} or
#'   \code{btergm}/\code{mtergm} object, this can be either a single network or
#'   a list of networks. By default, the (list of) network(s) provided as the
#'   left-hand side of the (T)ERGM formula is used.
#' @param type If \code{type = "tie"} is used, probabilities at the edge level
#'   are computed. For example, what is the probability of a specific node
#'   \code{i} to be connected to a specific node \code{j} given the rest of the
#'   network and given the model? If \code{type = "dyad"} is used, probabilities
#'   at the dyad level are computed. For example, what is the probability that
#'   node \code{i} is connected to node \code{j} but not vice-versa, or what is
#'   the probability that nodes \code{i} and \code{j} and mutually connected in
#'   a directed network? If \code{type = "node"} is used, probabilities at the
#'   node level are computed. For example, what is the probability that node
#'   \code{i} is connected to a set of three other \code{j} nodes given the rest
#'   of the network and the model?
#' @param i A single (sender) node \code{i} or a set of (sender) nodes \code{i}.
#'   If \code{type = "node"} is used, this can be more than one node and should
#'   be provided as a vector. The \code{i} argument can be either provided as
#'   the index of the node in the sociomatrix (e.g., the fourth node would be
#'   \code{i = 4}) or the row name of the node in the sociomatrix (e.g.,
#'   \code{i = "Peter"}). If more than one node is provided and
#'   \code{type = "node"}, there can be only one (receiver) node \code{j}. The
#'   \code{i} and \code{j} arguments are used to specify for which nodes
#'   probabilities should be computed. For example, what is the probability that
#'   \code{i = 4} is connected to \code{j = 7}?
#' @param j A single (receiver) node \code{j} or a set of (receiver) nodes
#'   \code{j}. If \code{type = "node"} is used, this can be more than one node
#'   and should be provided as a vector. The \code{j} argument can be either
#'   provided as the index of the node in the sociomatrix (e.g., the fourth node
#'   would be \code{j = 4}) or the column name of the node in the sociomatrix
#'   (e.g., \code{j = "Mary"}). If more than one node is provided and
#'   \code{type = "node"}, there can be only one (sender) node \code{i}. The
#'   \code{i} and \code{j} arguments are used to specify for which nodes
#'   probabilities should be computed. For example, what is the probability that
#'   \code{i = 4} is connected to \code{j = 7}?
#' @param t A vector of (numerical) time steps for which the probabilities
#'   should be computed. This only applies to \code{btergm} amd \code{mtergm}
#'   objects because \code{ergm} objects are by definition based on a single
#'   time step. By default, all available time steps are used. It is, for
#'   example, possible to compute probabilities only for a single time step by
#'   specifying, e.g., \code{t = 5} in order to compute probabilities for the
#'   fifth response network.
#' @param ... Further arguments to be passed on to subroutines.
#'
#' @references
#' Desmarais, Bruce A. and Skyler J. Cranmer (2012): Micro-Level Interpretation
#' of Exponential Random Graph Models with Application to Estuary Networks.
#' \emph{Policy Studies Journal} 40(3): 402--434.
#' \doi{10.1111/j.1541-0072.2012.00459.x}.
#'
#' Leifeld, Philip, Skyler J. Cranmer and Bruce A. Desmarais (2017): Temporal
#' Exponential Random Graph Models with btergm: Estimation and Bootstrap
#' Confidence Intervals. \emph{Journal of Statistical Software} 83(6): 1--36.
#' \doi{10.18637/jss.v083.i06}.
#'
#' Czarna, Anna Z., Philip Leifeld, Magdalena Smieja, Michael Dufner and Peter
#' Salovey (2016): Do Narcissism and Emotional Intelligence Win Us Friends?
#' Modeling Dynamics of Peer Popularity Using Inferential Network Analysis.
#' \emph{Personality and Social Psychology Bulletin} 42(11): 1588--1599.
#' \doi{10.1177/0146167216666265}.
#'
#' @examples
#' ##### The following example is a TERGM adaptation of the #####
#' ##### dyad-level example provided in figure 5(c) on page #####
#' ##### 424 of Desmarais and Cranmer (2012) in the PSJ. At #####
#' ##### each time step, it compares dyadic probabilities   #####
#' ##### (no tie, unidirectional tie, and reciprocal tie    #####
#' ##### probability) between a fitted model and a model    #####
#' ##### where the reciprocity effect is fixed at 0 based   #####
#' ##### on 20 randomly selected dyads per time step. The   #####
#' ##### results are visualized using a grouped bar plot.   #####
#'
#' \dontrun{
#'   # create toy dataset and fit a model
#'   networks <- list()
#'   for (i in 1:3) {           # create 3 random networks with 10 actors
#'     mat <- matrix(rbinom(100, 1, 0.25), nrow = 10, ncol = 10)
#'     diag(mat) <- 0           # loops are excluded
#'     nw <- network(mat)       # create network object
#'     networks[[i]] <- nw      # add network to the list
#'   }
#'   fit <- btergm(networks ~ edges + istar(2) + mutual, R = 200)
#'
#'   # extract coefficients and create null hypothesis vector
#'   null <- coef(fit)  # estimated coefs
#'   null[3] <- 0       # set mutual term = 0
#'
#'   # sample 20 dyads per time step and compute probability ratios
#'   probabilities <- matrix(nrow = 9, ncol = length(networks))
#'   # nrow = 9 because three probabilities + upper and lower CIs
#'   colnames(probabilities) <- paste("t =", 1:length(networks))
#'   for (t in 1:length(networks)) {
#'     d <- dim(as.matrix(networks[[t]]))  # how many row and column nodes?
#'     size <- d[1] * d[2]                 # size of the matrix
#'     nw <- matrix(1:size, nrow = d[1], ncol = d[2])
#'     nw <- nw[lower.tri(nw)]             # sample only from lower triangle b/c
#'     samp <- sample(nw, 20)              # dyadic probabilities are symmetric
#'     prob.est.00 <- numeric(0)
#'     prob.est.01 <- numeric(0)
#'     prob.est.11 <- numeric(0)
#'     prob.null.00 <- numeric(0)
#'     prob.null.01 <- numeric(0)
#'     prob.null.11 <- numeric(0)
#'     for (k in 1:20) {
#'       i <- arrayInd(samp[k], d)[1, 1]   # recover 'i's and 'j's from sample
#'       j <- arrayInd(samp[k], d)[1, 2]
#'       # run interpretation function with estimated coefs and mutual = 0:
#'       int.est <- interpret(fit, type = "dyad", i = i, j = j, t = t)
#'       int.null <- interpret(fit, coefficients = null, type = "dyad",
#'                             i = i, j = j, t = t)
#'       prob.est.00 <- c(prob.est.00, int.est[[1]][1, 1])
#'       prob.est.11 <- c(prob.est.11, int.est[[1]][2, 2])
#'       mean.est.01 <- (int.est[[1]][1, 2] + int.est[[1]][2, 1]) / 2
#'       prob.est.01 <- c(prob.est.01, mean.est.01)
#'       prob.null.00 <- c(prob.null.00, int.null[[1]][1, 1])
#'       prob.null.11 <- c(prob.null.11, int.null[[1]][2, 2])
#'       mean.null.01 <- (int.null[[1]][1, 2] + int.null[[1]][2, 1]) / 2
#'       prob.null.01 <- c(prob.null.01, mean.null.01)
#'     }
#'     prob.ratio.00 <- prob.est.00 / prob.null.00  # ratio of est. and null hyp
#'     prob.ratio.01 <- prob.est.01 / prob.null.01
#'     prob.ratio.11 <- prob.est.11 / prob.null.11
#'     probabilities[1, t] <- mean(prob.ratio.00)   # mean estimated 00 tie prob
#'     probabilities[2, t] <- mean(prob.ratio.01)   # mean estimated 01 tie prob
#'     probabilities[3, t] <- mean(prob.ratio.11)   # mean estimated 11 tie prob
#'     ci.00 <- t.test(prob.ratio.00, conf.level = 0.99)$conf.int
#'     ci.01 <- t.test(prob.ratio.01, conf.level = 0.99)$conf.int
#'     ci.11 <- t.test(prob.ratio.11, conf.level = 0.99)$conf.int
#'     probabilities[4, t] <- ci.00[1]              # lower 00 conf. interval
#'     probabilities[5, t] <- ci.01[1]              # lower 01 conf. interval
#'     probabilities[6, t] <- ci.11[1]              # lower 11 conf. interval
#'     probabilities[7, t] <- ci.00[2]              # upper 00 conf. interval
#'     probabilities[8, t] <- ci.01[2]              # upper 01 conf. interval
#'     probabilities[9, t] <- ci.11[2]              # upper 11 conf. interval
#'   }
#'
#'   # create barplots from probability ratios and CIs
#'   require("gplots")
#'   bp <- barplot2(probabilities[1:3, ], beside = TRUE, plot.ci = TRUE,
#'                  ci.l = probabilities[4:6, ], ci.u = probabilities[7:9, ],
#'                  col = c("tan", "tan2", "tan3"), ci.col = "grey40",
#'                  xlab = "Dyadic tie values", ylab = "Estimated Prob./Null Prob.")
#'   mtext(1, at = bp, text = c("(0,0)", "(0,1)", "(1,1)"), line = 0, cex = 0.5)
#'
#'
#'   ##### The following examples illustrate the behavior of  #####
#'   ##### the interpret function with undirected and/or      #####
#'   ##### bipartite graphs with or without structural zeros. #####
#'
#'   library("statnet")
#'   library("btergm")
#'
#'   # micro-level interpretation for undirected network with structural zeros
#'   set.seed(12345)
#'   mat <- matrix(rbinom(400, 1, 0.1), nrow = 20, ncol = 20)
#'   mat[1, 5] <- 1
#'   mat[10, 7] <- 1
#'   mat[15, 3] <- 1
#'   mat[18, 4] < 1
#'   nw <- network(mat, directed = FALSE, bipartite = FALSE)
#'   cv <- matrix(rnorm(400), nrow = 20, ncol = 20)
#'   offsetmat <- matrix(rbinom(400, 1, 0.1), nrow = 20, ncol = 20)
#'   offsetmat[1, 5] <- 1
#'   offsetmat[10, 7] <- 1
#'   offsetmat[15, 3] <- 1
#'   offsetmat[18, 4] < 1
#'   model <- ergm(nw ~ edges + kstar(2) + edgecov(cv) + offset(edgecov(offsetmat)),
#'                 offset.coef = -Inf)
#'   summary(model)
#'
#'   # tie-level interpretation (note that dyad interpretation would not make any
#'   # sense in an undirected network):
#'   interpret(model, type = "tie", i = 1, j = 2)  # 0.28 (= normal dyad)
#'   interpret(model, type = "tie", i = 1, j = 5)  # 0.00 (= structural zero)
#'
#'   # node-level interpretation; note the many 0 probabilities due to the
#'   # structural zeros; also note the warning message that the probabilities may
#'   # be slightly imprecise because -Inf needs to be approximated by some large
#'   # negative number (-9e8):
#'   interpret(model, type = "node", i = 1, j = 3:5)
#'
#'   # repeat the same exercise for a directed network
#'   nw <- network(mat, directed = TRUE, bipartite = FALSE)
#'   model <- ergm(nw ~ edges + istar(2) + edgecov(cv) + offset(edgecov(offsetmat)),
#'                 offset.coef = -Inf)
#'   interpret(model, type = "tie", i = 1, j = 2)  # 0.13 (= normal dyad)
#'   interpret(model, type = "tie", i = 1, j = 5)  # 0.00 (= structural zero)
#'   interpret(model, type = "dyad", i = 1, j = 2)  # results for normal dyad
#'   interpret(model, type = "dyad", i = 1, j = 5)  # results for i->j struct. zero
#'   interpret(model, type = "node", i = 1, j = 3:5)
#'
#'   # micro-level interpretation for bipartite graph with structural zeros
#'   set.seed(12345)
#'   mat <- matrix(rbinom(200, 1, 0.1), nrow = 20, ncol = 10)
#'   mat[1, 5] <- 1
#'   mat[10, 7] <- 1
#'   mat[15, 3] <- 1
#'   mat[18, 4] < 1
#'   nw <- network(mat, directed = FALSE, bipartite = TRUE)
#'   cv <- matrix(rnorm(200), nrow = 20, ncol = 10)  # some covariate
#'   offsetmat <- matrix(rbinom(200, 1, 0.1), nrow = 20, ncol = 10)
#'   offsetmat[1, 5] <- 1
#'   offsetmat[10, 7] <- 1
#'   offsetmat[15, 3] <- 1
#'   offsetmat[18, 4] < 1
#'   model <- ergm(nw ~ edges + b1star(2) + edgecov(cv)
#'                 + offset(edgecov(offsetmat)), offset.coef = -Inf)
#'   summary(model)
#'
#'   # tie-level interpretation; note the index for the second mode starts with 21
#'   interpret(model, type = "tie", i = 1, j = 21)
#'
#'   # dyad-level interpretation does not make sense because network is undirected;
#'   # node-level interpretation prints warning due to structural zeros, but
#'   # computes the correct probabilities (though slightly imprecise because -Inf
#'   # is approximated by some small number:
#'   interpret(model, type = "node", i = 1, j = 21:25)
#'
#'   # compute all dyadic probabilities
#'   dyads <- edgeprob(model)
#'   dyads
#' }
#'
#' @docType methods
#' @aliases interpret-methods
#' @family interpretation
#' @export
setGeneric("interpret", function(object, ...) standardGeneric("interpret"),
           package = "btergm")

#' Internal helper function for doing the actual interpretation computations
#'
#' Internal helper function for doing the actual interpretation computations.
#'
#' Takes the output of \code{\link{tergmprepare}} and a few parameters to
#' produce micro-level interpretation results. This is an internal worker
#' function. Users should call the high-level functions \code{interpret} and
#' \code{edgeprob} instead.
#'
#' @param l List of networks, with possibly just one network in this list, and a
#'   formula (\code{l$form}, with temporal indices, as prepared by
#'   \code{\link{tergmprepare}}).
#' @param type Type string. Can be \code{"tie"}, \code{"dyad"}, or
#'   \code{"node"}.
#'
#' @importFrom utils combn
#' @importFrom network is.network is.directed
#'
#' @noRd
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

#' @noRd
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

#' @noRd
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

#' @describeIn interpret Interpret method for \code{ergm} objects
setMethod("interpret", signature = className("ergm", "ergm"),
    definition = interpret.ergm)

#' @describeIn interpret Interpret method for \code{btergm} objects
setMethod("interpret", signature = className("btergm", "btergm"),
    definition = interpret.btergm)

#' @describeIn interpret Interpret method for \code{mtergm} objects
setMethod("interpret", signature = className("mtergm", "btergm"),
    definition = interpret.btergm)

#' Create all predicted tie probabilities using MPLE
#'
#' Create all predicted tie probabilities using MPLE.
#'
#' For a given (T)ERGM, return a data frame with all predicted edge
#' probabilities along with the design matrix of the MPLE logit model, based
#' on the estimated coefficients and the design matrix, for all time points,
#' along with \code{i}, \code{j}, and \code{t} variables indicating where the
#' respective dyad is located.
#'
#' \code{edgeprob} is a convenience function that creates a data frame with all
#' dyads in the ERGM or TERGM along with their edge probabilities and their
#' predictor values (i.e., change statistics). This is useful for creating
#' marginal effects plots or contrasting multiple groups of dyads. This function
#' works faster than the \code{\link{interpret}} function.
#'
#' @param verbose Print details?
#' @inheritParams interpret
#'
#' @return The first variable in the resulting data frame contains the edge
#' value (i.e., the dependent variable, which is usually binary). The next
#' variables contain all the predictors from the ERGM or TERGM (i.e., the change
#' statistics). The next five variables contain the indices of the sender (i),
#' the receiver (j), the time step (t), the vertex id of i (i.name), and the
#' vertex id of j (j.name). These five variables serve to identify the dyad. The
#' last variable contains the computed edge probabilities.
#'
#' @aliases edgeprob
#' @family interpretation
#'
#' @importFrom network network is.network is.bipartite
#' @export
edgeprob <- function (object, verbose = FALSE) {
  if ("ergm" %in% class(object)) {
    tergm <- FALSE
  } else if ("btergm" %in% class(object) || "mtergm" %in% class(object)) {
    tergm <- TRUE
  } else {
    stop(paste("The edgeprob function is only applicable to ergm, btergm, and",
               "mtergm objects."))
  }
  if (any(grepl("^gw.+#\\d{1,}$", names(coef(object))))) {
    stop("MPLE-based (T)ERGMs with variable GW* decay are currently not supported.")
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
    if (any(grepl("^gw.+#\\d{1,}$", colnames(mpli$predictor)))) {
      stop("MPLE-based (T)ERGMs with variable GW* decay are currently not supported.")
    }
    mpli$predictor <- cbind(mpli$predictor, i)
    # TODO: the previous line might say something like: "Note: Term nodeofactor("edu") skipped because it contributes no statistics." when a covariate is full of zeros.

    # retrieve row and column names
    mat.i <- as.matrix(l$networks[[i]])
    i_ind <- mpli$predictor[, which(colnames(mpli$predictor) == "edgecov.imat")]
    j_ind <- mpli$predictor[, which(colnames(mpli$predictor) == "edgecov.jmat")]
    i.names <- rownames(mat.i)[i_ind]
    if (isTRUE(l$bipartite)) {
      j.names <- colnames(mat.i)[j_ind - nrow(mat.i)]
    } else {
      j.names <- colnames(mat.i)[j_ind]
    }

    if (is.null(dyads)) {
      dyads <- cbind(as.data.frame(mpli$predictor), i.names, j.names)
    } else { # before rbind can be use, we need to check if some covariate levels were absent either before or now
      newColumns <- which(sapply(colnames(mpli$predictor), function(x) x %in% colnames(dyads)) == FALSE)
      if (length(newColumns) > 0) {
        for (j in 1:length(newColumns)) {
          dyads <- cbind(dyads[, 1:(newColumns[j] - 1)], 0, dyads[, newColumns[j]:ncol(dyads)])
          colnames(dyads)[newColumns[j]] <- names(newColumns)[j]
        }
      }
      notPresentAnymore <- which(sapply(colnames(dyads)[1:(ncol(dyads) - 2)], function(x) x %in% colnames(mpli$predictor)) == FALSE)
      if (length(notPresentAnymore) > 0) {
        for (j in 1:length(notPresentAnymore)) {
          mpli$predictor <- cbind(mpli$predictor[, 1:(notPresentAnymore[j] - 1)], 0, mpli$predictor[, notPresentAnymore[j]:ncol(mpli$predictor)])
          colnames(mpli$predictor)[notPresentAnymore[j]] <- names(notPresentAnymore)[j]
        }
      }
      dyads <- rbind(dyads, cbind(as.data.frame(mpli$predictor), i.names, j.names))
    }
  }
  term.names <- colnames(dyads)[1:(ncol(dyads) - 5)]
  term.names <- c(term.names, "i", "j", "t", "i.name", "j.name")
  colnames(dyads) <- term.names
  dyads <- cbind(Y, dyads)
  colnames(dyads)[1] <- "tie"
  class(dyads$i) <- "integer"
  class(dyads$j) <- "integer"
  class(dyads$t) <- "integer"
  cf <- coef(object)
  cf.length <- length(cf)
  cf <- cf[!cf %in% c(Inf, -Inf)]
  if (length(cf) != cf.length) {
    warning(paste("There are structural zeros or ones. For these dyads, the",
                  "predicted probabilities are not valid and must be manually replaced",
                  "by 0 or 1, respectively."))
  }
  cbcoef <- cbind(cf)
  chgstat <- dyads[, 2:(ncol(dyads) - 5)]

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
  result <- c(1 / (1 + exp(-lp)))
  dyads <- cbind(dyads, result)
  colnames(dyads)[ncol(dyads)] <- "probability"
  dyads <- dyads[order(dyads$t, dyads$i, dyads$j), ]
  rownames(dyads) <- NULL
  return(dyads)
}

#' Plot marginal effects for two-way interactions in (T)ERGMs
#'
#' Plot marginal effects for two-way interactions in (T)ERGMs.
#'
#' The \code{marginalplot} function creates marginal effects plots for ERGMs
#' with interaction effects. The user has to supply the \code{ergm} object and
#' the coefficient names of the first main variable, the second main variable,
#' and the interaction term as stored in the coefficients vector inside the
#' \code{ergm} object. It is possible to draw continuous curves or discrete
#' error bars depending on the nature of the data (using the \code{point}
#' argument). The distribution of the second (conditioning) variable can be
#' plotted at the bottom of the viewport using the \code{rug} argument.
#'
#' The resulting marginal effects plot is a \code{ggplot2} plot. This means it
#' can be extended by plotting additional elements and using themes.
#'
#' @param model An \code{ergm} object as generated by the \pkg{ergm} package.
#'   Note that marginal effects plots cannot be created for \code{btergm}
#'   objects because the variance-covariance matrix is not valid. However, it
#'   should be possible to apply the \code{marginalplot} function to
#'   MCMC-MLE-estimated TERGMs because the \code{ergm} object is stored in the
#'   \code{ergm} slot of an \code{mtergm} object. To do this, supply the
#'   \code{ergm} object instead of the \code{mtergm} object (e.g.,
#'   \code{marginalplot(mtergmobject@ergm)}).
#' @param var1 Name of the first main variable. This is the focal variable.
#' @param var2 Name of the second main variable. This is the conditioning
#'   variable.
#' @param inter Name of the interaction effect.
#' @param ci Significance level.
#' @param rug Display the distribution of the conditioning variable at the
#'   bottom of the plot?
#' @param point Display error bars for the levels of the conditioning variable
#'   (instead of a continuous curve)?
#' @param structzeromat An optional matrix object which indicates dyads that
#'   should be deleted prior to the calculation of the confidence interval for
#'   the marginal effect curve. This is useful when such a matrix was used to
#'   indicate structural zeros during estimation. In this event, the dyads
#'   characterized by structural zeros are not allowed to be tied, therefore
#'   they should be removed from the set of dyads used for the calculation of
#'   marginal effects. The matrix should contain ones for structural zeros and
#'   zeros for entries that should be used.
#' @param zeroline Draw a horizontal line to indicate zero for the first main
#' variable?
#' @param color Color of the curve, confidence interval, and distribution.
#' @param xlab Axis label for the second (conditioning) variable.
#' @param ylab Axis label for the first (focal) variable.
#'
#' @examples
#' \dontrun{
#' # data preparation
#' data("florentine")
#' n <- network.size(flobusiness)
#' wealth <- get.vertex.attribute(flobusiness, "wealth")
#' priorates <- get.vertex.attribute(flobusiness, "priorates")
#' wealth.icov <- matrix(rep(wealth, n), ncol = n, byrow = TRUE)
#' priorates.icov <- matrix(rep(priorates, n), ncol = n, byrow = TRUE)
#' interac <- wealth.icov * priorates.icov
#'
#' # estimate model with interaction effect
#' model <- ergm(flobusiness ~ edges + esp(1) + edgecov(wealth.icov)
#'                 + edgecov(priorates.icov) + edgecov(interac))
#'
#' # plot the interaction (note the additional optional ggplot2 elements)
#' marginalplot(model, var1 = "edgecov.wealth.icov",
#'              var2 = "edgecov.priorates.icov", inter = "edgecov.interac",
#'              color = "darkred", rug = TRUE, point = FALSE,
#'              xlab = "Priorates", ylab = "Wealth") +
#'   ggplot2::theme_bw() +
#'   ggplot2::ggtitle("Interaction effect")
#' }
#'
#' @family interpretation
#' @export
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
    gp <- ggplot2::ggplot(data = dta, ggplot2::aes(x = v2, y = delta1)) +
      ggplot2::geom_errorbar(data = dta,
                             ggplot2::aes(x = v2, ymin = lower, ymax = upper),
                             color = color) +
      ggplot2::geom_point(data = dta, ggplot2::aes(x = v2, y = delta1)) +
      ggplot2::ylab(ylab) +
      ggplot2::xlab(xlab)
  } else {
    gp <- ggplot2::ggplot(data = dta, ggplot2::aes(x = v2, y = delta1)) +
      ggplot2::geom_line(color = color) +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = "lower", ymax = "upper"),
                  alpha = 0.15,
                  fill = color) +
      ggplot2::ylab(ylab) +
      ggplot2::xlab(xlab)
  }

  # add distribution to the plot
  if (rug == TRUE) {
    gp <- gp + ggplot2::geom_rug(data = ep[var2],
                                 ggplot2::aes_string(x = var2),
                                 inherit.aes = FALSE,
                                 sides = "b",
                                 col = color,
                                 alpha = 0.1)
  }

  # add horizontal line for 0
  if (zeroline == TRUE) {
    gp <- gp + ggplot2::geom_hline(yintercept = 0, linetype = "dashed")
  }

  return(gp)
}