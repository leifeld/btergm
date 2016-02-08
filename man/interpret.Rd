\name{interpret}
\docType{methods}
\alias{interpret}
\alias{interpret.ergm}
\alias{interpret.btergm}
\alias{interpret.mtergm}
\alias{interpretation}
\alias{interpret-methods}
\alias{interpret,ergm-method}
\alias{interpret,btergm-method}
\alias{interpret,mtergm-method}
\title{Interpretation functions for ergm and btergm objects}
\description{Interpretation functions for ergm and btergm objects.}
\usage{
\S4method{interpret}{ergm}(object, formula = object$formula, 
    coefficients = coef(object), target = eval(parse(text = 
    deparse(formula[[2]]))), type = "tie", i, j, ...)

\S4method{interpret}{btergm}(object, formula = getformula(object), 
    coefficients = coef(object), target = NULL, type = "tie", i, j, 
    t = 1:object@time.steps, ...)

\S4method{interpret}{mtergm}(object, formula = getformula(object), 
    coefficients = coef(object), target = NULL, type = "tie", i, j, 
    t = 1:object@time.steps, ...)
}
\details{
The \code{interpret} function facilitates interpretation of ERGMs and TERGMs at the micro level via block Gibbs sampling, as described in Desmarais and Cranmer (2012). There are generic methods for \code{ergm} objects, \code{btergm} objects, and \code{mtergm} objects. The function can be used to interpret these models at the tie or edge level, dyad level, and block level.

For example, what is the probability that two specific nodes i (the sender) and node j (the receiver) are connected given the rest of the network and given the model? Or what is the probability that any two nodes are tied at t = 2 if they were tied (or disconnected) at t = 1 (i.e., what is the amount of tie stability)? These tie- or edge-level questions can be answered if the \code{type = "tie"} argument is used.

Another example: What is the probability that node i has a tie to node j but not vice-versa? Or that i and j maintain a reciprocal tie? Or that they are disconnected? How much more or less likely are i and j reciprocally connected if the \code{mutual} term in the model is fixed at 0 (compared to the model that includes the estimated parameter for reciprocity)? See example below. These dyad-level questions can be answered if the \code{type = "dyad"} argument is used.

Or what is the probability that a specific node i is connected to nodes j1 and j2 but not to j5 and j7? And how likely is any node i to be connected to exactly four j nodes? These node-level questions (focusing on the ties of node i or node j) can be answered by using the \code{type = "node"} argument.
}
\arguments{
\item{object}{ An \code{ergm} or \code{btergm} object. }
\item{formula}{ The formula to be used for computing probabilities. By default, the formula embedded in the model object is retrieved and used. }
\item{coefficients}{ The estimates on which probabilities should be based. By default, the coefficients from the model object are retrieved and used. Custom coefficients can be handed over, for example, in order to compare versions of the model where the reciprocity term is fixed at 0 versus versions of the model where the reciprocity term is left as in the empirical result. This is one of the examples described in Desmarais and Cranmer (2012).}
\item{target}{ The response network on which probabilities are based. Depending on whether the function is applied to an \code{ergm} or \code{btergm}/\code{mtergm} object, this can be either a single network or a list of networks. By default, the (list of) network(s) provided as the left-hand side of the (T)ERGM formula is used.}
\item{type}{ If \code{type = "tie"} is used, probabilities at the edge level are computed. For example, what is the probability of a specific node i to be connected to a specific node j given the rest of the network and given the model? If \code{type = "dyad"} is used, probabilities at the dyad level are computed. For example, what is the probability that node i is connected to node j but not vice-versa, or what is the probability that nodes i and j and mutually connected in a directed network? If \code{type = "node"} is used, probabilities at the node level are computed. For example, what is the probability that node i is connected to a set of three other j nodes given the rest of the network and the model? }
\item{i}{ A single (sender) node i or a set of (sender) nodes i. If \code{type = "node"} is used, this can be more than one node and should be provided as a vector. The \code{i} argument can be either provided as the index of the node in the sociomatrix (e.g., the fourth node would be \code{i = 4}) or the row name of the node in the sociomatrix (e.g., \code{i = "Peter"}). If more than one node is provided and \code{type = "node"}, there can be only one (receiver) node j. The \code{i} and \code{j} arguments are used to specify for which nodes probabilities should be computed. For example, what is the probability that \code{i = 4} is connected to \code{i = 7}? }
\item{j}{ A single (receiver) node j or a set of (receiver) nodes j. If \code{type = "node"} is used, this can be more than one node and should be provided as a vector. The \code{j} argument can be either provided as the index of the node in the sociomatrix (e.g., the fourth node would be \code{j = 4}) or the row name of the node in the sociomatrix (e.g., \code{j = "Mary"}). If more than one node is provided and \code{type = "node"}, there can be only one (sender) node i. The \code{i} and \code{j} arguments are used to specify for which nodes probabilities should be computed. For example, what is the probability that \code{i = 4} is connected to \code{i = 7}? }
\item{t}{ A vector of (numerical) time steps for which the probabilities should be computed. This only applies to \code{btergm} objects because \code{ergm} objects are by definition based on a single time step. By default, all available time steps are used. It is, for example, possible to compute probabilities only for a single time step by specifying, e.g., \code{t = 5} in order to compute probabilities for the fifth response network. }
\item{...}{ Further arguments to be handed over to subroutines. }
}
\references{
Desmarais, Bruce A. and Skyler J. Cranmer (2012):
Micro-Level Interpretation of Exponential Random Graph Models with
Application to Estuary Networks.
\emph{The Policy Studies Journal} 40(3): 402--434.
}
\seealso{
\link{btergm-package} \link{btergm} \link{timesteps.btergm}
}
\examples{
##### The following example is a TERGM adaptation of the #####
##### dyad-level example provided in figure 5(c) on page #####
##### 424 of Desmarais and Cranmer (2012) in the PSJ. At #####
##### each time step, it compares dyadic probabilities   #####
##### (no tie, unidirectional tie, and reciprocal tie    #####
##### probability) between a fitted model and a model    #####
##### where the reciprocity effect is fixed at 0 based   #####
##### on 20 randomly selected dyads per time step. The   #####
##### results are visualized using a grouped bar plot.   #####

\dontrun{
# create toy dataset and fit a model
networks <- list()
for (i in 1:3) {           # create 3 random networks with 10 actors
  mat <- matrix(rbinom(100, 1, 0.25), nrow = 10, ncol = 10)
  diag(mat) <- 0           # loops are excluded
  nw <- network(mat)       # create network object
  networks[[i]] <- nw      # add network to the list
}
fit <- btergm(networks ~ edges + istar(2) + mutual, R = 200)

# extract coefficients and create null hypothesis vector
null <- coef(fit)  # estimated coefs
null[3] <- 0       # set mutual term = 0

# sample 20 dyads per time step and compute probability ratios
probabilities <- matrix(nrow = 9, ncol = length(networks))
# nrow = 9 because three probabilities + upper and lower CIs
colnames(probabilities) <- paste("t =", 1:length(networks))
for (t in 1:length(networks)) {
  d <- dim(as.matrix(networks[[t]]))  # how many row and column nodes?
  size <- d[1] * d[2]                 # size of the matrix
  nw <- matrix(1:size, nrow = d[1], ncol = d[2])
  nw <- nw[lower.tri(nw)]             # sample only from lower triangle b/c
  samp <- sample(nw, 20)              # dyadic probabilities are symmetric
  prob.est.00 <- numeric(0)
  prob.est.01 <- numeric(0)
  prob.est.11 <- numeric(0)
  prob.null.00 <- numeric(0)
  prob.null.01 <- numeric(0)
  prob.null.11 <- numeric(0)
  for (k in 1:20) {
    i <- arrayInd(samp[k], d)[1, 1]   # recover 'i's and 'j's from sample
    j <- arrayInd(samp[k], d)[1, 2]
    # run interpretation function with estimated coefs and mutual = 0:
    int.est <- interpret(fit, type = "dyad", i = i, j = j, t = t)
    int.null <- interpret(fit, coefficients = null, type = "dyad", 
        i = i, j = j, t = t)
    prob.est.00 <- c(prob.est.00, int.est[[1]][1, 1])
    prob.est.11 <- c(prob.est.11, int.est[[1]][2, 2])
    mean.est.01 <- (int.est[[1]][1, 2] + int.est[[1]][2, 1]) / 2
    prob.est.01 <- c(prob.est.01, mean.est.01)
    prob.null.00 <- c(prob.null.00, int.null[[1]][1, 1])
    prob.null.11 <- c(prob.null.11, int.null[[1]][2, 2])
    mean.null.01 <- (int.null[[1]][1, 2] + int.null[[1]][2, 1]) / 2
    prob.null.01 <- c(prob.null.01, mean.null.01)
  }
  prob.ratio.00 <- prob.est.00 / prob.null.00  # ratio of est. and null hyp
  prob.ratio.01 <- prob.est.01 / prob.null.01
  prob.ratio.11 <- prob.est.11 / prob.null.11
  probabilities[1, t] <- mean(prob.ratio.00)   # mean estimated 00 tie prob
  probabilities[2, t] <- mean(prob.ratio.01)   # mean estimated 01 tie prob
  probabilities[3, t] <- mean(prob.ratio.11)   # mean estimated 11 tie prob
  ci.00 <- t.test(prob.ratio.00, conf.level = 0.99)$conf.int
  ci.01 <- t.test(prob.ratio.01, conf.level = 0.99)$conf.int
  ci.11 <- t.test(prob.ratio.11, conf.level = 0.99)$conf.int
  probabilities[4, t] <- ci.00[1]              # lower 00 conf. interval
  probabilities[5, t] <- ci.01[1]              # lower 01 conf. interval
  probabilities[6, t] <- ci.11[1]              # lower 11 conf. interval
  probabilities[7, t] <- ci.00[2]              # upper 00 conf. interval
  probabilities[8, t] <- ci.01[2]              # upper 01 conf. interval
  probabilities[9, t] <- ci.11[2]              # upper 11 conf. interval
}

# create barplots from probability ratios and CIs
require("gplots")
bp <- barplot2(probabilities[1:3, ], beside = TRUE, plot.ci = TRUE, 
    ci.l = probabilities[4:6, ], ci.u = probabilities[7:9, ], 
    col = c("tan", "tan2", "tan3"), ci.col = "grey40", 
    xlab = "Dyadic tie values", ylab = "Estimated Prob./Null Prob.")
mtext(1, at = bp, text = c("(0,0)", "(0,1)", "(1,1)"), line = 0, cex = 0.5)
}
}
