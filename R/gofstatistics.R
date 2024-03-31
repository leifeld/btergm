#' Statistics for goodness-of-fit assessment of network models
#' 
#' Statistics for goodness-of-fit assessment of network models.
#' 
#' These functions can be plugged into the \code{statistics} argument of the
#' \code{gof} methods in order to compare observed with simulated networks (see
#' the \link{gof-methods} help page). There are three types of statistics:
#' \enumerate{
#'   \item Univariate statistics, which aggregate a network into a single
#'     quantity. For example, modularity measures or density. The distribution
#'     of statistics can be displayed using histograms, density plots, and
#'     median bars. Univariate statistics take a sparse matrix (\code{mat})
#'     as an argument and return a single numeric value that summarize a network
#'     matrix.
#'   \item Multivariate statistics, which aggregate a network into a vector of
#'     quantities. For example, the distribution of geodesic distances, edgewise
#'     shared partners, or indegree. These statistics typically have multiple
#'     values, e.g., esp(1), esp(2), esp(3) etc. The results can be displayed
#'     using multiple boxplots for simulated networks and a black curve for the
#'     observed network(s). Multivariate statistics take a sparse matrix
#'     (\code{mat}) as an argument and return a vector of numeric values that
#'     summarize a network matrix.
#'   \item Tie prediction statistics, which predict dyad states the observed
#'     network(s) by the dyad states in the simulated networks. For example,
#'     receiver operating characteristics (ROC) or precision-recall curves (PR)
#'     of simulated networks based on the model, or ROC or PR predictions of
#'     community co-membership matrices of the simulated vs. the observed
#'     network(s). Tie prediction statistics take a list of simulated sparse
#'     network matrices and another list of observed sparse network matrices
#'     (possibly containing only a single sparse matrix) as arguments and return
#'     a \code{rocpr}, \code{roc}, or \code{pr} object (as created by the
#'     \link{rocpr} function).
#' }
#' 
#' Users can create their own statistics for use with the \code{gof} methods. To
#' do so, one needs to write a function that accepts and returns the respective
#' objects described in the enumeration above. It is advisable to look at the
#' definitions of some of the existing functions to add custom functions. It is
#' also possible to add an attribute called \code{label} to the return object,
#' which describes what is being returned by the function. This label will be
#' used as a descriptive label in the plot and for verbose output during
#' computations. The examples section contains an example of a custom user
#' statistic. Note that all statistics \emph{must} contain the \code{...}
#' argument to ensure that custom arguments of other statistics do not cause an
#' error.
#' 
#' To aid the development of custom statistics, the helper function
#' \code{comemb} is available: it accepts a vector of community memberships and
#' converts it to a co-membership matrix. This function is also used internally
#' by statistics like \code{walktrap.roc} and others.
#' 
#' @param vec A vector of community memberships in order to create a community
#'   co-membership matrix.
#' @param mat A sparse network matrix as created by the \code{Matrix} function
#'   in the \pkg{Matrix} package.
#' @param sim A list of simulated networks. Each element in the list should be a
#'   sparse matrix as created by the \code{\link[Matrix]{Matrix}} function in
#'   the \pkg{Matrix} package.
#' @param obs A list of observed (= target) networks. Each element in the list
#'   should be a sparse matrix as created by the \code{\link[Matrix]{Matrix}}
#'   function in the \pkg{Matrix} package.
#' @param roc Compute receiver-operating characteristics (ROC)?
#' @param pr Compute precision-recall curve (PR)?
#' @param joint Merge all time steps into a single big prediction task and
#'   compute predictive fit (instead of computing GOF for all time steps
#'   separately)?
#' @param pr.impute In some cases, the first precision value of the
#'   precision-recall curve is undefined. The \code{pr.impute} argument serves
#'   to impute this missing value to ensure that the AUC-PR value is not
#'   severely biased. Possible values are \code{"no"} for no imputation,
#'   \code{"one"} for using a value of \code{1.0}, \code{"second"} for using the
#'   next (= adjacent) precision value, \code{"poly1"} for fitting a straight
#'   line through the remaining curve to predict the first value, \code{"poly2"}
#'   for fitting a second-order polynomial curve etc. until \code{"poly9"}.
#'   Warning: this is a pragmatic solution. Please double-check whether the
#'   imputation makes sense. This can be checked by plotting the resulting
#'   object and using the \code{pr.poly} argument to plot the predicted curve on
#'   top of the actual PR curve.
#' @param ... Additional arguments. This must be present in all auxiliary GOF
#'   statistics.
#'
#' @references
#' Leifeld, Philip, Skyler J. Cranmer and Bruce A. Desmarais (2018): Temporal
#' Exponential Random Graph Models with btergm: Estimation and Bootstrap
#' Confidence Intervals. \emph{Journal of Statistical Software} 83(6): 1--36.
#' \doi{10.18637/jss.v083.i06}.
#' 
#' @examples
#' # To see how these statistics are used, look at the examples section of 
#' # ?"gof-methods". The following example illustrates how custom 
#' # statistics can be created. Suppose one is interested in the density 
#' # of a network. Then a univariate statistic can be created as follows.
#' 
#' dens <- function(mat, ...) {        # univariate: one argument
#'   mat <- as.matrix(mat)             # sparse matrix -> normal matrix
#'   d <- sna::gden(mat)               # compute the actual statistic
#'   attributes(d)$label <- "Density"  # add a descriptive label
#'   return(d)                         # return the statistic
#' }
#' 
#' # Note that the '...' argument must be present in all statistics. 
#' # Now the statistic can be used in the statistics argument of one of 
#' # the gof methods.
#' 
#' # For illustrative purposes, let us consider an existing statistic, the 
#' # indegree distribution, a multivariate statistic. It also accepts a 
#' # single argument. Note that the sparse matrix is converted to a 
#' # normal matrix object when it is used. First, statnet's summary 
#' # method is used to compute the statistic. Names are attached to the 
#' # resulting vector for the different indegree values. Then the vector 
#' # is returned.
#' 
#' ideg <- function(mat, ...) {
#'   d <- summary(mat ~ idegree(0:(nrow(mat) - 1)))
#'   names(d) <- 0:(length(d) - 1)
#'   attributes(d)$label <- "Indegree"
#'   return(d)
#' }
#' 
#' # See the gofstatistics.R file in the package for more complex examples.
#' 
#' @name gof-statistics
#' @aliases gof-statistics gofstatistics
#' 
#' @importFrom Matrix Matrix
#' @importFrom ROCR performance prediction
#' @importFrom igraph graph.adjacency spinglass.community fastgreedy.community
#'   cluster_louvain modularity edge.betweenness.community optimal.community
#'   walktrap.community
#' @importFrom sna triad.census gden symmetrize
#' @importFrom network summary.network
#' @importFrom ergm ergm.geodistdist
NULL

#' @describeIn gof-statistics Multivariate GOF statistic: dyad-wise shared
#'   partner distribution
#' @export
dsp <- function(mat, ...) {
  d <- summary(mat ~ dsp(0:(nrow(mat) - 2)))
  names(d) <- 0:(length(d) - 1)
  attributes(d)$label <- "Dyad-wise shared partners"
  return(d)
}

#' @describeIn gof-statistics Multivariate GOF statistic: edge-wise shared
#'   partner distribution
#' @export
esp <- function(mat, ...) {
  d <- summary(mat ~ esp(0:(nrow(mat) - 2)))
  names(d) <- 0:(length(d) - 1)
  attributes(d)$label <- "Edge-wise shared partners"
  return(d)
}

#' @describeIn gof-statistics Multivariate GOF statistic: non-edge-wise shared
#'   partner distribution
#' @export
nsp <- function(mat, ...) {
  d <- summary(mat ~ nsp(0:(nrow(mat) - 2)))
  names(d) <- 0:(length(d) - 1)
  attributes(d)$label <- "Non-edge-wise shared partners"
  return(d)
}

#' @describeIn gof-statistics Multivariate GOF statistic: degree distribution
#' @export
deg <- function(mat, ...) {
  d <- summary(network(as.matrix(mat), directed = FALSE) ~ degree(0:(nrow(mat) 
      - 1)))
  names(d) <- 0:(length(d) - 1)
  attributes(d)$label <- "Degree"
  return(d)
}

#' @describeIn gof-statistics Multivariate GOF statistic: degree distribution
#'   for the first mode
#' @export
b1deg <- function(mat, ...) {
  d <- summary(network(as.matrix(mat), directed = FALSE, bipartite = TRUE) ~ 
      b1degree(0:nrow(mat)))
  names(d) <- 0:(length(d)- 1)
  attributes(d)$label <- "Degree (first mode)"
  return(d)
}

#' @describeIn gof-statistics Multivariate GOF statistic: degree distribution
#'   for the second mode
#' @export
b2deg <- function(mat, ...) {
  d <- summary(network(as.matrix(mat), directed = FALSE, bipartite = TRUE) ~ 
      b2degree(0:ncol(mat)))
  names(d) <- 0:(length(d)- 1)
  attributes(d)$label <- "Degree (second mode)"
  return(d)
}

#' @describeIn gof-statistics Multivariate GOF statistic: outdegree distribution
#' @export
odeg <- function(mat, ...) {
  d <- summary(network(as.matrix(mat), directed = TRUE) ~ odegree(0:(nrow(mat) 
      - 1)))
  names(d) <- 0:(length(d) - 1)
  attributes(d)$label <- "Outdegree"
  return(d)
}

#' @describeIn gof-statistics Multivariate GOF statistic: indegree distribution
#' @export
ideg <- function(mat, ...) {
  d <- summary(network(as.matrix(mat), directed = TRUE) ~ idegree(0:(nrow(mat) 
      - 1)))
  names(d) <- 0:(length(d) - 1)
  attributes(d)$label <- "Indegree"
  return(d)
}

#' @describeIn gof-statistics Multivariate GOF statistic: k-star distribution
#' @export
kstar <- function(mat, ...) {
  d <- summary(network(as.matrix(mat), directed = FALSE) ~ kstar(0:(nrow(mat) 
      - 1)))
  names(d) <- 0:(length(d) - 1)
  attributes(d)$label <- "k-star"
  return(d)
}

#' @describeIn gof-statistics Multivariate GOF statistic: k-star distribution
#'   for the first mode
#' @export
b1star <- function(mat, ...) {
  d <- summary(network(as.matrix(mat), directed = FALSE, bipartite = TRUE) ~ 
      b1star(0:nrow(mat)))
  names(d) <- 0:(length(d) - 1)
  attributes(d)$label <- "k-star (first mode)"
  return(d)
}

#' @describeIn gof-statistics Multivariate GOF statistic: k-star distribution
#'   for the second mode
#' @export
b2star <- function(mat, ...) {
  d <- summary(network(as.matrix(mat), directed = FALSE, bipartite = TRUE) ~ 
      b2star(0:nrow(mat)))
  names(d) <- 0:(length(d) - 1)
  attributes(d)$label <- "k-star (second mode)"
  return(d)
}

#' @describeIn gof-statistics Multivariate GOF statistic: outgoing k-star
#'   distribution
#' @export
ostar <- function(mat, ...) {
  d <- summary(network(as.matrix(mat), directed = TRUE) ~ ostar(0:(nrow(mat) 
      - 1)))
  names(d) <- 0:(length(d) - 1)
  attributes(d)$label <- "Outgoing k-star"
  return(d)
}

#' @describeIn gof-statistics Multivariate GOF statistic: incoming k-star
#'   distribution
#' @export
istar <- function(mat, ...) {
  d <- summary(network(as.matrix(mat), directed = TRUE) ~ istar(0:(nrow(mat) 
      - 1)))
  names(d) <- 0:(length(d) - 1)
  attributes(d)$label <- "Incoming k-star"
  return(d)
}

#' @describeIn gof-statistics Multivariate GOF statistic: k-cycle distribution
#' @export
kcycle <- function(mat, ...) {
  d <- summary(mat ~ cycle(0:(nrow(mat) - 1)))
  names(d) <- 0:(length(d) - 1)
  attributes(d)$label <- "Cycle"
  return(d)
}

#' @describeIn gof-statistics Multivariate GOF statistic: geodesic distance
#'   distribution
#' @export
geodesic <- function(mat, ...) {
  mat[is.na(mat)] <- 0
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
  g <- fillup(ergm::ergm.geodistdist(network(as.matrix(mat), directed = TRUE)), 
      nrow(mat) - 1)
  attributes(g)$label <- "Geodesic distances"
  return(g)
}

#' @describeIn gof-statistics Multivariate GOF statistic: triad census in
#'   directed networks
#'
#' @importFrom network network
#' @export
triad.directed <- function(mat, ...) {
  tr <- sna::triad.census(network::network(as.matrix(mat), directed = TRUE), 
      mode = "digraph")[1, ]
  attributes(tr)$label <- "Triad census"
  return(tr)
}

#' @describeIn gof-statistics Multivariate GOF statistic: triad census in
#'   undirected networks
#'
#' @importFrom network network
#' @export
triad.undirected <- function(mat, ...) {
  tr <- sna::triad.census(network::network(as.matrix(mat), directed = FALSE), 
      mode = "graph")[1, ]
  attributes(tr)$label <- "Triad census"
  return(tr)
}

#' @describeIn gof-statistics Helper function: create community co-membership
#'   matrix
#' @export
comemb <- function(vec) {
  comemb <- matrix(0, nrow = length(vec), ncol = length(vec))
  for (a in 1:length(vec)) {
    for (b in 1:length(vec)) {
      if (vec[a] == vec[b]) {
        comemb[a, b] <- 1
      } else {
        comemb[a, b] <- 0
      }
    }
  }
  return(comemb)
}

#' @describeIn gof-statistics Univariate GOF statistic: Walktrap modularity
#'   distribution
#' @export
walktrap.modularity <- function(mat, ...) {
  mat[is.na(mat)] <- 0
  if (is.mat.directed(as.matrix(mat))) {
    m <- "directed"
  } else {
    m <- "undirected"
  }
  if (sum(mat) == 0) {
    mod <- 0
  } else {
    g <- igraph::graph.adjacency(as.matrix(mat), mode = m)
    wt <- igraph::walktrap.community(g)
    mod <- igraph::modularity(wt)
  }
  attributes(mod)$label <- "Modularity (walktrap)"
  return(mod)
}

#' @describeIn gof-statistics Tie prediction GOF statistic: ROC of Walktrap
#'   community detection. Receiver-operating characteristics of predicting the
#'   community structure in the observed network(s) by the community structure
#'   in the simulated networks, as computed by the Walktrap algorithm.
#' @export
walktrap.roc <- function(sim, obs, ...) {
  for (i in 1:length(sim)) {
    sim[[i]][is.na(sim[[i]])] <- 0
  }
  for (i in 1:length(obs)) {
    obs[[i]][is.na(obs[[i]])] <- 0
  }
  fun <- function(x) {
    m <- is.mat.directed(as.matrix(x))
    if (m == TRUE) {
      m <- "directed"
    } else {
      m <- "undirected"
    }
    if (sum(x) == 0) {
      memb <- rep(1, nrow(x))
    } else {
      g <- igraph::graph.adjacency(as.matrix(x), mode = m)
      memb <- igraph::walktrap.community(g)$membership
    }
    return(comemb(memb))
  }
  sim <- lapply(sim, fun)
  obs <- lapply(obs, fun)
  object <- suppressMessages(rocpr(sim, obs)[-c(4, 5, 8, 9)])
  class(object) <- "roc"
  object$label <- "Walktrap community comembership prediction"
  attributes(object)$label <- object$label
  return(object)
}

#' @describeIn gof-statistics Tie prediction GOF statistic: PR of Walktrap
#'   community detection. Precision-recall curve for predicting the community
#'   structure in the observed network(s) by the community structure in the
#'   simulated networks, as computed by the Walktrap algorithm.
#' @export
walktrap.pr <- function(sim, obs, ...) {
  for (i in 1:length(sim)) {
    sim[[i]][is.na(sim[[i]])] <- 0
  }
  for (i in 1:length(obs)) {
    obs[[i]][is.na(obs[[i]])] <- 0
  }
  fun <- function(x) {
    m <- is.mat.directed(as.matrix(x))
    if (m == TRUE) {
      m <- "directed"
    } else {
      m <- "undirected"
    }
    if (sum(x) == 0) {
      memb <- rep(1, nrow(x))
    } else {
      g.obs <- igraph::graph.adjacency(as.matrix(x), mode = m)
      memb <- igraph::walktrap.community(g.obs)$membership
    }
    return(comemb(memb))
  }
  sim <- lapply(sim, fun)
  obs <- lapply(obs, fun)
  object <- suppressMessages(rocpr(sim, obs)[-c(2, 3, 6, 7)])
  class(object) <- "pr"
  object$label <- "Walktrap community comembership prediction"
  attributes(object)$label <- object$label
  return(object)
}

#' @describeIn gof-statistics Univariate GOF statistic: fast and greedy
#'   modularity distribution
#' @export
fastgreedy.modularity <- function(mat, ...) {
  mat[is.na(mat)] <- 0
  if (sum(mat) == 0) {
    mod <- 0
  } else {
    g <- igraph::graph.adjacency(as.matrix(mat), mode = "undirected")
    wt <- igraph::fastgreedy.community(g)
    mod <- igraph::modularity(wt)
  }
  attributes(mod)$label <- "Modularity (fast & greedy)"
  return(mod)
}

#' @describeIn gof-statistics Tie prediction GOF statistic: ROC of fast and
#'   greedy community detection. Receiver-operating characteristics of
#'   predicting the community structure in the observed network(s) by the
#'   community structure in the simulated networks, as computed by the fast and
#'   greedy algorithm. Only sensible with undirected networks.
#' @export
fastgreedy.roc <- function(sim, obs, ...) {
  for (i in 1:length(sim)) {
    sim[[i]][is.na(sim[[i]])] <- 0
  }
  for (i in 1:length(obs)) {
    obs[[i]][is.na(obs[[i]])] <- 0
  }
  fun <- function(x) {
    if (sum(x) == 0) {
      memb <- rep(1, nrow(x))
    } else {
      g <- igraph::graph.adjacency(as.matrix(x), mode = "undirected")
      memb <- igraph::fastgreedy.community(g)$membership
    }
    return(comemb(memb))
  }
  sim <- lapply(sim, fun)
  obs <- lapply(obs, fun)
  object <- suppressMessages(rocpr(sim, obs)[-c(4, 5, 8, 9)])
  class(object) <- "roc"
  object$label <- "Fast & greedy community comembership prediction"
  attributes(object)$label <- object$label
  return(object)
}

#' @describeIn gof-statistics Tie prediction GOF statistic: PR of fast and
#'   greedy community detection. Precision-recall curve for predicting the
#'   community structure in the observed network(s) by the community structure
#'   in the simulated networks, as computed by the fast and greedy algorithm.
#'   Only sensible with undirected networks.
#' @export
fastgreedy.pr <- function(sim, obs, ...) {
  for (i in 1:length(sim)) {
    sim[[i]][is.na(sim[[i]])] <- 0
  }
  for (i in 1:length(obs)) {
    obs[[i]][is.na(obs[[i]])] <- 0
  }
  fun <- function(x) {
    if (sum(x) == 0) {
      memb <- rep(1, nrow(x))
    } else {
      g <- igraph::graph.adjacency(as.matrix(x), mode = "undirected")
      memb <- igraph::fastgreedy.community(g)$membership
    }
    return(comemb(memb))
  }
  sim <- lapply(sim, fun)
  obs <- lapply(obs, fun)
  object <- suppressMessages(rocpr(sim, obs)[-c(2, 3, 6, 7)])
  class(object) <- "pr"
  object$label <- "Fast & greedy community comembership prediction"
  attributes(object)$label <- object$label
  return(object)
}

#' @describeIn gof-statistics Univariate GOF statistic: Louvain clustering
#'   modularity distribution
#' @export
louvain.modularity <- function(mat, ...) {
  if (is.mat.directed(as.matrix(mat))) {
    m <- "directed"
  } else {
    m <- "undirected"
  }
  mat[is.na(mat)] <- 0
  if (sum(mat) == 0) {
    mod <- 0
  } else {
    g <- igraph::graph.adjacency(as.matrix(mat), mode = m)
    clus <- igraph::cluster_louvain(g)
    mod <- igraph::modularity(clus)
  }
  attributes(mod)$label <- "Modularity (Louvain)"
  return(mod)
}

#' @describeIn gof-statistics Tie prediction GOF statistic: ROC of Louvain
#'   community detection. Receiver-operating characteristics of predicting the
#'   community structure in the observed network(s) by the community structure
#'   in the simulated networks, as computed by the Louvain algorithm.
#' @export
louvain.roc <- function(sim, obs, ...) {
  for (i in 1:length(sim)) {
    sim[[i]][is.na(sim[[i]])] <- 0
  }
  for (i in 1:length(obs)) {
    obs[[i]][is.na(obs[[i]])] <- 0
  }
  fun <- function(x) {
    m <- is.mat.directed(as.matrix(x))
    if (m == TRUE) {
      m <- "directed"
    } else {
      m <- "undirected"
    }
    if (sum(x) == 0) {
      memb <- rep(1, nrow(x))
    } else {
      g <- igraph::graph.adjacency(as.matrix(x), mode = m)
      memb <- igraph::cluster_louvain(g)$membership
    }
    return(comemb(memb))
  }
  sim <- lapply(sim, fun)
  obs <- lapply(obs, fun)
  object <- suppressMessages(rocpr(sim, obs)[-c(4, 5, 8, 9)])
  class(object) <- "roc"
  object$label <- "Louvain community comembership prediction"
  attributes(object)$label <- object$label
  return(object)
}

#' @describeIn gof-statistics Tie prediction GOF statistic: PR of Louvain
#'   community detection. Precision-recall curve for predicting the community
#'   structure in the observed network(s) by the community structure in the
#'   simulated networks, as computed by the Louvain algorithm.
#' @export
louvain.pr <- function(sim, obs, ...) {
  for (i in 1:length(sim)) {
    sim[[i]][is.na(sim[[i]])] <- 0
  }
  for (i in 1:length(obs)) {
    obs[[i]][is.na(obs[[i]])] <- 0
  }
  fun <- function(x) {
    m <- is.mat.directed(as.matrix(x))
    if (m == TRUE) {
      m <- "directed"
    } else {
      m <- "undirected"
    }
    if (sum(x) == 0) {
      memb <- rep(1, nrow(x))
    } else {
      g <- igraph::graph.adjacency(as.matrix(x), mode = m)
      memb <- igraph::cluster_louvain(g)$membership
    }
    return(comemb(memb))
  }
  sim <- lapply(sim, fun)
  obs <- lapply(obs, fun)
  object <- suppressMessages(rocpr(sim, obs)[-c(2, 3, 6, 7)])
  class(object) <- "pr"
  object$label <- "Louvain community comembership prediction"
  attributes(object)$label <- object$label
  return(object)
}

#' @describeIn gof-statistics Univariate GOF statistic: maximal modularity
#'   distribution
#' @export
maxmod.modularity <- function(mat, ...) {
  mat[is.na(mat)] <- 0
  if (is.mat.directed(as.matrix(mat))) {
    m <- "directed"
  } else {
    m <- "undirected"
  }
  if (sum(mat) == 0) {
    mod <- 0
  } else {
    g <- igraph::graph.adjacency(as.matrix(mat), mode = m)
    wt <- igraph::optimal.community(g)
    mod <- igraph::modularity(wt)
  }
  attributes(mod)$label <- "Maximum modularity"
  return(mod)
}

#' @describeIn gof-statistics Tie prediction GOF statistic: ROC of maximal
#'   modularity community detection. Receiver-operating characteristics of
#'   predicting the community structure in the observed network(s) by the
#'   community structure in the simulated networks, as computed by the
#'   modularity maximization algorithm.
#' @export
maxmod.roc <- function(sim, obs, ...) {
  for (i in 1:length(sim)) {
    sim[[i]][is.na(sim[[i]])] <- 0
  }
  for (i in 1:length(obs)) {
    obs[[i]][is.na(obs[[i]])] <- 0
  }
  fun <- function(x) {
    m <- is.mat.directed(as.matrix(x))
    if (m == TRUE) {
      m <- "directed"
    } else {
      m <- "undirected"
    }
    if (sum(x) == 0) {
      memb <- rep(1, nrow(x))
    } else {
      g <- igraph::graph.adjacency(as.matrix(x), mode = m)
      memb <- igraph::optimal.community(g)$membership
    }
    return(comemb(memb))
  }
  sim <- lapply(sim, fun)
  obs <- lapply(obs, fun)
  object <- suppressMessages(rocpr(sim, obs)[-c(4, 5, 8, 9)])
  class(object) <- "roc"
  object$label <- "Maximum modularity community comembership prediction"
  attributes(object)$label <- object$label
  return(object)
}

#' @describeIn gof-statistics Tie prediction GOF statistic: PR of maximal
#'   modularity community detection. Precision-recall curve for predicting the
#'   community structure in the observed network(s) by the community structure
#'   in the simulated networks, as computed by the modularity maximization
#'   algorithm.
#' @export
maxmod.pr <- function(sim, obs, ...) {
  for (i in 1:length(sim)) {
    sim[[i]][is.na(sim[[i]])] <- 0
  }
  for (i in 1:length(obs)) {
    obs[[i]][is.na(obs[[i]])] <- 0
  }
  fun <- function(x) {
    m <- is.mat.directed(as.matrix(x))
    if (m == TRUE) {
      m <- "directed"
    } else {
      m <- "undirected"
    }
    if (sum(x) == 0) {
      memb <- rep(1, nrow(x))
    } else {
      g <- igraph::graph.adjacency(as.matrix(x), mode = m)
      memb <- igraph::optimal.community(g)$membership
    }
    return(comemb(memb))
  }
  sim <- lapply(sim, fun)
  obs <- lapply(obs, fun)
  object <- suppressMessages(rocpr(sim, obs)[-c(2, 3, 6, 7)])
  class(object) <- "pr"
  object$label <- "Maximum modularity community comembership prediction"
  attributes(object)$label <- object$label
  return(object)
}

#' @describeIn gof-statistics Univariate GOF statistic: edge betweenness
#'   modularity distribution
#' @export
edgebetweenness.modularity <- function(mat, ...) {
  mat[is.na(mat)] <- 0
  if (is.mat.directed(as.matrix(mat))) {
    m <- "directed"
  } else {
    m <- "undirected"
  }
  if (sum(mat) == 0) {
    mod <- 0
  } else {
    g <- igraph::graph.adjacency(as.matrix(mat), mode = m)
    eb <- igraph::edge.betweenness.community(g)
    mod <- igraph::modularity(eb)
  }
  attributes(mod)$label <- "Modularity (edge betweenness)"
  return(mod)
}

#' @describeIn gof-statistics Tie prediction GOF statistic: ROC of edge
#'   betweenness community detection. Receiver-operating characteristics of
#'   predicting the community structure in the observed network(s) by the
#'   community structure in the simulated networks, as computed by the
#'   Girvan-Newman edge betweenness community detection method.
#' @export
edgebetweenness.roc <- function(sim, obs, ...) {
  for (i in 1:length(sim)) {
    sim[[i]][is.na(sim[[i]])] <- 0
  }
  for (i in 1:length(obs)) {
    obs[[i]][is.na(obs[[i]])] <- 0
  }
  fun <- function(x) {
    m <- is.mat.directed(as.matrix(x))
    if (m == TRUE) {
      m <- "directed"
    } else {
      m <- "undirected"
    }
    if (sum(x) == 0) {
      memb <- rep(1, nrow(x))
    } else {
      g <- igraph::graph.adjacency(as.matrix(x), mode = m)
      memb <- igraph::edge.betweenness.community(g)$membership
    }
    return(comemb(memb))
  }
  sim <- lapply(sim, fun)
  obs <- lapply(obs, fun)
  object <- suppressMessages(rocpr(sim, obs)[-c(4, 5, 8, 9)])
  class(object) <- "roc"
  object$label <- "Edge betweenness community comembership prediction"
  attributes(object)$label <- object$label
  return(object)
}

#' @describeIn gof-statistics Tie prediction GOF statistic: PR of edge
#'   betweenness community detection. Precision-recall curve for predicting the
#'   community structure in the observed network(s) by the community structure
#'   in the simulated networks, as computed by the Girvan-Newman edge
#'   betweenness community detection method.
#' @export
edgebetweenness.pr <- function(sim, obs, ...) {
  for (i in 1:length(sim)) {
    sim[[i]][is.na(sim[[i]])] <- 0
  }
  for (i in 1:length(obs)) {
    obs[[i]][is.na(obs[[i]])] <- 0
  }
  fun <- function(x) {
    m <- is.mat.directed(as.matrix(x))
    if (m == TRUE) {
      m <- "directed"
    } else {
      m <- "undirected"
    }
    if (sum(x) == 0) {
      memb <- rep(1, nrow(x))
    } else {
      g <- igraph::graph.adjacency(as.matrix(x), mode = m)
      memb <- igraph::edge.betweenness.community(g)$membership
    }
    return(comemb(memb))
  }
  sim <- lapply(sim, fun)
  obs <- lapply(obs, fun)
  object <- suppressMessages(rocpr(sim, obs)[-c(2, 3, 6, 7)])
  class(object) <- "pr"
  object$label <- "Edge betweenness community comembership prediction"
  attributes(object)$label <- object$label
  return(object)
}

#' @describeIn gof-statistics Univariate GOF statistic: spinglass modularity
#'   distribution
#' @export
spinglass.modularity <- function(mat, ...) {
  mat[is.na(mat)] <- 0
  if (is.mat.directed(as.matrix(mat))) {
    m <- "directed"
  } else {
    m <- "undirected"
  }
  if (sum(mat) == 0) {
    mod <- 0
  } else {
    g <- igraph::graph.adjacency(as.matrix(mat), mode = m)
    eb <- igraph::spinglass.community(g)
    mod <- igraph::modularity(eb)
  }
  attributes(mod)$label <- "Modularity (spinglass)"
  return(mod)
}

#' @describeIn gof-statistics Tie prediction GOF statistic: ROC of spinglass
#'   community detection. Receiver-operating characteristics of predicting the
#'   community structure in the observed network(s) by the community structure
#'   in the simulated networks, as computed by the Spinglass algorithm.
#' @export
spinglass.roc <- function(sim, obs, ...) {
  for (i in 1:length(sim)) {
    sim[[i]][is.na(sim[[i]])] <- 0
  }
  for (i in 1:length(obs)) {
    obs[[i]][is.na(obs[[i]])] <- 0
  }
  fun <- function(x) {
    m <- is.mat.directed(as.matrix(x))
    if (m == TRUE) {
      m <- "directed"
    } else {
      m <- "undirected"
    }
    if (sum(x) == 0) {
      memb <- rep(1, nrow(x))
    } else {
      g <- igraph::graph.adjacency(as.matrix(x), mode = m)
      memb <- igraph::spinglass.community(g)$membership
    }
    return(comemb(memb))
  }
  sim <- lapply(sim, fun)
  obs <- lapply(obs, fun)
  object <- suppressMessages(rocpr(sim, obs)[-c(4, 5, 8, 9)])
  class(object) <- "roc"
  object$label <- "Spinglass community comembership prediction"
  attributes(object)$label <- object$label
  return(object)
}

#' @describeIn gof-statistics Tie prediction GOF statistic: PR of spinglass
#'   community detection. Precision-recall curve for predicting the community
#'   structure in the observed network(s) by the community structure in the
#'   simulated networks, as computed by the Spinglass algorithm.
#' @export
spinglass.pr <- function(sim, obs, ...) {
  for (i in 1:length(sim)) {
    sim[[i]][is.na(sim[[i]])] <- 0
  }
  for (i in 1:length(obs)) {
    obs[[i]][is.na(obs[[i]])] <- 0
  }
  fun <- function(x) {
    m <- is.mat.directed(as.matrix(x))
    if (m == TRUE) {
      m <- "directed"
    } else {
      m <- "undirected"
    }
    if (sum(x) == 0) {
      memb <- rep(1, nrow(x))
    } else {
      g <- igraph::graph.adjacency(as.matrix(x), mode = m)
      memb <- igraph::spinglass.community(g)$membership
    }
    return(comemb(memb))
  }
  sim <- lapply(sim, fun)
  obs <- lapply(obs, fun)
  object <- suppressMessages(rocpr(sim, obs)[-c(2, 3, 6, 7)])
  class(object) <- "pr"
  object$label <- "Spinglass community comembership prediction"
  attributes(object)$label <- object$label
  return(object)
}


#' AUC-PR function -- modified version of: https://github.com/ipa-tys/ROCR/pull/2
#' @noRd
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

#' @describeIn gof-statistics Tie prediction GOF statistic: ROC and PR curves.
#'   Receiver-operating characteristics (ROC) and precision-recall curve (PR).
#'   Prediction of the dyad states of the observed network(s) by the dyad states
#'   of the simulated networks.
#' @importFrom network network
#' @export
rocpr <- function(sim, obs, roc = TRUE, pr = TRUE, joint = TRUE, 
    pr.impute = "poly4", ...) {
  
  directed <- sapply(obs, is.mat.directed)
  twomode <- !sapply(obs, is.mat.onemode)
  
  # create random graphs with corresponding tie probability of each time step
  nsim <- length(sim) / length(obs)
  rg <- list()
  for (i in 1:length(obs)) {
    rn <- nrow(obs[[i]])
    cn <- ncol(obs[[i]])
    n <- rn * cn
    dens <- sna::gden(network::network(as.matrix(obs[[i]]), directed = directed[i], 
        bipartite = twomode[i]))
    rlist <- list()
    for (j in 1:nsim) {
      r <- matrix(rbinom(n = n, size = 1, prob = dens), nrow = rn, ncol = cn)
      if (twomode[[i]] == FALSE) {
        diag(r) <- 0
        if (directed[[i]] == FALSE) {
          r <- symmetrize(r, rule = "upper")
        }
      }
      rg[length(rg) + 1] <- list(Matrix::Matrix(r))
    }
  }
  
  # ROCR
  target.pr <- list()
  rgraph.pr <- list()
  target.y <- list()
  for (j in 1:length(obs)) {
    net <- obs[[j]]
    index.start <- j * nsim - nsim + 1 # sim 1-100 for obs 1, 101-200 for 2 etc.
    index.stop <- (j + 1) * nsim - nsim
    sums <- 0
    rg.sums <- 0
    for (i in index.start:index.stop) {
      sums <- sums + sim[[i]]
      rg.sums <- rg.sums + rg[[i]]
    }
    sums <- sums / length(sim)
    rg.sums <- rg.sums / length(rg)
    if (directed[[j]] == TRUE || twomode[[j]] == TRUE) {
      predicted <- c(as.matrix(sums))
      y <- c(as.matrix(net))
      rg.pr <- c(as.matrix(rg.sums))
    } else {
      predicted <- sums[lower.tri(sums)]
      y <- as.matrix(net)[lower.tri(as.matrix(net))]
      rg.pr <- rg.sums[lower.tri(rg.sums)]
    }
    predicted <- predicted[!is.na(y)]
    rg.pr <- rg.pr[!is.na(y)]
    y <- y[!is.na(y)]
    target.pr[[j]] <- predicted
    rgraph.pr[[j]] <- rg.pr
    target.y[[j]] <- y
  }
  if (joint == TRUE) {  # merge into one single prediction task
    target.pr <- unlist(target.pr)
    target.y <- unlist(target.y)
    rgraph.pr <- unlist(rgraph.pr)
  }
  pred <- prediction(target.pr, target.y)
  # print(unlist(performance(pred, measure = "auc")@y.values))
  # print(unlist(performance(pred, measure = "aucpr")@y.values))
  rocperf <- performance(pred, "tpr", "fpr")  # ROC curve
  prperf <- performance(pred, "ppv", "tpr")  # precision-recall curve
  rg.pred <- prediction(rgraph.pr, target.y)
  rg.roc <- performance(rg.pred, "tpr", "fpr")  # ROC curve
  rg.pr <- performance(rg.pred, "ppv", "tpr")  # precision-recall curve
  
  # impute the first PR value (which is sometimes NaN)
  for (j in 1:length(prperf@y.values)) {
    fp <- pred@fp[[j]]
    tp <- pred@tp[[j]]
    if (fp[1] == 0 & tp[1] == 0) {
      prperf@y.values[[j]][1] <- 1
    } else if (pr.impute == "no") {
      message(paste0("t = ", j, ": warning -- the first PR value was not ", 
        "imputed; this may lead to underestimated AUC-PR values."))
      # do nothing
    } else if (is.nan(prperf@y.values[[j]][1])) {
      if (pr.impute == "second") {
        message(paste0("t = ", j, ": imputing the first PR value by the ", 
            "next (= adjacent) value."))
        prperf@y.values[[j]][1] <- prperf@y.values[[j]][2]
      } else if (pr.impute == "one") {
        message(paste0("t = ", j, ": imputing the first PR value by the ", 
            "maximum value of 1."))
        prperf@y.values[[j]][1] <- 1
      } else if (grepl("^poly[1-9]", pr.impute)) {
        num <- as.numeric(substr(pr.impute, 5, 5))
        message(paste0("t = ", j, ": imputing the first PR value using a ", 
            "polynomial of order ", num, ". Check the results by plotting ",
            "the GOF object using the \"pr.poly = ", num, "\" argument."))
        p <- data.frame(poly(prperf@x.values[[j]], num, raw = TRUE))
        fit <- lm(prperf@y.values[[j]] ~ ., data = p)
        prperf@y.values[[j]][1] <- predict(fit, newdata = p[1, ])
      } else {
        message(paste0("t = ", j, ": PR imputation method not recognized. ", 
            "Not using any imputation."))
      }
      if (prperf@y.values[[j]][1] < 0) {
        prperf@y.values[[j]][0] <- 0
      }
      if (prperf@y.values[[j]][1] > 1) {
        prperf@y.values[[j]][1] <- 1
      }
    }
  }
  
  # impute the first PR value of random graph
  for (j in 1:length(rg.pr@y.values)) {
    fp <- rg.pred@fp[[j]]
    tp <- rg.pred@tp[[j]]
    if (fp[1] == 0 & tp[1] == 0) {
      rg.pr@y.values[[j]][1] <- 1
    } else if (pr.impute == "no") {
      message(paste0("t = ", j, ": warning -- the first PR value was not ", 
        "imputed; this may lead to underestimated AUC-PR values for the ", 
        "random graph."))
      # do nothing
    } else if (is.nan(rg.pr@y.values[[j]][1])) {
      if (pr.impute == "second") {
        message(paste0("t = ", j, ": imputing the first PR value by the ", 
            "next (= adjacent) value for the random graph."))
        rg.pr@y.values[[j]][1] <- rg.pr@y.values[[j]][2]
      } else if (pr.impute == "one") {
        message(paste0("t = ", j, ": imputing the first PR value by the ", 
            "maxmum value of 1 for the random graph."))
        rg.pr@y.values[[j]][1] <- 1
      } else if (grepl("^poly[1-9]", pr.impute)) {
        num <- as.numeric(substr(pr.impute, 5, 5))
        message(paste0("t = ", j, ": imputing the first PR value using a ", 
            "polynomial of order ", num, " for the random graph. Check the ", 
            "results by plotting the GOF object using the \"pr.poly = ", num, 
            "\" argument."))
        p <- data.frame(poly(rg.pr@x.values[[j]], num, raw = TRUE))
        fit <- lm(rg.pr@y.values[[j]] ~ ., data = p)
        rg.pr@y.values[[j]][1] <- predict(fit, newdata = p[1, ])
      } else {
        message(paste0("t = ", j, ": PR imputation method not recognized. ", 
            "Not using any imputation."))
      }
      if (rg.pr@y.values[[j]][1] < 0) {
        rg.pr@y.values[[j]][0] <- 0
      }
      if (rg.pr@y.values[[j]][1] > 1) {
        rg.pr@y.values[[j]][1] <- 1
      }
    }
  }
  
  auc.roc <- unlist(performance(pred, measure = "auc")@y.values)  # ROC-AUC
  auc.pr <- aucpr(pred, precision = prperf@y.values, recall = prperf@x.values)  # PR-AUC
  rg.auc.roc <- unlist(performance(rg.pred, measure = "auc")@y.values)
  rg.auc.pr <- aucpr(rg.pred, precision = rg.pr@y.values, 
      recall = rg.pr@x.values)
  
  rocpr <- list()
  rocpr$type <- "rocpr"
  rocpr$auc.roc <- auc.roc
  rocpr$auc.roc.rgraph <- rg.auc.roc
  rocpr$auc.pr <- auc.pr
  rocpr$auc.pr.rgraph <- rg.auc.pr
  rocpr$roc <- rocperf
  rocpr$roc.rgraph <- rg.roc
  rocpr$pr <- prperf
  rocpr$pr.rgraph <- rg.pr
  if (roc == TRUE && pr == FALSE) {
    rocpr <- rocpr[-c(4, 5, 8, 9)]
    attributes(rocpr)$label <- "Receiver-operating characteristics"
    class(rocpr) <- "roc"
  } else if (roc == FALSE && pr == TRUE) {
    rocpr <- rocpr[-c(2, 3, 6, 7)]
    attributes(rocpr)$label <- "Precision-recall curve"
    class(rocpr) <- "pr"
  } else {
    attributes(rocpr)$label <- "ROC and PR curve"
    class(rocpr) <- "rocpr"
  }
  return(rocpr)
}