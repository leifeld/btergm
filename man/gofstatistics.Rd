\name{gofstatistics}
\alias{gofstatistics}
\alias{gof-statistics}
\alias{statistics}
\alias{gofterms}
\alias{gof-terms}
\alias{comemb}
\alias{dsp}
\alias{esp}
\alias{nsp}
\alias{deg}
\alias{ideg}
\alias{odeg}
\alias{b1deg}
\alias{b2deg}
\alias{istar}
\alias{kstar}
\alias{ostar}
\alias{b1star}
\alias{b2star}
\alias{kcycle}
\alias{geodesic}
\alias{triad.directed}
\alias{triad.undirected}
\alias{pr}
\alias{roc}
\alias{rocpr}
\alias{rocprgof}
\alias{edgebetweenness.modularity}
\alias{edgebetweenness.pr}
\alias{edgebetweenness.roc}
\alias{fastgreedy.modularity}
\alias{fastgreedy.pr}
\alias{fastgreedy.roc}
\alias{maxmod.modularity}
\alias{maxmod.pr}
\alias{maxmod.roc}
\alias{spinglass.modularity}
\alias{spinglass.pr}
\alias{spinglass.roc}
\alias{walktrap.modularity}
\alias{walktrap.pr}
\alias{walktrap.roc}
\title{Statistics for goodness-of-fit assessment of network models}
\description{
Statistics for goodness-of-fit assessment of network models.
}
\details{
These functions can be plugged into the \code{statistics} argument of the \code{gof} methods in order to compare observed with simulated networks (see the \link{gof-methods} help page). There are three types of statistics:

(1) Univariate statistics, which aggregate a network into a single quantity. For example, modularity measures or density. The distribution of statistics can be displayed using histograms, density plots, and median bars. Univariate statistics take a sparse matrix (\code{mat}) as an argument and return a single numeric value that summarize a network matrix.

(2) Multivariate statistics, which aggregate a network into a vector of quantities. For example, the distribution of geodesic distances, edgewise shared partners, or indegree. These statistics typically have multiple values, e.g., esp(1), esp(2), esp(3) etc. The results can be displayed using multiple boxplots for simulated networks and a black curve for the observed network(s). Multivariate statistics take a sparse matrix (\code{mat}) as an argument and return a vector of numeric values that summarize a network matrix.

(3) Tie prediction statistics, which predict dyad states the observed network(s) by the dyad states in the simulated networks. For example, receiver operating characteristics (ROC) or precision-recall curves (PR) of simulated networks based on the model, or ROC or PR predictions of community co-membership matrices of the simulated vs. the observed network(s). Tie prediction statistics take a list of simulated sparse network matrices and another list of observed sparse network matrices (possibly containing only a single sparse matrix) as arguments and return a \code{rocpr}, \code{roc}, or \code{pr} object (as created by the respective functions \link{rocpr}, \link{rocprgof}, \link{roc}, and \link{pr}).

Users can create their own statistics for use with the code{gof} methods. To do so, one needs to write a function that accepts and returns the respective objects described in the enumeration above. It is advisable to look at the definitions of some of the existing functions to add custom functions. It is also possible to add an attribute called \code{label} to the return object, which describes what is being returned by the function. This label will be used as a descriptive label in the plot and for verbose output during computations. The examples section contains an example of a custom user statistic.

To aid the development of custom statistics, several helper functions are available: The \code{roc}, \code{pr}, and \code{rocpr} functions accept lists of simulated and observed sparse network matrices and compute ROC and precision recall curves as well as the area under the curve that can be used as network statistics. These functions are used internally for a number of functions related to community structure, where the community structure in the simulated networks is compared to the community structure in the observed network(s) by means of tie prediction. The \code{rocprgof} function provides the same functionality as the \code{rocpr} function, but it has an additional argument for controlling imputation of the first PR value. Another helper function is \code{comemb}, which accepts a vector of community memberships and converts it to a co-membership matrix. This function is also used internally by statistics like \code{walktrap.roc} and others.
}
\usage{
b1deg(mat)

b1star(mat)

b2deg(mat)

b2star(mat)

comemb(vec)

deg(mat)

dsp(mat)

edgebetweenness.modularity(mat)

edgebetweenness.pr(sim, obs)

edgebetweenness.roc(sim, obs)

esp(mat)

fastgreedy.modularity(mat)

fastgreedy.pr(sim, obs)

fastgreedy.roc(sim, obs)

geodesic(mat)

ideg(mat)

istar(mat)

kcycle(mat)

kstar(mat)

maxmod.modularity(mat)

maxmod.pr(sim, obs)

maxmod.roc(sim, obs)

nsp(mat)

odeg(mat)

ostar(mat)

pr(sim, obs)

roc(sim, obs)

rocpr(sim, obs)

rocprgof(sim, obs, pr.impute = "poly4")

spinglass.modularity(mat)

spinglass.pr(sim, obs)

spinglass.roc(sim, obs)

triad.directed(mat)

triad.undirected(mat)

walktrap.modularity(mat)

walktrap.pr(sim, obs)

walktrap.roc(sim, obs)
}
\section{Network statistics}{
The following built-in functions can be handed over to the \code{statistics} argument. See the usage section for their respective arguments.

(1) Univariate statistics:
\describe{

\item{\code{walktrap.modularity(mat)}}{
Modularity distribution as computed by the Walktrap algorithm.
}

\item{\code{fastgreedy.modularity(mat)}}{
Modularity distribution as computed by the fast and greedy algorithm. Only sensible with undirected networks.
}

\item{\code{maxmod.modularity(mat)}}{
Optimal modularity distribution.
}

\item{\code{edgebetweenness.modularity(mat)}}{
Modularity distribution as computed by the Girvan-Newman edge betweenness community detection method.
}

\item{\code{spinglass.modularity(mat)}}{
Modularity distribution as computed by the Spinglass algorithm.
}

}

(2) Multivariate statistics:
\describe{
\item{\code{dsp}}{
Dyad-wise shared partner distribution.
}

\item{\code{esp(mat)}}{
Edge-wise shared partner distribution.
}

\item{\code{nsp(mat)}}{
Non-edge-wise shared partner distribution.
}

\item{\code{deg(mat)}}{
Degree distribution (for undirected networks).
}

\item{\code{ideg(mat)}}{
Indegree distribution (for directed networks).
}

\item{\code{odeg(mat)}}{
Outdegree distribution (for directed networks).
}

\item{\code{b1deg(mat)}}{
Degree distribution (for the first mode in a two-mode network).
}

\item{\code{b2deg(mat)}}{
Degree distribution (for the second mode in a two-mode network).
}

\item{\code{kstar(mat)}}{
k-star distribution (for undirected networks).
}

\item{\code{istar(mat)}}{
in-star distribution (for directed networks).
}

\item{\code{ostar(mat)}}{
out-star distribution (for directed networks).
}

\item{\code{b1star(mat)}}{
k-star distribution (for the first mode in a two-mode network).
}

\item{\code{b2star(mat)}}{
k-star distribution (for the second mode in a two-mode network).
}

\item{\code{kcycle(mat)}}{
k-cycle distribution (for undirected networks).
}

\item{\code{geodesic(mat)}}{
Geodesic distance (or shortest path) distribution.
}

\item{\code{triad.directed(mat)}}{
Triad census (directed networks).
}

\item{\code{triad.undirected(mat)}}{
Triad census (undirected networks).
}

}

(3) Tie prediction statistics:
\describe{
\item{\code{walktrap.roc(sim, obs)}}{
Receiver-operating characteristics of predicting the community structure in the observed network(s) by the community structure in the simulated networks, as computed by the Walktrap algorithm.
}

\item{\code{walktrap.pr(sim, obs)}}{
Precision-recall curve for predicting the community structure in the observed network(s) by the community structure in the simulated networks, as computed by the Walktrap algorithm.
}

\item{\code{fastgreedy.roc(sim, obs)}}{
Receiver-operating characteristics of predicting the community structure in the observed network(s) by the community structure in the simulated networks, as computed by the fast and greedy algorithm. Only sensible with undirected networks.
}

\item{\code{fastgreedy.pr(sim, obs)}}{
Precision-recall curve for predicting the community structure in the observed network(s) by the community structure in the simulated networks, as computed by the fast and greedy algorithm. Only sensible with undirected networks.
}

\item{\code{maxmod.roc(sim, obs)}}{
Receiver-operating characteristics of predicting the community structure in the observed network(s) by the community structure in the simulated networks, as computed by the modularity maximization algorithm.
}

\item{\code{maxmod.pr(sim, obs)}}{
Precision-recall curve for predicting the community structure in the observed network(s) by the community structure in the simulated networks, as computed by the modularity maximization algorithm.
}

\item{\code{edgebetweenness.roc(sim, obs)}}{
Receiver-operating characteristics of predicting the community structure in the observed network(s) by the community structure in the simulated networks, as computed by the Girvan-Newman edge betweenness community detection method.
}

\item{\code{edgebetweenness.pr(sim, obs)}}{
Precision-recall curve for predicting the community structure in the observed network(s) by the community structure in the simulated networks, as computed by the Girvan-Newman edge betweenness community detection method.
}

\item{\code{spinglass.roc(sim, obs)}}{
Receiver-operating characteristics of predicting the community structure in the observed network(s) by the community structure in the simulated networks, as computed by the Spinglass algorithm.
}

\item{\code{spinglass.pr(sim, obs)}}{
Precision-recall curve for predicting the community structure in the observed network(s) by the community structure in the simulated networks, as computed by the Spinglass algorithm.
}

\item{\code{roc(sim, obs)}}{
Receiver-operating characteristics. Prediction of the dyad states of the observed network(s) by the dyad states of the simulated networks.
}

\item{\code{pr(sim, obs)}}{
Precision-recall curve. Prediction of the dyad states of the observed network(s) by the dyad states of the simulated networks.
}

\item{\code{rocpr(sim, obs)}}{
Both receiver-operating characteristics and precision-recall curve. Prediction of the dyad states of the observed network(s) by the dyad states of the simulated networks.
}
}
}
\arguments{
\item{vec}{ A vector of community memberships in order to create a community co-membership matrix. }
\item{mat}{ A sparse network matrix as created by the \code{Matrix} function in the \pkg{Matrix} package. }
\item{sim}{ A list of simulated networks. Each element in the list should be a sparse matrix as created by the \code{Matrix} function in the \pkg{Matrix} package. }
\item{obs}{ A list of observed (= target) networks. Each element in the list should be a sparse matrix as created by the \code{Matrix} function in the \pkg{Matrix} package. }
\item{pr.impute}{ In some cases, the first precision value of the precision-recall curve is undefined. The \code{pr.impute} argument serves to impute this missing value to ensure that the AUC-PR value is not severely biased. Possible values are \code{"no"} for no imputation, \code{"one"} for using a value of \code{1.0}, \code{"second"} for using the next (= adjacent) precision value, \code{"poly1"} for fitting a straight line through the remaining curve to predict the first value, \code{"poly2"} for fitting a second-order polynomial curve etc. until \code{"poly9"}. Warning: this is a pragmatic solution. Please double-check whether the imputation makes sense. This can be checked by plotting the resulting object and using the \code{pr.poly} argument to plot the predicted curve on top of the actual PR curve. }
}
\seealso{
\link{btergm-package} \link{gof} \link{gof-methods}
}
\examples{
# To see how these statistics are used, look at the examples section of 
# ?"gof-methods". The following example illustrates how custom 
# statistics can be created. Suppose one is interested in the density 
# of a network. Then a univariate statistic can be created as follows.

dens <- function(mat) {             # univariate: one argument
  mat <- as.matrix(mat)             # sparse matrix -> normal matrix
  d <- sna::gden(mat)               # compute the actual statistic
  attributes(d)$label <- "Density"  # add a descriptive label
  return(d)                         # return the statistic
}

# Now the statistic can be used in the statistics argument of one of 
# the gof methods.

# For illustrative purposes, let us consider an existing statistic, the 
# indegree distribution, a multivariate statistic. It also accepts a 
# single argument. Note that the sparse matrix is converted to a 
# normal matrix object when it is used. First, statnet's summary 
# method is used to compute the statistic. Names are attached to the 
# resulting vector for the different indegree values. Then the vector 
# is returned.

ideg <- function(mat) {
  d <- summary(mat ~ idegree(0:(nrow(mat) - 1)))
  names(d) <- 0:(length(d) - 1)
  attributes(d)$label <- "Indegree"
  return(d)
}

# See the gofstatistics.R file in the package for more complex examples.
}
\author{
Philip Leifeld (\url{http://www.philipleifeld.com})
}
\keyword{statistics}
\keyword{gof}
